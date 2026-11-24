import threading

import numpy as np


def get_kmers(sequence, k=3):
    """
    Split a protein sequence into overlapping k-mers.
    Example: ABCDE with k=3 -> ABC, BCD, CDE
    """
    sequence = sequence.upper().strip()

    if len(sequence) < k:
        return set()

    return {
        sequence[i:i + k]
        for i in range(len(sequence) - k + 1)
    }


def build_kmer_index(sequences, k=3, sizes=None):
    """
    Build an inverted index:
    k-mer -> set of sequence IDs containing that k-mer.

    If a dict is passed as `sizes`, it is filled with
    sequence ID -> number of unique k-mers in that sequence.
    """
    index = {}

    for sequence_id, sequence in sequences.items():
        kmers = get_kmers(sequence, k)
        if sizes is not None:
            sizes[sequence_id] = len(kmers)

        for kmer in kmers:
            if kmer not in index:
                index[kmer] = set()

            index[kmer].add(sequence_id)

    return index

def find_candidates(query_sequence, index, k=3, top_n=10, sizes=None):
    """
    Find the sequences that share the most k-mers with the query.

    Without `sizes` the ranking is the raw number of shared k-mers. That
    favours very long proteins (titin and similar share some k-mers with
    everything). With `sizes` (sequence ID -> number of unique k-mers) the
    ranking is the Jaccard similarity, shared / (query + sequence - shared),
    which is the same idea as the k-mer metric and treats long and short
    proteins fairly.
    """
    scores = {}
    query_kmers = get_kmers(query_sequence, k)

    for kmer in query_kmers:
        for sequence_id in index.get(kmer, set()):
            scores[sequence_id] = scores.get(sequence_id, 0) + 1

    if sizes is not None:
        n_query = len(query_kmers)
        scores = {
            sequence_id: shared / (n_query + sizes[sequence_id] - shared)
            for sequence_id, shared in scores.items()
        }

    ranked = sorted(
        scores.items(),
        key=lambda item: item[1],
        reverse=True
    )

    return ranked[:top_n]


def build_index_from_database(all_sequences, k=3):
    """
    Build a k-mer index from database rows:
    (gene_symbol, isoform_id, protein_sequence)
    """
    sequences = {}

    for symbol, isoform_id, sequence in all_sequences:
        sequence_id = (symbol, isoform_id)
        sequences[sequence_id] = sequence

    return build_kmer_index(sequences, k)


# The index takes several seconds to build for the whole database, so it is
# built once and reused by every search (rebuilt only if the database changes).
_search_index_cache = {}
_search_index_lock = threading.Lock()


def get_search_index(all_sequences, k=3):
    """
    Returns (index, sizes, lookup) for the given database rows, cached.
    lookup maps (symbol, isoform_id) -> sequence.
    """
    key = (k, len(all_sequences))
    with _search_index_lock:
        cached = _search_index_cache.get(key)
        if cached is None:
            sequences = {
                (symbol, isoform_id): sequence
                for symbol, isoform_id, sequence in all_sequences
            }
            sizes = {}
            index = build_kmer_index(sequences, k, sizes)
            cached = (index, sizes, sequences)
            _search_index_cache.clear()
            _search_index_cache[key] = cached
        return cached


def get_candidate_sequences(query_sequence, all_sequences, k=3, top_n=20, index=None):
    """
    Return the top candidate sequences using the k-mer index.
    Output format matches the database:
    (symbol, isoform_id, sequence)

    When no index is given, the cached index is used and candidates are ranked
    by Jaccard similarity (see find_candidates).
    """
    if index is None:
        if not get_kmers(query_sequence, k):
            # Query shorter than k: it has no k-mers, so the filter cannot
            # rank. Score everything, as get_dna_candidates does.
            return list(all_sequences)
        index, sizes, sequence_lookup = get_search_index(all_sequences, k)
        candidates = find_candidates(query_sequence, index, k, top_n, sizes)
        return [
            (symbol, isoform_id, sequence_lookup[(symbol, isoform_id)])
            for (symbol, isoform_id), score in candidates
        ]

    candidates = find_candidates(query_sequence, index, k, top_n)

    sequence_lookup = {
        (symbol, isoform_id): sequence
        for symbol, isoform_id, sequence in all_sequences
    }

    return [
        (symbol, isoform_id, sequence_lookup[(symbol, isoform_id)])
        for (symbol, isoform_id), score in candidates
    ]


# -----------------------------------------------------------------------
# Nucleotide (DNA) candidate filter
# -----------------------------------------------------------------------
#
# The protein index above stores Python sets, which needs about 480 MB for
# 18,741 proteins. mRNAs are about 7 times longer, so the same structure
# would need several GB. DNA has only 4 letters, so every k-mer (k <= 8)
# fits in a 16-bit integer: each sequence is stored as a sorted array of
# its distinct k-mer codes, and all arrays are kept end to end in one
# NumPy array. Counting the k-mers a query shares with every sequence is
# then a single vectorized lookup.

_DNA_CODE = np.full(256, -1, dtype=np.int8)
for _i, _c in enumerate('ACGT'):
    _DNA_CODE[ord(_c)] = _i
    _DNA_CODE[ord(_c.lower())] = _i


def dna_kmer_codes(sequence, k=7):
    """
    Returns the distinct k-mers of a DNA sequence as sorted integer codes
    (A=0, C=1, G=2, T=3, read as a base-4 number). Windows that contain a
    letter other than A, C, G or T (such as N) are skipped.
    """
    letters = _DNA_CODE[np.frombuffer(sequence.encode('ascii', 'replace'), dtype=np.uint8)]
    n = len(letters) - k + 1
    if n <= 0:
        return np.empty(0, dtype=np.uint32)
    bad = letters < 0
    values = np.where(bad, 0, letters).astype(np.uint32)
    codes = np.zeros(n, dtype=np.uint32)
    for j in range(k):
        codes = codes * 4 + values[j:j + n]
    if bad.any():
        bad_window = np.convolve(bad.astype(np.int32), np.ones(k, dtype=np.int32), 'valid') > 0
        codes = codes[~bad_window]
    return np.unique(codes)


class DnaKmerIndex:
    """Distinct k-mer codes of many DNA sequences, stored end to end."""

    def __init__(self, sequences, k=7):
        if 4 ** k > 65536:
            raise ValueError('DnaKmerIndex supports k <= 8')
        self.k = k
        parts = [dna_kmer_codes(seq, k).astype(np.uint16) for seq in sequences]
        self.sizes = np.array([len(p) for p in parts], dtype=np.int64)
        self.starts = np.zeros(len(parts), dtype=np.int64)
        if len(parts) > 1:
            self.starts[1:] = np.cumsum(self.sizes)[:-1]
        self.codes = np.concatenate(parts) if parts else np.empty(0, dtype=np.uint16)

    def jaccard(self, query_sequence):
        """Jaccard similarity of the query's k-mers with every sequence."""
        query = dna_kmer_codes(query_sequence, self.k)
        shared = np.zeros(len(self.sizes), dtype=np.int64)
        if len(query) == 0 or len(self.codes) == 0:
            return shared.astype(float)
        in_query = np.zeros(4 ** self.k, dtype=bool)
        in_query[query] = True
        hits = in_query[self.codes]
        filled = self.sizes > 0
        shared[filled] = np.add.reduceat(hits, self.starts[filled], dtype=np.int64)
        union = len(query) + self.sizes - shared
        return np.divide(shared, union, out=np.zeros(len(shared)), where=union > 0)


_dna_index_cache = {}
_dna_index_lock = threading.Lock()


def get_dna_index(all_sequences, k=7):
    """Returns the DnaKmerIndex for these database rows, built once and cached."""
    key = (k, len(all_sequences))
    with _dna_index_lock:
        cached = _dna_index_cache.get(key)
        if cached is None:
            cached = DnaKmerIndex([seq for _, _, seq in all_sequences], k)
            _dna_index_cache.clear()
            _dna_index_cache[key] = cached
        return cached


def get_dna_candidates(query_sequence, all_sequences, k=7, top_n=500):
    """
    Returns the top_n rows (symbol, isoform_id, sequence) whose k-mer sets
    are most similar to the query's (Jaccard), best first.
    """
    if len(all_sequences) <= top_n:
        return list(all_sequences)
    if len(dna_kmer_codes(query_sequence, k)) == 0:
        # Too short (or all N) to have any k-mer: the filter cannot rank,
        # so score everything. Such queries are cheap to align.
        return list(all_sequences)
    scores = get_dna_index(all_sequences, k).jaccard(query_sequence)
    top = np.argpartition(-scores, top_n - 1)[:top_n]
    top = top[np.argsort(-scores[top], kind='stable')]
    return [all_sequences[i] for i in top]
