import threading


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
