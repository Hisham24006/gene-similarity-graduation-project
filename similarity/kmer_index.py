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


def build_kmer_index(sequences, k=3):
    """
    Build an inverted index:
    k-mer -> set of sequence IDs containing that k-mer.
    """
    index = {}

    for sequence_id, sequence in sequences.items():
        for kmer in get_kmers(sequence, k):
            if kmer not in index:
                index[kmer] = set()

            index[kmer].add(sequence_id)

    return index

def find_candidates(query_sequence, index, k=3, top_n=10):
    """
    Find the sequences that share the most k-mers with the query.
    """
    scores = {}

    for kmer in get_kmers(query_sequence, k):
        for sequence_id in index.get(kmer, set()):
            scores[sequence_id] = scores.get(sequence_id, 0) + 1

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