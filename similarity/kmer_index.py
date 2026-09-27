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


def find_candidates_normalized(query_sequence, index, sequences, k=3, top_n=10):
    """
    Rank candidates using normalized k-mer overlap.
    """
    scores = {}
    query_kmers = set(get_kmers(query_sequence, k))

    for kmer in query_kmers:
        for sequence_id in index.get(kmer, set()):
            scores[sequence_id] = scores.get(sequence_id, 0) + 1

    normalized_scores = {}

    for sequence_id, shared_count in scores.items():
        candidate_kmers = set(get_kmers(sequences[sequence_id], k))

        union_size = len(query_kmers | candidate_kmers)

        normalized_scores[sequence_id] = (
            shared_count / union_size if union_size > 0 else 0
        )

    ranked = sorted(
        normalized_scores.items(),
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


def get_candidate_sequences(query_sequence, all_sequences, k=3, top_n=20):
    """
    Return the top candidate sequences using the k-mer index.
    Output format matches the database:
    (symbol, isoform_id, sequence)
    """
    index = build_index_from_database(all_sequences, k)
    candidates = find_candidates(query_sequence, index, k, top_n)

    sequence_lookup = {
        (symbol, isoform_id): sequence
        for symbol, isoform_id, sequence in all_sequences
    }

    return [
        (symbol, isoform_id, sequence_lookup[(symbol, isoform_id)])
        for (symbol, isoform_id), score in candidates
    ]

