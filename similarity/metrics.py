"""
metrics.py
----------
Similarity metrics for gene/protein sequence comparison.

Available metrics:
    - kmer_similarity(seq_a, seq_b, k)  : k-mer based similarity (fast, alignment-free)

Planned metrics (to be added):
    - edit_distance_similarity(seq_a, seq_b)
    - motif_similarity(seq_a, seq_b)
"""


# -----------------------------------------------------------------------
# K-mer Similarity
# -----------------------------------------------------------------------

def get_kmers(sequence, k):
    """
    Splits a sequence into all overlapping k-mers.
    Returns a set of unique k-mers.

    Example (k=3):
        "MEEPQ" -> {"MEE", "EEP", "EPQ"}
    """
    if len(sequence) < k:
        return set()
    return set(sequence[i:i+k] for i in range(len(sequence) - k + 1))


def kmer_similarity(seq_a, seq_b, k=3):
    """
    Computes the Jaccard similarity between two sequences based on their k-mers.

    Jaccard similarity = |shared k-mers| / |total unique k-mers|

    Returns a float between 0.0 (nothing shared) and 1.0 (identical k-mer sets).

    Args:
        seq_a : first sequence string
        seq_b : second sequence string
        k     : k-mer size (default 3 for proteins, use 7-11 for DNA)
    """
    kmers_a = get_kmers(seq_a, k)
    kmers_b = get_kmers(seq_b, k)

    if not kmers_a and not kmers_b:
        return 0.0

    intersection = kmers_a & kmers_b  # shared k-mers
    union        = kmers_a | kmers_b  # all unique k-mers

    return len(intersection) / len(union)
