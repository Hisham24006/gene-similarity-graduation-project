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


# -----------------------------------------------------------------------
# Edit Distance Similarity (Levenshtein)
# -----------------------------------------------------------------------

def edit_distance(seq_a, seq_b):
    """
    Computes the Levenshtein edit distance between two sequences
    using dynamic programming.

    Each insertion, deletion, or substitution costs 1.
    Returns the minimum number of edits to turn seq_a into seq_b.

    Example:
        edit_distance("MEEPQ", "MKKPQ") -> 2
        (substitute E->K at position 2, substitute E->K at position 3)
    """
    len_a, len_b = len(seq_a), len(seq_b)

    # Build a (len_a+1) x (len_b+1) grid
    # dp[i][j] = edit distance between seq_a[:i] and seq_b[:j]
    dp = list(range(len_b + 1))

    for i in range(1, len_a + 1):
        prev = dp[0]
        dp[0] = i
        for j in range(1, len_b + 1):
            temp = dp[j]
            if seq_a[i - 1] == seq_b[j - 1]:
                dp[j] = prev              # no edit needed
            else:
                dp[j] = 1 + min(
                    prev,                 # substitution
                    dp[j],                # deletion
                    dp[j - 1]            # insertion
                )
            prev = temp

    return dp[len_b]


def edit_distance_similarity(seq_a, seq_b):
    """
    Normalizes edit distance into a 0.0 to 1.0 similarity score.

    similarity = 1 - (edit_distance / max_possible_distance)

    Max possible distance = length of the longer sequence
    (replacing every character of the shorter with the longer).

    Returns 1.0 for identical sequences, 0.0 for completely different.

    Note: For very long sequences (>1000 AA) this can be slow.
    Consider using kmer_similarity for a faster approximation.
    """
    if not seq_a and not seq_b:
        return 1.0
    if not seq_a or not seq_b:
        return 0.0

    distance = edit_distance(seq_a, seq_b)
    max_distance = max(len(seq_a), len(seq_b))
    return 1.0 - (distance / max_distance)
