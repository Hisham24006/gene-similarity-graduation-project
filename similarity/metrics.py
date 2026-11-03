"""
metrics.py
----------
Similarity metrics for gene and protein sequence comparison.

Available metrics:
    - blosum_similarity(seq_a, seq_b)        : local alignment normalized by max self-score (proteins)
    - dna_similarity(seq_a, seq_b)           : local alignment for nucleotide sequences (DNA/RNA)
    - kmer_similarity(seq_a, seq_b, k)       : Jaccard similarity of k-mers (fast, alignment-free)
    - edit_distance_similarity(seq_a, seq_b) : normalized Levenshtein distance
"""

from Bio.Align import PairwiseAligner, substitution_matrices

# -----------------------------------------------------------------------
# Shared Pairwise Aligners
# -----------------------------------------------------------------------

# Protein Aligner (BLOSUM62)
protein_aligner = PairwiseAligner()
protein_aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
protein_aligner.mode = 'local'
protein_aligner.open_gap_score = -10
protein_aligner.extend_gap_score = -0.5

# DNA / Nucleotide Aligner
dna_aligner = PairwiseAligner()
dna_aligner.mode = 'local'
dna_aligner.match_score = 2
dna_aligner.mismatch_score = -3
dna_aligner.open_gap_score = -5
dna_aligner.extend_gap_score = -2


def blosum_similarity(seq_a, seq_b):
    """
    Normalized BLOSUM62 local alignment similarity score (0.0 to 1.0).
    Divides the local alignment score by the self-alignment score of seq_a.
    """
    score = protein_aligner.score(seq_a, seq_b)
    mx = protein_aligner.score(seq_a, seq_a)
    return max(0.0, score / mx) if mx > 0 else 0.0


def dna_similarity(seq_a, seq_b):
    """
    Normalized DNA local alignment score (0.0 to 1.0).
    Divides the local alignment score by the self-alignment score of seq_a.
    """
    score = dna_aligner.score(seq_a, seq_b)
    mx = dna_aligner.score(seq_a, seq_a)
    return max(0.0, score / mx) if mx > 0 else 0.0


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
    return set(sequence[i:i + k] for i in range(len(sequence) - k + 1))


def kmer_similarity(seq_a, seq_b, k=3):
    """
    Computes the Jaccard similarity between two sequences based on their k-mers.

    Jaccard similarity = |shared k-mers| / |total unique k-mers|

    Returns a float between 0.0 (nothing shared) and 1.0 (identical k-mer sets).

    Args:
        seq_a : first sequence string
        seq_b : second sequence string
        k     : k-mer size (default 3 for proteins, 7-11 for DNA)
    """
    kmers_a = get_kmers(seq_a, k)
    kmers_b = get_kmers(seq_b, k)

    if not kmers_a and not kmers_b:
        return 0.0

    intersection = kmers_a & kmers_b
    union = kmers_a | kmers_b

    return len(intersection) / len(union)


# -----------------------------------------------------------------------
# Edit Distance Similarity (Levenshtein)
# -----------------------------------------------------------------------

def edit_distance(seq_a, seq_b):
    """
    Computes the Levenshtein edit distance between two sequences
    using dynamic programming.

    Each insertion, deletion, or substitution costs 1.
    """
    len_a, len_b = len(seq_a), len(seq_b)

    # 1D DP array of size (len_b + 1)
    dp = list(range(len_b + 1))

    for i in range(1, len_a + 1):
        prev = dp[0]
        dp[0] = i
        for j in range(1, len_b + 1):
            temp = dp[j]
            if seq_a[i - 1] == seq_b[j - 1]:
                dp[j] = prev
            else:
                dp[j] = 1 + min(prev, dp[j], dp[j - 1])
            prev = temp

    return dp[len_b]


def edit_distance_similarity(seq_a, seq_b):
    """
    Normalizes edit distance into a 0.0 to 1.0 similarity score.

    similarity = 1 - (edit_distance / max_possible_distance)
    Max possible distance = length of the longer sequence.
    """
    if not seq_a and not seq_b:
        return 1.0
    if not seq_a or not seq_b:
        return 0.0

    distance = edit_distance(seq_a, seq_b)
    max_distance = max(len(seq_a), len(seq_b))
    return 1.0 - (distance / max_distance)
