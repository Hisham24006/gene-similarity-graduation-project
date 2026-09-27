import os
import sys
import time

# Allow access to the database folder
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'database'))

from database_manager import GeneDatabase
from kmer_index import get_candidate_sequences
from metrics import kmer_similarity, edit_distance_similarity
from Bio.Align import PairwiseAligner, substitution_matrices

DB_PATH = os.path.join(
    os.path.dirname(__file__), '..', 'database', 'gene_vault.db'
)

TOP_N = 20
K = 3

aligner = PairwiseAligner()
aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
aligner.mode = "local"
aligner.open_gap_score = -10
aligner.extend_gap_score = -0.5


def blosum_similarity(seq_a, seq_b):
    score = aligner.score(seq_a, seq_b)
    max_score = aligner.score(seq_a, seq_a)
    return max(0.0, score / max_score) if max_score > 0 else 0.0


def main():
    db = GeneDatabase(db_path=DB_PATH)
    all_seqs = db.get_all_sequences_with_isoform("protein")

    # Use the first real protein as the test query
    symbol, isoform_id, query_seq = all_seqs[0]

    print("=== Track A Filter Benchmark ===")
    print(f"Query: {symbol} | {isoform_id}")
    print(f"Database size: {len(all_seqs)} protein isoforms")
    print(f"Candidate limit: {TOP_N}")

    # Baseline: old method refines every protein
    start = time.perf_counter()

    for _, _, local_seq in all_seqs:
        blosum_similarity(query_seq, local_seq)
        kmer_similarity(query_seq, local_seq, k=K)
        edit_distance_similarity(query_seq, local_seq)

    baseline_time = time.perf_counter() - start
    full_search_space = list(all_seqs)

    # New Track A filtering
    start = time.perf_counter()
    candidates = get_candidate_sequences(
        query_seq,
        all_seqs,
        k=K,
        top_n=TOP_N
    )
    filter_time = time.perf_counter() - start

    # Refine only the filtered candidates
    start = time.perf_counter()

    for _, _, local_seq in candidates:
        blosum_similarity(query_seq, local_seq)
        kmer_similarity(query_seq, local_seq, k=K)
        edit_distance_similarity(query_seq, local_seq)

    refinement_time = time.perf_counter() - start
    new_total_time = filter_time + refinement_time

    reduction = (
        1 - (len(candidates) / len(full_search_space))
    ) * 100

    print("\n--- Results ---")
    print(f"Old refinement candidates: {len(full_search_space)}")
    print(f"New refinement candidates: {len(candidates)}")
    print(f"Search-space reduction: {reduction:.2f}%")
    print(f"K-mer filtering time: {filter_time * 1000:.3f} ms")
    print(f"Old full refinement time: {baseline_time:.3f} s")
    print(f"New refinement time: {refinement_time:.3f} s")
    print(f"New total time (filter + refine): {new_total_time:.3f} s")
    print(f"Speedup: {baseline_time / new_total_time:.2f}x")

    print("\nTop candidates:")
    for candidate_symbol, candidate_isoform, _ in candidates[:5]:
        print(f"  {candidate_symbol} | {candidate_isoform}")

    db.close()


if __name__ == "__main__":
    main()