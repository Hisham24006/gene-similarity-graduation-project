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

TOP_N = 80
K = 3
TEST_QUERIES = 10

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

    print("=== Track A Multi-Query Filter Benchmark ===")
    print(f"Database size: {len(all_seqs)} protein isoforms")
    print(f"Candidate limit: {TOP_N}")
    print(f"Queries tested: {TEST_QUERIES}")

    # Pick queries spread across the database
    step = max(1, len(all_seqs) // TEST_QUERIES)
    test_seqs = all_seqs[::step][:TEST_QUERIES]

    speedups = []
    recalls = []

    for query_number, (symbol, isoform_id, query_seq) in enumerate(test_seqs, 1):
        print(f"\n--- Query {query_number}: {symbol} | {isoform_id} ---")

        # OLD METHOD:
        # Refine every protein and save BLOSUM scores for quality testing
        full_ranked = []

        start = time.perf_counter()

        for candidate_symbol, candidate_isoform, local_seq in all_seqs:
            blosum_score = blosum_similarity(query_seq, local_seq)
            kmer_similarity(query_seq, local_seq, k=K)
            edit_distance_similarity(query_seq, local_seq)

            full_ranked.append(
                (candidate_symbol, candidate_isoform, blosum_score)
            )

        baseline_time = time.perf_counter() - start

        # NEW METHOD:
        # Use k-mer filtering to select candidates
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

        speedup = baseline_time / new_total_time
        speedups.append(speedup)

        # Quality check:
        # How many of the BLOSUM full-search top 10 survived filtering?
        full_ranked.sort(key=lambda x: x[2], reverse=True)
        full_top_10 = full_ranked[:10]

        candidate_ids = {
            (candidate_symbol, candidate_isoform)
            for candidate_symbol, candidate_isoform, _ in candidates
        }

        preserved = sum(
            1
            for candidate_symbol, candidate_isoform, _ in full_top_10
            if (candidate_symbol, candidate_isoform) in candidate_ids
        )

        recall = (preserved / 10) * 100
        recalls.append(recall)

        print(f"Old full refinement: {baseline_time:.3f} s")
        print(f"New filter + refinement: {new_total_time:.3f} s")
        print(f"Speedup: {speedup:.2f}x")
        print(f"Recall@10: {recall:.2f}%")
        

    print("\n=== Overall Results ===")
    print(f"Average speedup: {sum(speedups) / len(speedups):.2f}x")
    print(f"Average Recall@10: {sum(recalls) / len(recalls):.2f}%")
    print(f"Lowest Recall@10: {min(recalls):.2f}%")

    db.close()

if __name__ == "__main__":
    main()
