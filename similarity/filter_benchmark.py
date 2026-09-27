import os
import sys
import time

# Allow access to the database folder
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'database'))

from database_manager import GeneDatabase
from kmer_index import get_candidate_sequences

DB_PATH = os.path.join(
    os.path.dirname(__file__), '..', 'database', 'gene_vault.db'
)

TOP_N = 20
K = 3


def main():
    db = GeneDatabase(db_path=DB_PATH)
    all_seqs = db.get_all_sequences_with_isoform("protein")

    # Use the first real protein as the test query
    symbol, isoform_id, query_seq = all_seqs[0]

    print("=== Track A Filter Benchmark ===")
    print(f"Query: {symbol} | {isoform_id}")
    print(f"Database size: {len(all_seqs)} protein isoforms")
    print(f"Candidate limit: {TOP_N}")

    # Baseline: number of proteins the old search would refine
    start = time.perf_counter()
    full_search_space = list(all_seqs)
    baseline_time = time.perf_counter() - start

    # New Track A filtering
    start = time.perf_counter()
    candidates = get_candidate_sequences(
        query_seq,
        all_seqs,
        k=K,
        top_n=TOP_N
    )
    filter_time = time.perf_counter() - start

    reduction = (
        1 - (len(candidates) / len(full_search_space))
    ) * 100

    print("\n--- Results ---")
    print(f"Old refinement candidates: {len(full_search_space)}")
    print(f"New refinement candidates: {len(candidates)}")
    print(f"Search-space reduction: {reduction:.2f}%")
    print(f"K-mer filtering time: {filter_time * 1000:.3f} ms")

    print("\nTop candidates:")
    for candidate_symbol, candidate_isoform, _ in candidates[:5]:
        print(f"  {candidate_symbol} | {candidate_isoform}")

    db.close()


if __name__ == "__main__":
    main()