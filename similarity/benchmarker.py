"""
benchmarker.py
--------------
Benchmarks all similarity metrics in terms of:
    - Speed    : how long each metric takes to run
    - Agreement: how well each metric agrees with BLOSUM62 (Pearson correlation)

Output:
    - Printed benchmark table
    - results/benchmark_speed.png     : bar chart of average time per comparison
    - results/benchmark_agreement.png : correlation of each metric vs BLOSUM62

Usage:
    python benchmarker.py
"""

import os
import sys
import time
import ssl
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

if not os.environ.get('PYTHONHTTPSVERIFY', '') and getattr(ssl, '_create_unverified_context', None):
    ssl._create_default_https_context = ssl._create_unverified_context

sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'database'))

from database_manager import GeneDatabase
from metrics import kmer_similarity, edit_distance_similarity
from motif_similarity import load_jaspar_motifs, fetch_promoter, cached_scan_motifs, jaccard_motif_similarity
from Bio.Align import PairwiseAligner, substitution_matrices

# -----------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------

DB_PATH     = os.path.join(os.path.dirname(__file__), '..', 'database', 'gene_vault.db')
RESULTS_DIR = os.path.join(os.path.dirname(__file__), '..', 'results')
KMER_K      = 3
N_PAIRS     = 50  # number of random sequence pairs to benchmark

# -----------------------------------------------------------------------
# BLOSUM62 setup
# -----------------------------------------------------------------------

aligner = PairwiseAligner()
aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
aligner.mode = 'local'
aligner.open_gap_score = -10
aligner.extend_gap_score = -0.5

def blosum_similarity(seq_a, seq_b):
    score = aligner.score(seq_a, seq_b)
    max_score = aligner.score(seq_a, seq_a)
    return max(0.0, score / max_score) if max_score > 0 else 0.0

# -----------------------------------------------------------------------
# Load data
# -----------------------------------------------------------------------

def load_benchmark_data():
    """
    Loads sequence pairs and motif TF sets from the database and cache.
    Returns:
        seq_pairs  : list of (symbol_a, seq_a, symbol_b, seq_b) tuples
        tf_map     : dict {symbol: set of TF names}
    """
    db       = GeneDatabase(db_path=DB_PATH)
    all_seqs = db.get_all_sequences_with_isoform("protein")
    db.close()

    # Use first isoform per gene
    seen = {}
    for symbol, isoform_id, seq in all_seqs:
        if symbol not in seen:
            seen[symbol] = seq

    genes = list(seen.keys())

    # Build sequence pairs (random sample)
    import random
    random.seed(42)
    all_pairs = [(a, b) for i, a in enumerate(genes) for b in genes[i+1:]]
    pairs     = random.sample(all_pairs, min(N_PAIRS, len(all_pairs)))
    seq_pairs = [(a, seen[a], b, seen[b]) for a, b in pairs]

    # Load TF sets from cache
    print("  Loading JASPAR motifs...")
    bio_motifs = load_jaspar_motifs()

    print("  Loading promoter TF sets...")
    tf_map = {}
    for symbol in genes:
        promoter        = fetch_promoter(symbol)
        tf_map[symbol]  = cached_scan_motifs(symbol, promoter, bio_motifs) if promoter else set()

    return seq_pairs, tf_map


# -----------------------------------------------------------------------
# Benchmarking
# -----------------------------------------------------------------------

def benchmark(seq_pairs, tf_map):
    """
    Runs all 4 metrics on each sequence pair, recording scores and times.
    Returns a dict with metric names as keys and lists of (score, time) as values.
    """
    metrics = {
        "BLOSUM62":      {"scores": [], "times": []},
        "K-mer (k=3)":   {"scores": [], "times": []},
        "Edit distance": {"scores": [], "times": []},
        "Motif (JASPAR)":{"scores": [], "times": []},
    }

    total = len(seq_pairs)
    print(f"\n  Benchmarking {total} sequence pairs...\n")

    for idx, (sym_a, seq_a, sym_b, seq_b) in enumerate(seq_pairs, 1):
        print(f"  [{idx}/{total}] {sym_a} vs {sym_b}")

        # BLOSUM62
        t0 = time.perf_counter()
        s  = blosum_similarity(seq_a, seq_b)
        metrics["BLOSUM62"]["scores"].append(s)
        metrics["BLOSUM62"]["times"].append(time.perf_counter() - t0)

        # K-mer
        t0 = time.perf_counter()
        s  = kmer_similarity(seq_a, seq_b, k=KMER_K)
        metrics["K-mer (k=3)"]["scores"].append(s)
        metrics["K-mer (k=3)"]["times"].append(time.perf_counter() - t0)

        # Edit distance
        t0 = time.perf_counter()
        s  = edit_distance_similarity(seq_a, seq_b)
        metrics["Edit distance"]["scores"].append(s)
        metrics["Edit distance"]["times"].append(time.perf_counter() - t0)

        # Motif
        t0 = time.perf_counter()
        s  = jaccard_motif_similarity(
            tf_map.get(sym_a, set()),
            tf_map.get(sym_b, set())
        )
        metrics["Motif (JASPAR)"]["scores"].append(s)
        metrics["Motif (JASPAR)"]["times"].append(time.perf_counter() - t0)

    return metrics


# -----------------------------------------------------------------------
# Results table
# -----------------------------------------------------------------------

def print_table(metrics):
    """Prints a formatted benchmark summary table."""
    blosum_scores = np.array(metrics["BLOSUM62"]["scores"])

    print("\n" + "=" * 72)
    print(f"  {'Metric':<20} {'Avg time (ms)':>15} {'Min (ms)':>10} {'Max (ms)':>10} {'vs BLOSUM62':>12}")
    print("=" * 72)

    for name, data in metrics.items():
        times  = np.array(data["times"]) * 1000  # convert to ms
        scores = np.array(data["scores"])
        avg_t  = np.mean(times)
        min_t  = np.min(times)
        max_t  = np.max(times)

        # Pearson correlation with BLOSUM62
        if name == "BLOSUM62":
            corr = 1.0
        else:
            corr = np.corrcoef(blosum_scores, scores)[0, 1]

        print(f"  {name:<20} {avg_t:>15.2f} {min_t:>10.2f} {max_t:>10.2f} {corr:>12.4f}")

    print("=" * 72)
    print("  vs BLOSUM62 = Pearson correlation (1.0 = perfect agreement)")
    print("=" * 72 + "\n")


# -----------------------------------------------------------------------
# Visualizations
# -----------------------------------------------------------------------

def plot_speed(metrics, output_dir):
    """Bar chart comparing average time per comparison for each metric."""
    os.makedirs(output_dir, exist_ok=True)

    names  = list(metrics.keys())
    times  = [np.mean(metrics[n]["times"]) * 1000 for n in names]
    colors = ["#4C72B0", "#55A868", "#C44E52", "#EF9F27"]

    fig, ax = plt.subplots(figsize=(9, 5))
    bars = ax.bar(names, times, color=colors, width=0.5)

    # Value labels on bars
    for bar, t in zip(bars, times):
        ax.text(bar.get_x() + bar.get_width() / 2,
                bar.get_height() + max(times) * 0.01,
                f"{t:.2f} ms", ha="center", va="bottom", fontsize=9)

    ax.set_title("Average time per comparison by metric", fontsize=13, pad=12)
    ax.set_ylabel("Time (milliseconds)", fontsize=11)
    ax.set_xlabel("Metric", fontsize=11)
    ax.set_yscale("log")  # log scale because BLOSUM62 is much slower
    ax.grid(axis="y", linestyle="--", alpha=0.4)
    ax.yaxis.set_major_formatter(ticker.ScalarFormatter())

    fig.tight_layout()
    path = os.path.join(output_dir, "benchmark_speed.png")
    fig.savefig(path, dpi=150)
    plt.close(fig)
    print(f"  Speed chart saved   → {path}")


def plot_agreement(metrics, output_dir):
    """
    Scatter plots comparing each metric's scores against BLOSUM62.
    One subplot per non-BLOSUM metric.
    """
    os.makedirs(output_dir, exist_ok=True)

    blosum_scores = np.array(metrics["BLOSUM62"]["scores"])
    other_metrics = {k: v for k, v in metrics.items() if k != "BLOSUM62"}
    colors        = ["#55A868", "#C44E52", "#EF9F27"]

    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    fig.suptitle("Metric agreement vs BLOSUM62", fontsize=13, y=1.02)

    for ax, (name, data), color in zip(axes, other_metrics.items(), colors):
        scores = np.array(data["scores"])
        corr   = np.corrcoef(blosum_scores, scores)[0, 1]

        ax.scatter(blosum_scores, scores, color=color, alpha=0.7, edgecolors="white", s=60)

        # Trend line
        m, b   = np.polyfit(blosum_scores, scores, 1)
        x_line = np.linspace(0, 1, 100)
        ax.plot(x_line, m * x_line + b, color="black", linewidth=1, linestyle="--", alpha=0.6)

        ax.set_title(f"{name}\nr = {corr:.4f}", fontsize=11)
        ax.set_xlabel("BLOSUM62 score", fontsize=10)
        ax.set_ylabel(f"{name} score", fontsize=10)
        ax.set_xlim(-0.05, 1.05)
        ax.set_ylim(-0.05, 1.05)
        ax.grid(linestyle="--", alpha=0.3)

    fig.tight_layout()
    path = os.path.join(output_dir, "benchmark_agreement.png")
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Agreement chart saved → {path}")


# -----------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------

if __name__ == "__main__":
    print("\n=== Gene Similarity Metric Benchmarker ===\n")
    print(f"  Pairs to test : {N_PAIRS}")
    print(f"  K-mer k       : {KMER_K}")

    print("\n--- Loading data ---")
    seq_pairs, tf_map = load_benchmark_data()

    print("\n--- Running benchmarks ---")
    metrics = benchmark(seq_pairs, tf_map)

    print("\n--- Results ---")
    print_table(metrics)

    print("--- Generating charts ---")
    plot_speed(metrics, RESULTS_DIR)
    plot_agreement(metrics, RESULTS_DIR)

    print("\nDone.")
