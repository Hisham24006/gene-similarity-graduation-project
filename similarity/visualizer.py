"""
visualizer.py
-------------
Generates visualizations for gene similarity search results.

Functions:
    plot_bar_chart(results, query_gene, output_dir)
        Top 10 most similar genes — all 3 metrics side by side.

    plot_heatmap(db, metric_fn, seq_type, output_dir)
        All genes vs all genes similarity heatmap (one metric at a time).
"""

import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import seaborn as sns


# -----------------------------------------------------------------------
# Bar Chart — Top 10 results for a query gene
# -----------------------------------------------------------------------

def plot_bar_chart(results, query_gene, output_dir="results"):
    """
    Plots a grouped bar chart of the top 10 most similar genes
    showing BLOSUM62, K-mer, and Edit distance scores side by side.

    Args:
        results    : list of (symbol, isoform_id, blosum, kmer, edit, avg)
                     sorted by average score descending
        query_gene : name of the query gene (used in title and filename)
        output_dir : folder to save the chart image
    """
    os.makedirs(output_dir, exist_ok=True)

    # Take top 10, skip self-matches (score == 1.0 on all metrics)
    top = results[:10]

    labels   = [f"{s}\n{i[:12]}" for s, i, b, k, e, a in top]
    blosum   = [b for s, i, b, k, e, a in top]
    kmer     = [k for s, i, b, k, e, a in top]
    edit     = [e for s, i, b, k, e, a in top]

    x      = np.arange(len(labels))
    width  = 0.25

    fig, ax = plt.subplots(figsize=(14, 6))

    bars1 = ax.bar(x - width, blosum, width, label="BLOSUM62",     color="#4C72B0")
    bars2 = ax.bar(x,         kmer,   width, label="K-mer (k=3)",  color="#55A868")
    bars3 = ax.bar(x + width, edit,   width, label="Edit distance", color="#C44E52")

    ax.set_title(f"Top 10 similar genes to {query_gene}", fontsize=14, pad=15)
    ax.set_xlabel("Gene | Isoform", fontsize=11)
    ax.set_ylabel("Similarity score (0–1)", fontsize=11)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=8)
    ax.set_ylim(0, 1.1)
    ax.yaxis.set_major_formatter(ticker.FormatStrFormatter("%.2f"))
    ax.legend(fontsize=10)
    ax.grid(axis="y", linestyle="--", alpha=0.4)

    fig.tight_layout()

    path = os.path.join(output_dir, f"bar_{query_gene}.png")
    fig.savefig(path, dpi=150)
    plt.close(fig)
    print(f"  Bar chart saved → {path}")
    return path


# -----------------------------------------------------------------------
# Heatmap — All genes vs all genes
# -----------------------------------------------------------------------

def plot_heatmap(all_seqs, metric_fn, metric_name, query_gene, output_dir="results"):
    """
    Computes and plots a similarity heatmap of all genes in the database
    against each other using the provided metric function.

    Uses the first isoform per gene to keep the matrix manageable.

    Args:
        all_seqs    : list of (symbol, isoform_id, sequence) from the database
        metric_fn   : similarity function (seq_a, seq_b) -> float
        metric_name : label for the metric (used in title and filename)
        query_gene  : used in filename to group outputs by query
        output_dir  : folder to save the heatmap image
    """
    os.makedirs(output_dir, exist_ok=True)

    # Use only the first isoform per gene symbol to keep matrix clean
    seen = {}
    for symbol, isoform_id, seq in all_seqs:
        if symbol not in seen:
            seen[symbol] = seq

    genes = sorted(seen.keys())
    n     = len(genes)

    print(f"  Computing {metric_name} heatmap ({n}x{n} matrix)...")

    # Build similarity matrix
    matrix = np.zeros((n, n))
    for i, gene_a in enumerate(genes):
        for j, gene_b in enumerate(genes):
            if i == j:
                matrix[i][j] = 1.0
            elif j > i:
                score = metric_fn(seen[gene_a], seen[gene_b])
                matrix[i][j] = score
                matrix[j][i] = score  # symmetric

    fig, ax = plt.subplots(figsize=(16, 14))

    sns.heatmap(
        matrix,
        xticklabels=genes,
        yticklabels=genes,
        cmap="YlOrRd",
        vmin=0,
        vmax=1,
        ax=ax,
        square=True,
        linewidths=0.3,
        linecolor="white",
        cbar_kws={"label": "Similarity (0–1)", "shrink": 0.8}
    )

    ax.set_title(f"Gene similarity heatmap — {metric_name}", fontsize=14, pad=15)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=90, fontsize=7)
    ax.set_yticklabels(ax.get_yticklabels(), rotation=0,  fontsize=7)

    fig.tight_layout()

    safe_name = metric_name.lower().replace(" ", "_").replace("(", "").replace(")", "")
    path = os.path.join(output_dir, f"heatmap_{safe_name}.png")
    fig.savefig(path, dpi=150)
    plt.close(fig)
    print(f"  Heatmap saved     → {path}")
    return path
