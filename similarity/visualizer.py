"""
visualizer.py
-------------
Generates visualizations for gene similarity search results.

Functions:
    plot_bar_chart(results, query_gene, output_dir)
        Top 10 most similar sequences, one bar per metric that was run.

    plot_heatmap(db, metric_fn, seq_type, output_dir)
        All genes vs all genes similarity heatmap (one metric at a time).
"""

import os
import re
import numpy as np
import matplotlib
matplotlib.use("Agg")  # file-only backend: charts are drawn in a web-app background thread
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import seaborn as sns


# -----------------------------------------------------------------------
# Bar Chart — Top 10 results for a query gene
# -----------------------------------------------------------------------

BAR_METRICS = [
    ("blosum", "Alignment (BLOSUM62 / DNA)", "#4C72B0"),
    ("kmer",   "K-mer (k=3)",                "#55A868"),
    ("edit",   "Edit distance",              "#C44E52"),
    ("motif",  "Motif (JASPAR)",             "#DD8452"),
]


def safe_filename(name):
    """Replaces characters that are not allowed in file names (e.g. '|')."""
    return re.sub(r"[^A-Za-z0-9._-]+", "_", str(name)).strip("_") or "query"


def plot_bar_chart(results, query_gene, output_dir="results"):
    """
    Plots a grouped bar chart of the top 10 most similar sequences,
    with one bar per metric that was run.

    Args:
        results    : list of dicts {'gene', 'isoform', 'scores': {metric: score}}
                     in ranked order (the format returned by the web app)
        query_gene : name of the query gene (used in title and filename)
        output_dir : folder to save the chart image
    """
    os.makedirs(output_dir, exist_ok=True)

    # Skip self-matches (identical sequence: 1.0 on every metric), keep top 10
    top = [
        r for r in results
        if not (r["scores"] and all(v >= 1.0 for v in r["scores"].values()))
    ][:10]
    if not top:
        print("  Bar chart skipped: no results besides self-matches.")
        return None

    metrics = [m for m in BAR_METRICS if any(m[0] in r["scores"] for r in top)]
    if not metrics:
        return None

    labels = [f"{r['gene']}\n{r['isoform'][:12]}" for r in top]
    x      = np.arange(len(labels))
    width  = 0.8 / len(metrics)

    fig, ax = plt.subplots(figsize=(14, 6))

    for n, (key, name, color) in enumerate(metrics):
        values = [r["scores"].get(key, 0.0) for r in top]
        offset = (n - (len(metrics) - 1) / 2) * width
        ax.bar(x + offset, values, width, label=name, color=color)

    ax.set_title(f"Top {len(top)} similar sequences to {query_gene}", fontsize=14, pad=15)
    ax.set_xlabel("Gene | Isoform", fontsize=11)
    ax.set_ylabel("Similarity score (0–1)", fontsize=11)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=8)
    ax.set_ylim(0, 1.1)
    ax.yaxis.set_major_formatter(ticker.FormatStrFormatter("%.2f"))
    ax.legend(fontsize=10)
    ax.grid(axis="y", linestyle="--", alpha=0.4)

    fig.tight_layout()

    path = os.path.join(output_dir, f"bar_{safe_filename(query_gene)}.png")
    fig.savefig(path, dpi=150)
    plt.close(fig)
    print(f"  Bar chart saved -> {path}")
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
    print(f"  Heatmap saved     -> {path}")
    return path
