# Effective and Efficient Gene Similarity Search

A bioinformatics tool for finding similar genes using multiple similarity metrics, built as a junior graduation project.

Given a gene of interest, the tool searches a local database of 60 human genes and ranks them by similarity using four different metrics — from fast alignment-free methods to biologically rich motif-based comparison.

---

## Project Structure

```
gene-similarity-graduation-project/
├── data/
│   ├── proteins.fasta          # Original protein sequences (legacy)
│   ├── proteins727.fasta       # Original protein sequences (legacy)
│   ├── promoter_cache.json     # Cached Ensembl promoter sequences
│   └── tf_cache.json           # Cached JASPAR TF binding site sets
├── database/
│   ├── build_db.py             # Fetches sequences and builds gene_vault.db
│   ├── database_manager.py     # SQLite read/write interface
│   └── gene_vault.db           # Local gene database (60 genes, 3 sources)
├── similarity/
│   ├── prototype.py            # Main search interface
│   ├── metrics.py              # K-mer and edit distance similarity
│   ├── motif_similarity.py     # JASPAR motif-based similarity
│   ├── visualizer.py           # Bar charts and heatmaps
│   └── benchmarker.py          # Speed and accuracy benchmarks
├── results/                    # Generated visualizations (git-ignored)
├── requirements.txt
└── .gitignore
```

---

## Database

The local database (`gene_vault.db`) contains **60 human genes** across four biological categories:

| Category | Examples |
|---|---|
| Cancer-related | TP53, BRCA1, BRCA2, KRAS, EGFR, ERBB2 |
| Cell signaling | AKT1, MTOR, JAK2, STAT3, MAPK1 |
| Transcription factors | MYC, E2F1, SP1, FOXP3, GATA1 |
| DNA repair | MLH1, RAD51, CHEK1, PALB2, ATM |

For each gene, the database stores:
- **Protein sequences** from NCBI RefSeq and UniProt (up to 5 isoforms each)
- **DNA/mRNA sequences** from NCBI RefSeq and Ensembl (up to 5 isoforms each)
- **Metadata**: description, sequence length, source

---

## Similarity Metrics

| Metric | Method | Speed | Biology-aware |
|---|---|---|---|
| BLOSUM62 | Local pairwise alignment | Slow | Yes — amino acid substitution likelihood |
| K-mer (k=3) | Jaccard similarity of shared 3-mers | Very fast | No |
| Edit distance | Normalized Levenshtein distance | Medium | No |
| Motif (JASPAR) | Shared transcription factor binding sites | Fast (cached) | Yes — functional regulation |

The motif metric is particularly meaningful — two genes sharing the same transcription factor binding sites in their promoters are regulated similarly, even if their protein sequences look very different.

---

## Installation

```bash
git clone https://github.com/Hisham24006/gene-similarity-graduation-project.git
cd gene-similarity-graduation-project
pip install -r requirements.txt
```

---

## Usage

### Run the similarity search

```bash
python similarity/prototype.py
```

You will be prompted to either:
1. **Search by gene name** — the tool fetches the sequence from NCBI automatically
2. **Paste a raw sequence** — supports plain amino acid sequences and FASTA format

The tool then compares your query against all 302 protein isoforms in the database using all four metrics and displays a ranked results table.

### Rebuild or expand the database

```bash
python database/build_db.py
```

Add new gene symbols to the `GENES` list in `build_db.py` before running. Already-stored genes are skipped automatically.

### Run benchmarks

```bash
python similarity/benchmarker.py
```

Tests 50 random gene pairs and outputs a speed and agreement comparison table, plus two chart images in `results/`.

---

## Results

After running the similarity search, three visualization files are saved to `results/`:

- **`bar_<GENE>.png`** — top 10 most similar genes with all 4 metric scores side by side
- **`heatmap_k-mer_k=3.png`** — 60×60 gene similarity heatmap using k-mer
- **`heatmap_edit_distance.png`** — 60×60 gene similarity heatmap using edit distance
- **`heatmap_motif_jaspar.png`** — 60×60 gene similarity heatmap using motif similarity
- **`benchmark_speed.png`** — average time per comparison for each metric (log scale)
- **`benchmark_agreement.png`** — scatter plots showing correlation of each metric with BLOSUM62

---

## Dependencies

```
biopython==1.86
requests==2.32.5
numpy==2.4.3
matplotlib==3.10.8
seaborn==0.13.2
pyjaspar
```

Install with:
```bash
pip install -r requirements.txt
```

---

## Data Sources

- **NCBI** — protein and mRNA sequences via Entrez API
- **UniProt** — reviewed protein sequences via REST API
- **Ensembl** — DNA sequences and promoter regions via REST API
- **JASPAR 2024** — transcription factor binding profiles via pyJASPAR

---

## Authors

Hisham Alsaadi — University of Sharjah
