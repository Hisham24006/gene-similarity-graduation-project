# Effective and Efficient Gene Similarity Search

A bioinformatics tool for finding similar genes using multiple similarity metrics, built as a junior graduation project at the University of Sharjah.

Given a gene of interest, the system searches a local database of 60 human genes and ranks them by similarity using four complementary metrics — from fast alignment-free methods to biologically rich motif-based comparison. A Flask web application provides an interactive interface for all features.

---

## Project Structure

```
gene-similarity-graduation-project/
├── data/
│   ├── proteins.fasta            # Legacy protein sequences
│   ├── proteins727.fasta         # Legacy protein sequences
│   ├── promoter_cache.json       # Cached Ensembl promoter sequences (1000 bp per gene)
│   ├── tf_cache.json             # Cached JASPAR TF binding site sets per gene
│   └── recent_searches.json      # Recent search history (auto-generated)
├── database/
│   ├── build_db.py               # Fetches sequences and builds gene_vault.db
│   ├── database_manager.py       # SQLite read/write interface
│   └── gene_vault.db             # Local gene database (60 genes, 3 sources)
├── similarity/
│   ├── prototype.py              # Command-line similarity search interface
│   ├── metrics.py                # K-mer and edit distance implementations
│   ├── motif_similarity.py       # JASPAR motif scanning and TF set similarity
│   ├── visualizer.py             # Bar charts and heatmaps
│   └── benchmarker.py            # Speed and accuracy benchmark script
├── webapp/
│   ├── app.py                    # Flask application (routes + REST API)
│   └── templates/
│       ├── base.html             # Shared layout and Apple-inspired design system
│       ├── search.html           # Similarity search page
│       ├── database.html         # Gene database browser
│       ├── visualize.html        # Visualization charts page
│       └── benchmark.html        # Metric benchmark page
├── results/                      # Generated chart images (git-ignored)
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
- **Metadata**: description, sequence length, source, isoform ID

---

## Similarity Metrics

| Metric | Method | Speed | Biology-aware |
|---|---|---|---|
| BLOSUM62 | Local pairwise alignment | Slow | Yes — amino acid substitution likelihood |
| K-mer (k=3) | Jaccard similarity of shared 3-mers | Very fast | No |
| Edit distance | Normalized Levenshtein distance | Medium | No |
| Motif (JASPAR) | Shared transcription factor binding sites | Fast (cached) | Yes — functional co-regulation |

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

### Run the Web Application (recommended)

```bash
python webapp/app.py
```

Then open your browser to `http://localhost:5000`

The web app has four pages:

- **Search** — Enter a gene name or paste a sequence (protein or DNA). Toggle which metrics to run, set result count, and view ranked results with color-coded score bars. Includes a filter box to search within results and a recent searches panel.
- **Database** — Browse all 60 genes with isoform counts, descriptions, and a sequence viewer that shows full FASTA sequences with a copy button.
- **Visualizations** — View bar charts and heatmaps generated automatically after each search.
- **Benchmark** — Run a speed and accuracy comparison across all four metrics on a configurable number of random gene pairs.

### Run the Command-Line Search

```bash
python similarity/prototype.py
```

### Rebuild or Expand the Database

Add new gene symbols to the `GENES` list in `database/build_db.py`, then run:

```bash
python database/build_db.py
```

Already-stored genes are skipped automatically.

### Run Benchmarks (standalone)

```bash
python similarity/benchmarker.py
```

---

## Dependencies

```
biopython==1.86
requests==2.32.5
numpy==2.4.3
matplotlib==3.10.8
seaborn==0.13.2
pyjaspar
flask
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

## Team

| Name | ID | Specialization | Role |
|---|---|---|---|
| Hisham Alsadi | U23100356 | Biomedical Informatics | Similarity metrics (BLOSUM62, k-mer, edit distance, JASPAR) |
| Yousef Ahmed Ibrahim | U22105542 | Biomedical Informatics | Database layer (SQLite schema, multi-source fetching, isoforms) |
| Obaid Alabdouli | U23105535 | Computer Science | Web application (Flask backend, REST API) |
| Ahmed Alsamiri | U23103109 | Computer Science | Web application (frontend pages, UI design) |

**Supervised by:** Prof. Osman Abul  
**Co-supervised by:** Prof. Naved  
**University of Sharjah — 2025/2026**