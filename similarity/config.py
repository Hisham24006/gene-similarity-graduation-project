"""
config.py
---------
Centralized configuration, paths, and environment settings.
"""

import os
import ssl
from Bio import Entrez

# ---------------------------------------------------------------------------
# SSL Certificate Handling (needed on some academic / corporate networks)
# ---------------------------------------------------------------------------
if not os.environ.get('PYTHONHTTPSVERIFY', '') and getattr(ssl, '_create_unverified_context', None):
    ssl._create_default_https_context = ssl._create_unverified_context

# ---------------------------------------------------------------------------
# Project Directories and File Paths
# ---------------------------------------------------------------------------
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_DIR = os.path.join(BASE_DIR, 'data')
DB_DIR = os.path.join(BASE_DIR, 'database')
SIM_DIR = os.path.join(BASE_DIR, 'similarity')
RESULTS_DIR = os.path.join(BASE_DIR, 'results')

DB_PATH = os.path.join(DB_DIR, 'gene_vault.db')
PROMOTER_CACHE_PATH = os.path.join(DATA_DIR, 'promoter_cache.json')
TF_CACHE_PATH = os.path.join(DATA_DIR, 'tf_cache.json')
RECENT_SEARCHES_PATH = os.path.join(DATA_DIR, 'recent_searches.json')

# ---------------------------------------------------------------------------
# NCBI & Similarity Parameters
# ---------------------------------------------------------------------------
ENTREZ_EMAIL = os.environ.get("ENTREZ_EMAIL", "hishamalsaadi06@gmail.com")
Entrez.email = ENTREZ_EMAIL

KMER_K = 3      # protein k-mers (20^3 = 8,000 possible 3-mers)
KMER_K_DNA = 7  # nucleotide k-mers (4^7 = 16,384 possible; 4^3 = 64 is too few)

# Two-stage search: a k-mer filter keeps this many candidates, and only they
# are scored with the selected metrics (see kmer_index.py).
PROTEIN_CANDIDATES = 80
# 8-mers (65,536 possible) separate relatives far better than 7-mers: an mRNA
# contains about a fifth of all 7-mers, so unrelated mRNAs share many by chance.
# On 12 test queries, 500 candidates kept 13 of the 14 hits scoring >= 0.1.
DNA_FILTER_K = 8          # k-mer size of the nucleotide filter index
DNA_CANDIDATES = 500
PROMOTER_BP = 1000
PSEUDOCOUNT = 0.5
SCORE_THRESHOLD_FRACTION = 0.80

# BLOSUM62 gap costs offered by NCBI BLAST, as (existence, extension).
# A gap of length k costs existence + k * extension, as in BLAST.
BLOSUM62_GAP_COSTS = [
    (11, 2), (10, 2), (9, 2), (8, 2), (7, 2), (6, 2),
    (13, 1), (12, 1), (11, 1), (10, 1), (9, 1),
]
DEFAULT_GAP_COST = (11, 1)  # NCBI BLAST default for BLOSUM62

# The 20 standard amino acids plus:
#   U = selenocysteine (e.g. TXNRD1, GPX1); scored as C by metrics.for_blosum
#   X = unknown, B = D or N, Z = E or Q (ambiguity codes BLOSUM62 scores)
# A trailing '*' (stop) is removed from pasted queries before this check.
VALID_AA = set("ACDEFGHIKLMNPQRSTUVWYXBZacdefghiklmnpqrstuvwyxbz")
VALID_DNA = set("ATGCatgcNn")

