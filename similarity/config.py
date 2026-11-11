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

VALID_AA = set("ACDEFGHIKLMNPQRSTVWYacdefghiklmnpqrstvwy")
VALID_DNA = set("ATGCatgcNn")

