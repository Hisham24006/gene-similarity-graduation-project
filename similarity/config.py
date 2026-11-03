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

KMER_K = 3
PROMOTER_BP = 1000
PSEUDOCOUNT = 0.5
SCORE_THRESHOLD_FRACTION = 0.80

VALID_AA = set("ACDEFGHIKLMNPQRSTVWYacdefghiklmnpqrstvwy")
VALID_DNA = set("ATGCatgcNn")

