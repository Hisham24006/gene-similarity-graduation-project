"""
motif_similarity.py
-------------------
Computes gene similarity based on shared transcription factor binding motifs
using JASPAR motif profiles via pyJASPAR and Biopython motif scanning.

Pipeline:
    1. Load all human TF motifs from JASPAR via pyJASPAR
    2. Fetch promoter region (1000 bp upstream) from Ensembl for each gene
    3. Scan promoter against all JASPAR motifs using Biopython PSSMs
    4. Compare two genes by their shared TF binding site sets (Jaccard similarity)

Functions:
    load_jaspar_motifs()                 -> list of Biopython motif objects
    fetch_promoter(gene_symbol)          -> DNA string or None
    scan_motifs(promoter_seq, motifs)    -> set of TF names that match
    jaccard_motif_similarity(tfs_a, tfs_b) -> float 0.0 to 1.0
    build_motif_cache(gene_list)         -> dict {symbol: set of TF names}
"""

import os
import ssl
import time
import json
import requests
from pyjaspar import jaspardb
from Bio.Seq import Seq

# -----------------------------------------------------------------------
# SSL fix
# -----------------------------------------------------------------------
if not os.environ.get('PYTHONHTTPSVERIFY', '') and getattr(ssl, '_create_unverified_context', None):
    ssl._create_default_https_context = ssl._create_unverified_context

# -----------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------

PROMOTER_BP              = 1000   # bp upstream of TSS to fetch
PSEUDOCOUNT              = 0.5    # avoids log(0) in PSSM
SCORE_THRESHOLD_FRACTION = 0.80   # hit must score >= 80% of max possible

# Cache file paths
CACHE_DIR      = os.path.join(os.path.dirname(__file__), '..', 'data')
PROMOTER_CACHE = os.path.join(CACHE_DIR, 'promoter_cache.json')
TF_CACHE       = os.path.join(CACHE_DIR, 'tf_cache.json')


# -----------------------------------------------------------------------
# JASPAR motif loading via pyJASPAR
# -----------------------------------------------------------------------

def load_jaspar_motifs():
    """
    Loads all human TF motifs from JASPAR 2024 using pyJASPAR.
    Returns a list of (tf_name, pssm, max_score) tuples ready for scanning.
    """
    print("  [JASPAR] Loading human TF motifs via pyJASPAR...")

    jdb      = jaspardb(release="JASPAR2024")
    motif_list = jdb.fetch_motifs(
        collection="CORE",
        tax_group=["vertebrates"]
    )

    bio_motifs = []
    for m in motif_list:
        try:
            pssm = m.counts.normalize(pseudocounts=PSEUDOCOUNT).log_odds()

            # Max possible score = best base at each position
            max_score = sum(
                max(pssm[base][i] for base in "ACGT")
                for i in range(m.length)
            )

            if max_score > 0:
                bio_motifs.append((m.name, pssm, max_score))

        except Exception:
            continue

    print(f"  [JASPAR] Loaded {len(bio_motifs)} motif models.")
    return bio_motifs


# -----------------------------------------------------------------------
# Promoter fetching
# -----------------------------------------------------------------------

def fetch_promoter(gene_symbol):
    """
    Fetches the promoter region (PROMOTER_BP bases upstream of TSS)
    for a human gene using the Ensembl REST API.
    Returns a DNA string (uppercase) or None if not found.
    Caches results to data/promoter_cache.json.
    """
    os.makedirs(CACHE_DIR, exist_ok=True)

    cache = {}
    if os.path.exists(PROMOTER_CACHE):
        with open(PROMOTER_CACHE, 'r') as f:
            cache = json.load(f)

    if gene_symbol in cache:
        return cache[gene_symbol]

    print(f"  [Ensembl] Fetching promoter for {gene_symbol}...")
    try:
        # Step 1: get gene coordinates
        url_lookup = (
            f"https://rest.ensembl.org/lookup/symbol/homo_sapiens/{gene_symbol}"
            f"?content-type=application/json"
        )
        res = requests.get(url_lookup, timeout=15)
        if res.status_code != 200:
            print(f"  [Ensembl] Gene not found: {gene_symbol}")
            return None

        info   = res.json()
        chrom  = info.get("seq_region_name")
        strand = info.get("strand")
        tss    = info.get("start") if strand == 1 else info.get("end")

        if not chrom or not tss:
            return None

        # Step 2: fetch upstream sequence
        if strand == 1:
            seq_start = max(1, tss - PROMOTER_BP)
            seq_end   = tss
        else:
            seq_start = tss
            seq_end   = tss + PROMOTER_BP

        url_seq = (
            f"https://rest.ensembl.org/sequence/region/human"
            f"/{chrom}:{seq_start}..{seq_end}:{strand}"
            f"?content-type=application/json"
        )
        res_seq = requests.get(url_seq, timeout=15)
        if res_seq.status_code != 200:
            print(f"  [Ensembl] Could not fetch sequence for {gene_symbol}")
            return None

        seq = res_seq.json().get("seq", "").upper()
        if seq:
            cache[gene_symbol] = seq
            with open(PROMOTER_CACHE, 'w') as f:
                json.dump(cache, f)
            print(f"  [Ensembl] Got promoter for {gene_symbol} ({len(seq)} bp)")
            return seq

    except Exception as e:
        print(f"  [Ensembl] Error for {gene_symbol}: {e}")

    return None


# -----------------------------------------------------------------------
# Motif scanning
# -----------------------------------------------------------------------

def scan_motifs(promoter_seq, bio_motifs):
    """
    Scans a promoter sequence against all JASPAR motifs using Biopython PSSMs.
    A TF is considered to bind if its score exceeds SCORE_THRESHOLD_FRACTION
    of the max possible score at any position on either strand.
    Returns a set of TF names found in the promoter.
    """
    if not promoter_seq or len(promoter_seq) < 10:
        return set()

    clean_seq = ''.join(c for c in promoter_seq.upper() if c in "ACGT")
    if len(clean_seq) < 10:
        return set()

    seq    = Seq(clean_seq)
    seq_rc = seq.reverse_complement()
    found  = set()

    for tf_name, pssm, max_score in bio_motifs:
        threshold = SCORE_THRESHOLD_FRACTION * max_score
        motif_len = pssm.length

        if len(clean_seq) < motif_len:
            continue

        hit = False

        # Scan forward strand
        for i in range(len(clean_seq) - motif_len + 1):
            subseq = str(seq[i:i + motif_len])
            if any(c not in "ACGT" for c in subseq):
                continue
            try:
                score = sum(pssm[subseq[j]][j] for j in range(motif_len))
                if score >= threshold:
                    hit = True
                    break
            except Exception:
                continue

        # Scan reverse strand if not found
        if not hit:
            for i in range(len(clean_seq) - motif_len + 1):
                subseq = str(seq_rc[i:i + motif_len])
                if any(c not in "ACGT" for c in subseq):
                    continue
                try:
                    score = sum(pssm[subseq[j]][j] for j in range(motif_len))
                    if score >= threshold:
                        hit = True
                        break
                except Exception:
                    continue

        if hit:
            found.add(tf_name)

    return found


# -----------------------------------------------------------------------
# Motif similarity
# -----------------------------------------------------------------------

def cached_scan_motifs(gene_symbol, promoter_seq, bio_motifs):
    """
    Scans motifs for a gene and caches the result to tf_cache.json.
    On subsequent runs, loads from cache instead of re-scanning.
    Returns a set of TF names.
    """
    os.makedirs(CACHE_DIR, exist_ok=True)

    # Load existing TF cache
    tf_cache = {}
    if os.path.exists(TF_CACHE):
        with open(TF_CACHE, 'r') as f:
            tf_cache = json.load(f)

    # Return cached result if available
    if gene_symbol in tf_cache:
        return set(tf_cache[gene_symbol])

    # Scan and cache
    tfs = scan_motifs(promoter_seq, bio_motifs)
    tf_cache[gene_symbol] = list(tfs)  # sets aren't JSON serializable
    with open(TF_CACHE, 'w') as f:
        json.dump(tf_cache, f)

    return tfs
    """
    Jaccard similarity between two sets of TF names.
    similarity = |shared TFs| / |total unique TFs|
    Returns 0.0 if both sets are empty, 1.0 if identical.
    """
    if not tfs_a and not tfs_b:
        return 0.0
    if not tfs_a or not tfs_b:
        return 0.0
    return len(tfs_a & tfs_b) / len(tfs_a | tfs_b)

def jaccard_motif_similarity(tfs_a, tfs_b):
    """
    Jaccard similarity between two sets of TF names.
    similarity = |shared TFs| / |total unique TFs|
    Returns 0.0 if both sets are empty, 1.0 if identical.
    """
    if not tfs_a and not tfs_b:
        return 0.0
    if not tfs_a or not tfs_b:
        return 0.0
    return len(tfs_a & tfs_b) / len(tfs_a | tfs_b)

# -----------------------------------------------------------------------
# Batch cache builder
# -----------------------------------------------------------------------

def build_motif_cache(gene_list):
    """
    Fetches promoters and scans motifs for all genes in the list.
    Returns a dict: {gene_symbol: set of TF names}
    Promoters are cached automatically to data/promoter_cache.json.
    """
    print("\n--- Loading JASPAR motifs ---")
    bio_motifs = load_jaspar_motifs()

    print(f"\n--- Scanning promoters for {len(gene_list)} genes ---")
    motif_map = {}

    for i, symbol in enumerate(gene_list, 1):
        print(f"[{i}/{len(gene_list)}] {symbol}")
        promoter = fetch_promoter(symbol)
        time.sleep(0.3)

        if promoter:
            tfs = cached_scan_motifs(symbol, promoter, bio_motifs)
            motif_map[symbol] = tfs
            print(f"  Found {len(tfs)} TF binding sites")
        else:
            motif_map[symbol] = set()
            print(f"  No promoter found, skipping.")

    return motif_map


# -----------------------------------------------------------------------
# Quick test
# -----------------------------------------------------------------------

if __name__ == "__main__":
    test_genes = ["TP53", "BRCA1", "KRAS", "MYC"]

    motif_map = build_motif_cache(test_genes)

    print("\n--- Motif similarity matrix ---")
    print(f"{'':10}", end="")
    for g in test_genes:
        print(f"{g:>10}", end="")
    print()

    for g1 in test_genes:
        print(f"{g1:<10}", end="")
        for g2 in test_genes:
            score = jaccard_motif_similarity(motif_map[g1], motif_map[g2])
            print(f"{score:>10.4f}", end="")
        print()