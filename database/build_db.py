"""
build_db.py
-----------
Fetches DNA and protein sequences + metadata for each gene from
NCBI, UniProt, and Ensembl, then stores everything in gene_vault.db.

Schema:
    genes     — id, symbol, organism, description
    sequences — id, gene_id, type, source, isoform_id, sequence, length, description

Run this script to build or expand your local database:
    python build_db.py

Safe to re-run — already stored isoforms are skipped automatically.
"""

import ssl
import os
import time
import requests
from io import StringIO
from Bio import Entrez, SeqIO
from database_manager import GeneDatabase

# -----------------------------------------------------------------------
# SSL fix (needed in some university/corporate networks)
# -----------------------------------------------------------------------
if not os.environ.get('PYTHONHTTPSVERIFY', '') and getattr(ssl, '_create_unverified_context', None):
    ssl._create_default_https_context = ssl._create_unverified_context

# -----------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------
Entrez.email = "hishamalsaadi06@gmail.com"

MAX_ISOFORMS = 5  # How many isoforms to fetch per gene per source

GENES = [
    # --- Original 20 ---
    "TP53",  "BRCA1", "KRAS",  "EGFR",  "HBB",
    "GAPDH", "MDM2",  "HRAS",  "MYC",   "AKT1",
    "CDK2",  "CDK4",  "MTOR",  "PTEN",  "RB1",
    "VEGFA", "INS",   "ALB",   "TNF",   "IL6",

    # --- Cancer related ---
    "BRCA2", "ATM",   "BRAF",  "PIK3CA","ERBB2",
    "NRAS",  "ABL1",  "BCL2",  "CDKN2A","FGFR1",

    # --- Cell signaling ---
    "SRC",   "JAK2",  "STAT3", "RAF1",  "MAP2K1",
    "MAP2K2","MAPK1", "MAPK3", "AKT2",  "AKT3",

    # --- Transcription factors ---
    "MYOD1", "E2F1",  "SP1",   "FOXP3", "GATA1",
    "TP63",  "TP73",  "PCNA",  "CCND1", "CCNE1",

    # --- DNA repair ---
    "MLH1",  "MSH2",  "RAD51", "XRCC1", "ERCC1",
    "CHEK1", "CHEK2", "FANCA", "PALB2", "NBN"
]

ORGANISM = "Homo sapiens"

# -----------------------------------------------------------------------
# Fetching functions
# -----------------------------------------------------------------------

def fetch_proteins_ncbi(symbol, max_isoforms=MAX_ISOFORMS):
    """
    Fetches up to max_isoforms RefSeq protein isoforms from NCBI.
    Returns list of (isoform_id, sequence, description).
    """
    results = []
    try:
        handle = Entrez.esearch(
            db="protein",
            term=f"{symbol}[Gene Name] AND {ORGANISM}[Organism] AND srcdb_refseq[PROP]",
            retmax=max_isoforms
        )
        record = Entrez.read(handle)
        handle.close()

        if not record["IdList"]:
            print(f"  [NCBI Protein] No results for {symbol}")
            return results

        ids = record["IdList"]
        handle = Entrez.efetch(db="protein", id=",".join(ids), rettype="fasta", retmode="text")
        fasta_data = handle.read()
        handle.close()

        for seq_record in SeqIO.parse(StringIO(fasta_data), "fasta"):
            results.append((seq_record.id, str(seq_record.seq), seq_record.description))

    except Exception as e:
        print(f"  [NCBI Protein] Error for {symbol}: {e}")

    return results


def fetch_dna_ncbi(symbol, max_isoforms=MAX_ISOFORMS):
    """
    Fetches up to max_isoforms RefSeq mRNA isoforms from NCBI.
    Returns list of (isoform_id, sequence, description).
    """
    results = []
    try:
        handle = Entrez.esearch(
            db="nucleotide",
            term=f"{symbol}[Gene Name] AND {ORGANISM}[Organism] AND mRNA[Filter] AND srcdb_refseq[PROP]",
            retmax=max_isoforms
        )
        record = Entrez.read(handle)
        handle.close()

        if not record["IdList"]:
            print(f"  [NCBI DNA]     No results for {symbol}")
            return results

        ids = record["IdList"]
        handle = Entrez.efetch(db="nucleotide", id=",".join(ids), rettype="fasta", retmode="text")
        fasta_data = handle.read()
        handle.close()

        for seq_record in SeqIO.parse(StringIO(fasta_data), "fasta"):
            results.append((seq_record.id, str(seq_record.seq), seq_record.description))

    except Exception as e:
        print(f"  [NCBI DNA]     Error for {symbol}: {e}")

    return results


def fetch_proteins_uniprot(symbol, max_isoforms=MAX_ISOFORMS):
    """
    Fetches up to max_isoforms reviewed UniProt protein isoforms.
    Returns list of (isoform_id, sequence, description).
    """
    results = []
    try:
        url = (
            f"https://rest.uniprot.org/uniprotkb/search"
            f"?query=gene_exact:{symbol}+AND+organism_id:9606+AND+reviewed:true"
            f"&format=fasta&size={max_isoforms}"
        )
        response = requests.get(url, timeout=15)

        if response.status_code == 200 and response.text.strip():
            for seq_record in SeqIO.parse(StringIO(response.text), "fasta"):
                results.append((seq_record.id, str(seq_record.seq), seq_record.description))
        else:
            print(f"  [UniProt]      No results for {symbol}")

    except Exception as e:
        print(f"  [UniProt]      Error for {symbol}: {e}")

    return results


def fetch_dna_ensembl(symbol):
    """
    Fetches the primary transcript DNA sequence from Ensembl.
    Returns (ensembl_id, sequence, description) or (None, None, None).
    """
    try:
        url_lookup = (
            f"https://rest.ensembl.org/lookup/symbol/homo_sapiens/{symbol}"
            f"?content-type=application/json"
        )
        res = requests.get(url_lookup, timeout=15)
        if res.status_code != 200:
            print(f"  [Ensembl]      No gene ID found for {symbol}")
            return None, None, None

        data = res.json()
        ensembl_id = data.get("id")
        display_name = data.get("display_name", symbol)
        if not ensembl_id:
            return None, None, None

        url_seq = (
            f"https://rest.ensembl.org/sequence/id/{ensembl_id}"
            f"?content-type=application/json"
        )
        res_seq = requests.get(url_seq, timeout=15)
        if res_seq.status_code == 200:
            seq = res_seq.json().get("seq")
            if seq:
                description = f"{display_name} | {ensembl_id} | Homo sapiens"
                return ensembl_id, seq, description

        print(f"  [Ensembl]      No sequence for {symbol}")
        return None, None, None

    except Exception as e:
        print(f"  [Ensembl]      Error for {symbol}: {e}")
        return None, None, None


# -----------------------------------------------------------------------
# Main build loop
# -----------------------------------------------------------------------

def build_database():
    db = GeneDatabase(db_path="gene_vault.db")

    print("=" * 55)
    print("  Building gene_vault.db")
    print(f"  {len(GENES)} genes  |  up to {MAX_ISOFORMS} isoforms per source")
    print(f"  Sources: NCBI Protein, NCBI DNA, UniProt, Ensembl")
    print("=" * 55)

    for i, symbol in enumerate(GENES, 1):
        print(f"\n[{i}/{len(GENES)}] {symbol}")

        # ---- NCBI Protein isoforms ----
        existing = db.count_isoforms(symbol, "protein", "NCBI")
        if existing >= MAX_ISOFORMS:
            print(f"  [NCBI Protein] {existing} isoforms already in DB, skipping.")
        else:
            isoforms = fetch_proteins_ncbi(symbol)
            time.sleep(0.4)
            saved, gene_desc = 0, "N/A"
            for iso_id, seq, desc in isoforms:
                if not db.isoform_exists(symbol, "protein", "NCBI", iso_id):
                    if gene_desc == "N/A" and desc:
                        gene_desc = desc
                    gene_id = db.save_gene(symbol, ORGANISM, gene_desc)
                    db.save_sequence(gene_id, "protein", "NCBI", seq,
                                     isoform_id=iso_id, description=desc)
                    saved += 1
            print(f"  [NCBI Protein] Saved {saved} new isoform(s)" if saved
                  else f"  [NCBI Protein] All isoforms already in DB.")

        # ---- NCBI DNA isoforms ----
        existing = db.count_isoforms(symbol, "dna", "NCBI")
        if existing >= MAX_ISOFORMS:
            print(f"  [NCBI DNA]     {existing} isoforms already in DB, skipping.")
        else:
            isoforms = fetch_dna_ncbi(symbol)
            time.sleep(0.4)
            saved = 0
            for iso_id, seq, desc in isoforms:
                if not db.isoform_exists(symbol, "dna", "NCBI", iso_id):
                    gene_id = db.save_gene(symbol, ORGANISM, "N/A")
                    db.save_sequence(gene_id, "dna", "NCBI", seq,
                                     isoform_id=iso_id, description=desc)
                    saved += 1
            print(f"  [NCBI DNA]     Saved {saved} new isoform(s)" if saved
                  else f"  [NCBI DNA]     All isoforms already in DB.")

        # ---- UniProt Protein isoforms ----
        existing = db.count_isoforms(symbol, "protein", "UniProt")
        if existing >= MAX_ISOFORMS:
            print(f"  [UniProt]      {existing} isoforms already in DB, skipping.")
        else:
            isoforms = fetch_proteins_uniprot(symbol)
            time.sleep(0.3)
            saved = 0
            for iso_id, seq, desc in isoforms:
                if not db.isoform_exists(symbol, "protein", "UniProt", iso_id):
                    gene_id = db.save_gene(symbol, ORGANISM, "N/A")
                    db.save_sequence(gene_id, "protein", "UniProt", seq,
                                     isoform_id=iso_id, description=desc)
                    saved += 1
            print(f"  [UniProt]      Saved {saved} new isoform(s)" if saved
                  else f"  [UniProt]      All isoforms already in DB.")

        # ---- Ensembl DNA ----
        existing = db.count_isoforms(symbol, "dna", "Ensembl")
        if existing > 0:
            print(f"  [Ensembl]      Already in DB, skipping.")
        else:
            ens_id, seq, desc = fetch_dna_ensembl(symbol)
            time.sleep(0.3)
            if ens_id and seq:
                gene_id = db.save_gene(symbol, ORGANISM, "N/A")
                db.save_sequence(gene_id, "dna", "Ensembl", seq,
                                 isoform_id=ens_id, description=desc)
                print(f"  [Ensembl]      Saved ({len(seq)} bp)")
            else:
                print(f"  [Ensembl]      No sequence found.")

    print("\n" + "=" * 55)
    db.summary()
    db.close()
    print("  gene_vault.db is ready.")
    print("=" * 55)


if __name__ == "__main__":
    build_database()
