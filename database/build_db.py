"""
build_db.py
-----------
Fetches DNA sequences, protein sequences, and metadata for each gene
from NCBI and UniProt, then stores everything in gene_vault.db.

Run this script ONCE to build your local database:
    python build_db.py

After it finishes, gene_vault.db will appear in the database/ folder.
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

GENES = [
    "TP53", "BRCA1", "KRAS", "EGFR", "HBB",
    "GAPDH", "MDM2", "HRAS", "MYC", "AKT1",
    "CDK2", "CDK4", "MTOR", "PTEN", "RB1",
    "VEGFA", "INS", "ALB", "TNF", "IL6"
]

ORGANISM = "Homo sapiens"

# -----------------------------------------------------------------------
# Fetching functions
# -----------------------------------------------------------------------

def fetch_protein_ncbi(symbol):
    """
    Fetches the top RefSeq protein sequence for a human gene from NCBI.
    Returns (sequence_str, description, length) or (None, None, None).
    """
    try:
        handle = Entrez.esearch(
            db="protein",
            term=f"{symbol}[Gene Name] AND {ORGANISM}[Organism] AND srcdb_refseq[PROP]",
            retmax=1
        )
        record = Entrez.read(handle)
        handle.close()

        if not record["IdList"]:
            print(f"  [NCBI Protein] No results for {symbol}")
            return None, None, None

        pid = record["IdList"][0]
        handle = Entrez.efetch(db="protein", id=pid, rettype="fasta", retmode="text")
        fasta_data = handle.read()
        handle.close()

        seq_record = SeqIO.read(StringIO(fasta_data), "fasta")
        seq_str = str(seq_record.seq)
        description = seq_record.description
        return seq_str, description, len(seq_str)

    except Exception as e:
        print(f"  [NCBI Protein] Error for {symbol}: {e}")
        return None, None, None


def fetch_dna_ncbi(symbol):
    """
    Fetches the top RefSeq mRNA (nucleotide) sequence for a human gene from NCBI.
    Returns (sequence_str, length) or (None, None).
    """
    try:
        handle = Entrez.esearch(
            db="nucleotide",
            term=f"{symbol}[Gene Name] AND {ORGANISM}[Organism] AND mRNA[Filter] AND srcdb_refseq[PROP]",
            retmax=1
        )
        record = Entrez.read(handle)
        handle.close()

        if not record["IdList"]:
            print(f"  [NCBI DNA]     No results for {symbol}")
            return None, None

        nid = record["IdList"][0]
        handle = Entrez.efetch(db="nucleotide", id=nid, rettype="fasta", retmode="text")
        fasta_data = handle.read()
        handle.close()

        seq_record = SeqIO.read(StringIO(fasta_data), "fasta")
        seq_str = str(seq_record.seq)
        return seq_str, len(seq_str)

    except Exception as e:
        print(f"  [NCBI DNA]     Error for {symbol}: {e}")
        return None, None


def fetch_protein_uniprot(symbol):
    """
    Fetches the top reviewed (Swiss-Prot) protein sequence for a human gene from UniProt.
    Returns sequence_str or None.
    """
    try:
        url = (
            f"https://rest.uniprot.org/uniprotkb/search"
            f"?query=gene:{symbol}+AND+organism_id:9606+AND+reviewed:true"
            f"&format=fasta&size=1"
        )
        response = requests.get(url, timeout=10)

        if response.status_code == 200 and response.text.strip():
            seq_record = SeqIO.read(StringIO(response.text), "fasta")
            return str(seq_record.seq)
        else:
            print(f"  [UniProt]      No results for {symbol}")
            return None

    except Exception as e:
        print(f"  [UniProt]      Error for {symbol}: {e}")
        return None


# -----------------------------------------------------------------------
# Main build loop
# -----------------------------------------------------------------------

def build_database():
    db = GeneDatabase(db_path="gene_vault.db")

    print("=" * 50)
    print("  Building gene_vault.db")
    print(f"  {len(GENES)} genes to process")
    print("=" * 50)

    for i, symbol in enumerate(GENES, 1):
        print(f"\n[{i}/{len(GENES)}] Processing {symbol}...")

        # Check which sequences are missing (not an all-or-nothing skip)
        need_ncbi_protein = not db.sequence_exists(symbol, "protein", "NCBI")
        need_ncbi_dna     = not db.sequence_exists(symbol, "dna",     "NCBI")
        need_uniprot      = not db.sequence_exists(symbol, "protein", "UniProt")

        if not any([need_ncbi_protein, need_ncbi_dna, need_uniprot]):
            print(f"  All sequences already present, skipping.")
            continue

        # 1. Fetch protein from NCBI (only if missing)
        protein_seq, description, protein_len = None, None, None
        if need_ncbi_protein:
            protein_seq, description, protein_len = fetch_protein_ncbi(symbol)
            time.sleep(0.4)
        else:
            print(f"  [NCBI Protein] Already in DB, skipping.")

        # 2. Fetch DNA from NCBI (only if missing)
        dna_seq, dna_len = None, None
        if need_ncbi_dna:
            dna_seq, dna_len = fetch_dna_ncbi(symbol)
            time.sleep(0.4)
        else:
            print(f"  [NCBI DNA]     Already in DB, skipping.")

        # 3. Fetch protein from UniProt (only if missing)
        uniprot_seq = None
        if need_uniprot:
            uniprot_seq = fetch_protein_uniprot(symbol)
            time.sleep(0.2)
        else:
            print(f"  [UniProt]      Already in DB, skipping.")

        # 4. Save gene metadata (INSERT OR IGNORE — safe to call again)
        gene_id = db.save_gene(
            symbol=symbol,
            organism=ORGANISM,
            description=description or "N/A",
            dna_length=dna_len or 0,
            protein_length=protein_len or 0
        )

        # 5. Save any newly fetched sequences
        if protein_seq:
            db.save_sequence(gene_id, "protein", "NCBI", protein_seq)
            print(f"  [NCBI Protein] Saved ({protein_len} AA)")

        if dna_seq:
            db.save_sequence(gene_id, "dna", "NCBI", dna_seq)
            print(f"  [NCBI DNA]     Saved ({dna_len} bp)")

        if uniprot_seq:
            db.save_sequence(gene_id, "protein", "UniProt", uniprot_seq)
            print(f"  [UniProt]      Saved ({len(uniprot_seq)} AA)")

    print("\n" + "=" * 50)
    db.summary()
    db.close()
    print("  gene_vault.db is ready.")
    print("=" * 50)


if __name__ == "__main__":
    build_database()
