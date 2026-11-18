"""
import_all_human.py
-------------------
Imports EVERY reviewed human protein from UniProt (~20,000 genes) into
gene_vault.db, with an automatic category from UniProt keywords.

Usage (from the database/ folder):
    python import_all_human.py              # all reviewed human genes
    python import_all_human.py --limit 500  # test with the first 500 only
"""
import argparse
import csv
import io
import os
import sqlite3
import ssl

import requests

from database_manager import GeneDatabase

if not os.environ.get('PYTHONHTTPSVERIFY', '') and getattr(ssl, '_create_unverified_context', None):
    ssl._create_default_https_context = ssl._create_unverified_context

HERE = os.path.dirname(os.path.abspath(__file__))
DB_PATH = os.path.join(HERE, "gene_vault.db")
ORGANISM = "Homo sapiens"

URL = "https://rest.uniprot.org/uniprotkb/stream"
PARAMS = {
    "query": "organism_id:9606 AND reviewed:true",
    "format": "tsv",
    "fields": "accession,gene_primary,protein_name,sequence,keyword",
}

# First match wins, so put the most specific categories first.
CATEGORY_RULES = [
    ("DNA repair",            "DNA repair"),
    ("Cancer-related",        "Proto-oncogene"),
    ("Cancer-related",        "Tumor suppressor"),
    ("Transcription factors", "Transcription regulation"),
    ("Cell signaling",        "Kinase"),
    ("Cell signaling",        "Receptor"),
]


def categorize(keywords):
    for category, kw in CATEGORY_RULES:
        if kw in keywords:
            return category
    return "Other"


def main(limit):
    print("Downloading reviewed human proteins from UniProt...")
    r = requests.get(URL, params=PARAMS, timeout=600)
    r.raise_for_status()
    rows = csv.reader(io.StringIO(r.text), delimiter="\t")
    next(rows)  # header

    db = GeneDatabase(db_path=DB_PATH)
    categories, added, skipped = {}, 0, 0

    for n, row in enumerate(rows, 1):
        if limit and n > limit:
            break
        if len(row) < 5:
            continue
        accession, symbol, name, seq, keywords = row[:5]
        symbol = symbol.strip()
        if not symbol or len(seq) < 10:
            skipped += 1
            continue
        if db.isoform_exists(symbol, "protein", "UniProt", accession):
            skipped += 1
            continue

        gene_id = db.save_gene(symbol, ORGANISM, name or "N/A")
        db.save_sequence(gene_id, "protein", "UniProt", seq,
                         isoform_id=accession, description=name)
        categories[symbol] = categorize(keywords)
        added += 1
        if added % 1000 == 0:
            print(f"  {added} proteins imported...")

    db.summary()
    db.close()

    # Add / fill the category column (separate connection, after db.close())
    conn = sqlite3.connect(DB_PATH, timeout=30)
    cols = [c[1] for c in conn.execute("PRAGMA table_info(genes)")]
    if "category" not in cols:
        conn.execute("ALTER TABLE genes ADD COLUMN category TEXT")
    conn.executemany(
        "UPDATE genes SET category=? WHERE symbol=? AND (category IS NULL OR category='')",
        [(cat, sym) for sym, cat in categories.items()])
    conn.execute("CREATE INDEX IF NOT EXISTS idx_genes_symbol ON genes(symbol)")
    conn.execute("CREATE INDEX IF NOT EXISTS idx_seq_gene ON sequences(gene_id)")
    conn.commit()
    conn.close()

    print(f"\nDone. Added {added} proteins, skipped {skipped}.")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=0, help="only import the first N (for testing)")
    main(ap.parse_args().limit)