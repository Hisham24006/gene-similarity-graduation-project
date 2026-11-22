"""
import_dna_human.py
-----------------
Adds one DNA (cDNA) sequence per gene to gene_vault.db, using Ensembl's
human transcript file. For each gene already in the database it keeps
the longest protein-coding transcript.

Usage (from the database/ folder):
    python import_dna_human.py                  # every gene in the database
    python import_dna_human.py --limit 200      # test with 200 genes
    python import_dna_human.py --max-len 10000  # skip transcripts longer than this
"""
import argparse
import gzip
import os
import re
import ssl

import requests
from Bio import SeqIO

from database_manager import GeneDatabase

if not os.environ.get('PYTHONHTTPSVERIFY', '') and getattr(ssl, '_create_unverified_context', None):
    ssl._create_default_https_context = ssl._create_unverified_context

HERE = os.path.dirname(os.path.abspath(__file__))
DB_PATH = os.path.join(HERE, "gene_vault.db")
DATA_DIR = os.path.normpath(os.path.join(HERE, "..", "data"))
GZ_PATH = os.path.join(DATA_DIR, "Homo_sapiens.GRCh38.cdna.all.fa.gz")

URL = ("https://ftp.ensembl.org/pub/current_fasta/homo_sapiens/cdna/"
       "Homo_sapiens.GRCh38.cdna.all.fa.gz")

SYMBOL_RE = re.compile(r"gene_symbol:(\S+)")


def download():
    if os.path.exists(GZ_PATH):
        print(f"Using file already downloaded: {GZ_PATH}")
        return
    os.makedirs(DATA_DIR, exist_ok=True)
    print("Downloading Ensembl human transcripts (one time, a few minutes)...")
    tmp = GZ_PATH + ".part"
    with requests.get(URL, stream=True, timeout=(30, 600)) as r:
        r.raise_for_status()
        with open(tmp, "wb") as f:
            for chunk in r.iter_content(1 << 20):
                f.write(chunk)
    os.replace(tmp, GZ_PATH)
    print("Download complete.")


def main(limit, max_len):
    db = GeneDatabase(db_path=DB_PATH)
    symbols = {row[0] for row in db.conn.execute("SELECT symbol FROM genes")}
    print(f"{len(symbols)} genes in the database.")

    download()

    print("Reading transcripts and picking the best one per gene...")
    best = {}
    with gzip.open(GZ_PATH, "rt") as f:
        for rec in SeqIO.parse(f, "fasta"):
            desc = rec.description
            m = SYMBOL_RE.search(desc)
            if not m:
                continue
            sym = m.group(1)
            if sym not in symbols:
                continue
            if "transcript_biotype:protein_coding" not in desc:
                continue
            n = len(rec.seq)
            if n < 50 or n > max_len:
                continue
            if sym not in best or n > best[sym][3]:
                best[sym] = (rec.id, str(rec.seq).upper(), desc, n)

    print(f"Found DNA for {len(best)} genes.")

    added = skipped = 0
    for sym in sorted(best):
        if limit and added >= limit:
            break
        iso_id, seq, desc, n = best[sym]
        if db.isoform_exists(sym, "dna", "Ensembl", iso_id):
            skipped += 1
            continue
        row = db.conn.execute("SELECT id FROM genes WHERE symbol=?", (sym,)).fetchone()
        if not row:
            continue
        label = f"{sym} | {iso_id} | longest protein-coding transcript cDNA | Homo sapiens"
        db.save_sequence(row[0], "dna", "Ensembl", seq,
                         isoform_id=iso_id, description=label)
        added += 1
        if added % 1000 == 0:
            print(f"  {added} DNA sequences added...")

    db.conn.commit()
    db.summary()
    db.close()
    print(f"\nDone. Added {added} DNA sequences, skipped {skipped}.")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=0, help="only add this many genes (for testing)")
    ap.add_argument("--max-len", type=int, default=10000,
                    help="skip transcripts longer than this many bases (default 10000)")
    args = ap.parse_args()
    main(args.limit, args.max_len)