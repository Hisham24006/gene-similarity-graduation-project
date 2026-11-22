"""
trim_db.py
----------
Shrinks gene_vault.db below a target size by deleting whole genes
(protein + DNA together), largest first. The original 60 genes are kept.

Usage (from the database/ folder):
    python trim_db.py              # shrink to about 90 MB
    python trim_db.py --target 80  # shrink to about 80 MB
"""
import argparse
import os
import sqlite3

HERE = os.path.dirname(os.path.abspath(__file__))
DB_PATH = os.path.join(HERE, "gene_vault.db")

PROTECTED = {
    "TP53", "BRCA1", "KRAS", "EGFR", "HBB", "GAPDH", "MDM2", "HRAS", "MYC", "AKT1",
    "CDK2", "CDK4", "MTOR", "PTEN", "RB1", "VEGFA", "INS", "ALB", "TNF", "IL6",
    "BRCA2", "ATM", "BRAF", "PIK3CA", "ERBB2", "NRAS", "ABL1", "BCL2", "CDKN2A", "FGFR1",
    "SRC", "JAK2", "STAT3", "RAF1", "MAP2K1", "MAP2K2", "MAPK1", "MAPK3", "AKT2", "AKT3",
    "MYOD1", "E2F1", "SP1", "FOXP3", "GATA1", "TP63", "TP73", "PCNA", "CCND1", "CCNE1",
    "MLH1", "MSH2", "RAD51", "XRCC1", "ERCC1", "CHEK1", "CHEK2", "FANCA", "PALB2", "NBN",
}


def counts(conn):
    g = conn.execute("SELECT COUNT(*) FROM genes").fetchone()[0]
    p = conn.execute("SELECT COUNT(*) FROM sequences WHERE type='protein'").fetchone()[0]
    d = conn.execute("SELECT COUNT(*) FROM sequences WHERE type='dna'").fetchone()[0]
    return g, p, d


def main(target_mb):
    size_mb = os.path.getsize(DB_PATH) / 1048576
    print(f"Current size: {size_mb:.1f} MB, target: {target_mb} MB")
    if size_mb <= target_mb:
        print("Already small enough. Nothing to do.")
        return

    conn = sqlite3.connect(DB_PATH, timeout=30)
    print("Before: %d genes, %d protein, %d DNA" % counts(conn))

    rows = conn.execute("""
        SELECT g.id, g.symbol, COALESCE(SUM(LENGTH(s.sequence)), 0) AS bytes
        FROM genes g LEFT JOIN sequences s ON s.gene_id = g.id
        GROUP BY g.id
        ORDER BY bytes DESC
    """).fetchall()

    total_seq_bytes = sum(r[2] for r in rows)
    ratio = (size_mb * 1048576) / max(total_seq_bytes, 1)   # file bytes per sequence byte
    need_seq_bytes = ((size_mb - target_mb) * 1048576) / ratio

    chosen, removed = [], 0
    for gene_id, symbol, nbytes in rows:
        if removed >= need_seq_bytes:
            break
        if symbol in PROTECTED:
            continue
        chosen.append((gene_id,))
        removed += nbytes

    conn.executemany("DELETE FROM sequences WHERE gene_id = ?", chosen)
    conn.executemany("DELETE FROM genes WHERE id = ?", chosen)
    conn.commit()
    print(f"Deleted {len(chosen)} genes. Shrinking the file (this can take a minute)...")
    conn.execute("VACUUM")
    print("After:  %d genes, %d protein, %d DNA" % counts(conn))
    conn.close()

    print(f"New size: {os.path.getsize(DB_PATH) / 1048576:.1f} MB")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--target", type=float, default=90, help="target size in MB (default 90)")
    main(ap.parse_args().target)