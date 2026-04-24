import sqlite3


class GeneDatabase:
    """
    Handles all read/write operations on gene_vault.db.

    Schema:
        genes     — one row per gene: id, symbol, organism, description
        sequences — one row per isoform: id, gene_id, type, source,
                    isoform_id, sequence, length, description
    """

    def __init__(self, db_path="gene_vault.db"):
        self.conn = sqlite3.connect(db_path)
        self._create_tables()

    # ------------------------------------------------------------------
    # Setup
    # ------------------------------------------------------------------

    def _create_tables(self):
        cursor = self.conn.cursor()

        cursor.execute('''
            CREATE TABLE IF NOT EXISTS genes (
                id          INTEGER PRIMARY KEY AUTOINCREMENT,
                symbol      TEXT    UNIQUE NOT NULL,
                organism    TEXT,
                description TEXT
            )
        ''')

        cursor.execute('''
            CREATE TABLE IF NOT EXISTS sequences (
                id          INTEGER PRIMARY KEY AUTOINCREMENT,
                gene_id     INTEGER NOT NULL,
                type        TEXT    NOT NULL,
                source      TEXT    NOT NULL,
                isoform_id  TEXT,
                sequence    TEXT    NOT NULL,
                length      INTEGER,
                description TEXT,
                FOREIGN KEY (gene_id) REFERENCES genes(id)
            )
        ''')

        self.conn.commit()

    # ------------------------------------------------------------------
    # Writing
    # ------------------------------------------------------------------

    def save_gene(self, symbol, organism, description):
        """
        Insert a gene into the genes table.
        If already exists, updates description if it was previously N/A.
        Returns the gene's id.
        """
        cursor = self.conn.cursor()

        cursor.execute('''
            INSERT OR IGNORE INTO genes (symbol, organism, description)
            VALUES (?, ?, ?)
        ''', (symbol, organism, description))

        cursor.execute('''
            UPDATE genes SET
                description = CASE
                    WHEN description = 'N/A' AND ? != 'N/A' THEN ?
                    ELSE description
                END
            WHERE symbol = ?
        ''', (description, description, symbol))

        self.conn.commit()
        cursor.execute("SELECT id FROM genes WHERE symbol = ?", (symbol,))
        return cursor.fetchone()[0]

    def save_sequence(self, gene_id, seq_type, source, sequence,
                      isoform_id=None, description=None):
        """
        Insert a sequence row linked to a gene.
        Uniqueness is based on gene_id + type + source + isoform_id.
        Skips if that exact isoform already exists.
        """
        cursor = self.conn.cursor()
        cursor.execute('''
            SELECT id FROM sequences
            WHERE gene_id = ? AND type = ? AND source = ? AND isoform_id IS ?
        ''', (gene_id, seq_type, source, isoform_id))

        if cursor.fetchone() is None:
            cursor.execute('''
                INSERT INTO sequences
                    (gene_id, type, source, isoform_id, sequence, length, description)
                VALUES (?, ?, ?, ?, ?, ?, ?)
            ''', (gene_id, seq_type, source, isoform_id,
                  sequence, len(sequence), description))
            self.conn.commit()
            return True
        return False  # already existed

    # ------------------------------------------------------------------
    # Reading
    # ------------------------------------------------------------------

    def get_all_sequences(self, seq_type):
        """
        Returns list of (symbol, sequence) for all isoforms of a given type.
        Used by the similarity engine to load the full local database.
        """
        cursor = self.conn.cursor()
        cursor.execute('''
            SELECT g.symbol, s.sequence
            FROM sequences s
            JOIN genes g ON s.gene_id = g.id
            WHERE s.type = ?
        ''', (seq_type,))
        return cursor.fetchall()

    def get_all_sequences_with_isoform(self, seq_type):
        """
        Returns list of (symbol, isoform_id, sequence) for all isoforms of a given type.
        Used by the similarity engine when isoform-level detail is needed.
        """
        cursor = self.conn.cursor()
        cursor.execute('''
            SELECT g.symbol, s.isoform_id, s.sequence
            FROM sequences s
            JOIN genes g ON s.gene_id = g.id
            WHERE s.type = ?
        ''', (seq_type,))
        return cursor.fetchall()

    def get_isoform_sequence(self, symbol, seq_type, isoform_id):
        """Returns the sequence for a specific isoform. None if not found."""
        cursor = self.conn.cursor()
        cursor.execute('''
            SELECT s.sequence
            FROM sequences s
            JOIN genes g ON s.gene_id = g.id
            WHERE g.symbol = ? AND s.type = ? AND s.isoform_id = ?
            LIMIT 1
        ''', (symbol, seq_type, isoform_id))
        row = cursor.fetchone()
        return row[0] if row else None

    def get_gene_sequence(self, symbol, seq_type):
        """Returns the first sequence of a given type for a gene. None if not found."""
        cursor = self.conn.cursor()
        cursor.execute('''
            SELECT s.sequence
            FROM sequences s
            JOIN genes g ON s.gene_id = g.id
            WHERE g.symbol = ? AND s.type = ?
            LIMIT 1
        ''', (symbol, seq_type))
        row = cursor.fetchone()
        return row[0] if row else None

    def get_gene_info(self, symbol):
        """Returns (symbol, organism, description) for a gene. None if not found."""
        cursor = self.conn.cursor()
        cursor.execute('''
            SELECT symbol, organism, description
            FROM genes WHERE symbol = ?
        ''', (symbol,))
        return cursor.fetchone()

    def get_isoforms(self, symbol, seq_type=None):
        """
        Returns all isoforms for a gene as list of dicts.
        Optionally filter by type ('dna' or 'protein').
        Each dict has: isoform_id, type, source, length, description.
        """
        cursor = self.conn.cursor()
        if seq_type:
            cursor.execute('''
                SELECT s.isoform_id, s.type, s.source, s.length, s.description
                FROM sequences s
                JOIN genes g ON s.gene_id = g.id
                WHERE g.symbol = ? AND s.type = ?
                ORDER BY s.source, s.type
            ''', (symbol, seq_type))
        else:
            cursor.execute('''
                SELECT s.isoform_id, s.type, s.source, s.length, s.description
                FROM sequences s
                JOIN genes g ON s.gene_id = g.id
                WHERE g.symbol = ?
                ORDER BY s.source, s.type
            ''', (symbol,))

        rows = cursor.fetchall()
        return [
            {
                "isoform_id":  r[0],
                "type":        r[1],
                "source":      r[2],
                "length":      r[3],
                "description": r[4]
            }
            for r in rows
        ]

    def list_genes(self):
        """Returns a list of all gene symbols in the database."""
        cursor = self.conn.cursor()
        cursor.execute("SELECT symbol FROM genes ORDER BY symbol")
        return [row[0] for row in cursor.fetchall()]

    def gene_exists(self, symbol):
        """Returns True if a gene is already in the database."""
        cursor = self.conn.cursor()
        cursor.execute("SELECT id FROM genes WHERE symbol = ?", (symbol,))
        return cursor.fetchone() is not None

    def isoform_exists(self, symbol, seq_type, source, isoform_id):
        """Returns True if a specific isoform already exists."""
        cursor = self.conn.cursor()
        cursor.execute('''
            SELECT s.id FROM sequences s
            JOIN genes g ON s.gene_id = g.id
            WHERE g.symbol = ? AND s.type = ? AND s.source = ? AND s.isoform_id = ?
        ''', (symbol, seq_type, source, isoform_id))
        return cursor.fetchone() is not None

    def count_isoforms(self, symbol, seq_type, source):
        """Returns how many isoforms are stored for a gene + type + source."""
        cursor = self.conn.cursor()
        cursor.execute('''
            SELECT COUNT(*) FROM sequences s
            JOIN genes g ON s.gene_id = g.id
            WHERE g.symbol = ? AND s.type = ? AND s.source = ?
        ''', (symbol, seq_type, source))
        return cursor.fetchone()[0]

    # ------------------------------------------------------------------
    # Utility
    # ------------------------------------------------------------------

    def close(self):
        self.conn.close()

    def summary(self):
        """Prints a quick summary of what's in the database."""
        cursor = self.conn.cursor()
        cursor.execute("SELECT COUNT(*) FROM genes")
        gene_count = cursor.fetchone()[0]
        cursor.execute("SELECT COUNT(*) FROM sequences WHERE type = 'protein'")
        protein_count = cursor.fetchone()[0]
        cursor.execute("SELECT COUNT(*) FROM sequences WHERE type = 'dna'")
        dna_count = cursor.fetchone()[0]
        print(f"\n--- Database Summary ---")
        print(f"Genes stored   : {gene_count}")
        print(f"Protein seqs   : {protein_count}")
        print(f"DNA seqs       : {dna_count}")
        print(f"Total seqs     : {protein_count + dna_count}")
        print(f"------------------------\n")