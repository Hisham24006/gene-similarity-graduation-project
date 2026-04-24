import sqlite3


class GeneDatabase:
    """
    Handles all read/write operations on the local gene_vault.db SQLite database.
    Tables:
        genes     — one row per gene (symbol, organism, description)
        sequences — many rows per gene (DNA, protein, from different sources)
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
                description TEXT,
                dna_length  INTEGER,
                protein_length INTEGER
            )
        ''')

        cursor.execute('''
            CREATE TABLE IF NOT EXISTS sequences (
                id       INTEGER PRIMARY KEY AUTOINCREMENT,
                gene_id  INTEGER NOT NULL,
                type     TEXT    NOT NULL,   -- 'dna' or 'protein'
                source   TEXT    NOT NULL,   -- 'NCBI' or 'UniProt'
                sequence TEXT    NOT NULL,
                FOREIGN KEY (gene_id) REFERENCES genes(id)
            )
        ''')

        self.conn.commit()

    # ------------------------------------------------------------------
    # Writing
    # ------------------------------------------------------------------

    def save_gene(self, symbol, organism, description, dna_length, protein_length):
        """
        Insert a gene into the genes table.
        If the gene already exists, update any fields that were previously
        missing (N/A description or 0 lengths) with real values.
        Returns the gene's id.
        """
        cursor = self.conn.cursor()

        # Insert if not already there
        cursor.execute('''
            INSERT OR IGNORE INTO genes (symbol, organism, description, dna_length, protein_length)
            VALUES (?, ?, ?, ?, ?)
        ''', (symbol, organism, description, dna_length, protein_length))

        # Update fields that were previously placeholder values
        cursor.execute('''
            UPDATE genes SET
                description    = CASE WHEN description    = 'N/A' AND ? != 'N/A' THEN ? ELSE description    END,
                protein_length = CASE WHEN protein_length = 0     AND ? != 0     THEN ? ELSE protein_length END,
                dna_length     = CASE WHEN dna_length     = 0     AND ? != 0     THEN ? ELSE dna_length     END
            WHERE symbol = ?
        ''', (description, description, protein_length, protein_length, dna_length, dna_length, symbol))

        self.conn.commit()
        cursor.execute("SELECT id FROM genes WHERE symbol = ?", (symbol,))
        return cursor.fetchone()[0]

    def save_sequence(self, gene_id, seq_type, source, sequence):
        """
        Insert a sequence row linked to a gene.
        Skips if the exact same gene_id + type + source already exists.
        """
        cursor = self.conn.cursor()
        cursor.execute('''
            SELECT id FROM sequences
            WHERE gene_id = ? AND type = ? AND source = ?
        ''', (gene_id, seq_type, source))

        if cursor.fetchone() is None:
            cursor.execute('''
                INSERT INTO sequences (gene_id, type, source, sequence)
                VALUES (?, ?, ?, ?)
            ''', (gene_id, seq_type, source, sequence))
            self.conn.commit()

    # ------------------------------------------------------------------
    # Reading
    # ------------------------------------------------------------------

    def get_all_sequences(self, seq_type):
        """
        Returns a list of (symbol, sequence) for every gene of the given type.
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

    def get_gene_sequence(self, symbol, seq_type):
        """
        Returns the first sequence of a given type for a specific gene symbol.
        Returns None if not found.
        """
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
        """
        Returns metadata (symbol, organism, description, lengths) for a gene.
        Returns None if not found.
        """
        cursor = self.conn.cursor()
        cursor.execute('''
            SELECT symbol, organism, description, dna_length, protein_length
            FROM genes WHERE symbol = ?
        ''', (symbol,))
        return cursor.fetchone()

    def list_genes(self):
        """Returns a list of all gene symbols stored in the database."""
        cursor = self.conn.cursor()
        cursor.execute("SELECT symbol FROM genes ORDER BY symbol")
        return [row[0] for row in cursor.fetchall()]

    def gene_exists(self, symbol):
        """Returns True if a gene is already in the database."""
        cursor = self.conn.cursor()
        cursor.execute("SELECT id FROM genes WHERE symbol = ?", (symbol,))
        return cursor.fetchone() is not None

    def sequence_exists(self, symbol, seq_type, source):
        """
        Returns True if a specific sequence (type + source) already
        exists for a gene. Used to fill in missing sequences without
        re-fetching everything.
        """
        cursor = self.conn.cursor()
        cursor.execute('''
            SELECT s.id FROM sequences s
            JOIN genes g ON s.gene_id = g.id
            WHERE g.symbol = ? AND s.type = ? AND s.source = ?
        ''', (symbol, seq_type, source))
        return cursor.fetchone() is not None

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
        print(f"------------------------\n")