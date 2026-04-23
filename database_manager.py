import sqlite3

class GeneDatabase:
    def __init__(self, db_path="gene_vault.db"):
        self.conn = sqlite3.connect(db_path)
        self._create_tables()

    def _create_tables(self):
        cursor = self.conn.cursor()
        # Table for general gene info
        cursor.execute('''CREATE TABLE IF NOT EXISTS genes 
                          (id INTEGER PRIMARY KEY, symbol TEXT UNIQUE, organism TEXT)''')
        # Table for different sequence representations (DNA/Protein) from different sources
        cursor.execute('''CREATE TABLE IF NOT EXISTS sequences 
                          (gene_id INTEGER, type TEXT, source TEXT, sequence TEXT, 
                           FOREIGN KEY(gene_id) REFERENCES genes(id))''')
        self.conn.commit()

    def save_gene_data(self, symbol, organism, representations):
        cursor = self.conn.cursor()
        try:
            cursor.execute("INSERT OR IGNORE INTO genes (symbol, organism) VALUES (?, ?)", (symbol, organism))
            cursor.execute("SELECT id FROM genes WHERE symbol = ?", (symbol,))
            gene_id = cursor.fetchone()[0]
            
            for rep in representations:
                cursor.execute("INSERT INTO sequences (gene_id, type, source, sequence) VALUES (?, ?, ?, ?)",
                               (gene_id, rep['type'], rep['source'], rep['sequence']))
            self.conn.commit()
        except Exception as e:
            print(f"Database Error: {e}")

    def get_sequences(self, symbol, seq_type):
        cursor = self.conn.cursor()
        cursor.execute('''SELECT s.sequence, s.source FROM sequences s 
                          JOIN genes g ON s.gene_id = g.id 
                          WHERE g.symbol = ? AND s.type = ?''', (symbol, seq_type))
        return cursor.fetchall()
