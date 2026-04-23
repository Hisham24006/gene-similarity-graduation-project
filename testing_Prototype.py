import os
import ssl
import sqlite3
import requests
from io import StringIO
from Bio import Entrez, SeqIO
from Bio.Align import PairwiseAligner, substitution_matrices

# Force bypass SSL for NCBI/UniProt connection issues
if (not os.environ.get('PYTHONHTTPSVERIFY', '') and getattr(ssl, '_create_unverified_context', None)):
    ssl._create_default_https_context = ssl._create_unverified_context

# ==========================================
# 1. DATABASE MANAGER CLASS
# ==========================================
class GeneDatabase:
    def __init__(self, db_path="gene_vault.db"):
        self.conn = sqlite3.connect(db_path)
        self._create_tables()

    def _create_tables(self):
        cursor = self.conn.cursor()
        # Table for general gene info
        cursor.execute('''CREATE TABLE IF NOT EXISTS genes 
                          (id INTEGER PRIMARY KEY, symbol TEXT UNIQUE, organism TEXT)''')
        # Table for DNA and Protein sequences from various sources
        cursor.execute('''CREATE TABLE IF NOT EXISTS sequences 
                          (gene_id INTEGER, type TEXT, source TEXT, sequence TEXT, 
                           FOREIGN KEY(gene_id) REFERENCES genes(id))''')
        self.conn.commit()

    def save_gene(self, symbol, organism, representations):
        cursor = self.conn.cursor()
        cursor.execute("INSERT OR IGNORE INTO genes (symbol, organism) VALUES (?, ?)", (symbol, organism))
        cursor.execute("SELECT id FROM genes WHERE symbol = ?", (symbol,))
        gene_id = cursor.fetchone()[0]
        
        for rep in representations:
            # Avoid duplicate sequences for the same source/type
            cursor.execute('''INSERT INTO sequences (gene_id, type, source, sequence) 
                              VALUES (?, ?, ?, ?)''', 
                           (gene_id, rep['type'], rep['source'], rep['sequence']))
        self.conn.commit()

    def get_all_local_sequences(self, seq_type):
        """Retrieves all stored sequences of a certain type for comparison."""
        cursor = self.conn.cursor()
        cursor.execute('''SELECT g.symbol, s.sequence FROM sequences s 
                          JOIN genes g ON s.gene_id = g.id 
                          WHERE s.type = ?''', (seq_type,))
        return cursor.fetchall()

# ==========================================
# 2. MULTI-SOURCE FETCHING CLASS
# ==========================================
class BioFetcher:
    def __init__(self, email):
        Entrez.email = email

    def fetch_ncbi(self, symbol, db_type="protein"):
        """Fetches from NCBI (db='protein' or 'nucleotide')."""
        print(f"Searching NCBI {db_type} for {symbol}...")
        try:
            term = f"{symbol}[Gene Name] AND \"Homo sapiens\"[Organism] AND srcdb_refseq[PROP]"
            handle = Entrez.esearch(db=db_type, term=term, retmax=1)
            record = Entrez.read(handle)
            if record["IdList"]:
                fetch = Entrez.efetch(db=db_type, id=record["IdList"][0], rettype="fasta", retmode="text")
                seq_record = SeqIO.read(StringIO(fetch.read()), "fasta")
                return str(seq_record.seq)
        except Exception as e:
            print(f"NCBI Error for {symbol}: {e}")
        return None

    def fetch_uniprot(self, symbol):
        """Fetches Protein sequence from UniProt REST API."""
        print(f"Searching UniProt for {symbol}...")
        url = f"https://rest.uniprot.org/uniprotkb/search?query=gene:{symbol}+AND+organism_id:9606&format=fasta"
        try:
            response = requests.get(url)
            if response.status_code == 200 and response.text:
                seq_record = SeqIO.read(StringIO(response.text), "fasta")
                return str(seq_record.seq)
        except Exception as e:
            print(f"UniProt Error for {symbol}: {e}")
        return None

# ==========================================
# 3. SIMILARITY ENGINE CLASS
# ==========================================
class SimilarityEngine:
    def __init__(self):
        self.aligner = PairwiseAligner()
        self.aligner.mode = 'local'

    def configure(self, mode):
        if mode == "protein":
            self.aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
            self.aligner.open_gap_score = -10
            self.aligner.extend_gap_score = -0.5
        else: # DNA Mode
            self.aligner.match_score = 2
            self.aligner.mismatch_score = -3
            self.aligner.open_gap_score = -5

    def calculate_similarity(self, query_seq, target_seq):
        score = self.aligner.score(query_seq, target_seq)
        max_possible = self.aligner.score(query_seq, query_seq)
        return max(0, score / max_possible)

# ==========================================
# 4. MAIN EXECUTION LOOP
# ==========================================
if __name__ == "__main__":
    # Setup
    user_email = "hishamalsaadi06@gmail.com" # From your prototype.py
    db = GeneDatabase()
    fetcher = BioFetcher(user_email)
    engine = SimilarityEngine()

    # User Inputs
    target_symbol = input("Enter Gene Symbol (e.g., TP53): ").upper()
    search_mode = input("Compare by 'protein' or 'dna'? ").lower()
    
    engine.configure(search_mode)

    # 1. Data Acquisition
    # Check if we have this gene locally first
    local_results = db.get_all_local_sequences(search_mode)
    query_seq = next((seq for sym, seq in local_results if sym == target_symbol), None)

    if not query_seq:
        print(f"\n--- Initializing Multi-Source Fetch for {target_symbol} ---")
        p_ncbi = fetcher.fetch_ncbi(target_symbol, "protein")
        p_uniprot = fetcher.fetch_uniprot(target_symbol)
        d_ncbi = fetcher.fetch_ncbi(target_symbol, "nucleotide")

        # Save representations to SQLite
        reps = []
        if p_ncbi: reps.append({'type': 'protein', 'source': 'NCBI', 'sequence': p_ncbi})
        if p_uniprot: reps.append({'type': 'protein', 'source': 'UniProt', 'sequence': p_uniprot})
        if d_ncbi: reps.append({'type': 'dna', 'source': 'NCBI', 'sequence': d_ncbi})
        
        db.save_gene(target_symbol, "Homo sapiens", reps)
        
        # Set query sequence for this session
        query_seq = p_ncbi if search_mode == "protein" else d_ncbi

    # 2. Run Similarity Search against all genes in Local Database
    if query_seq:
        print(f"\n--- Running {search_mode.upper()} Similarity Search ---")
        all_stored = db.get_all_local_sequences(search_mode)
        results = []

        for symbol, seq in all_stored:
            sim_score = engine.calculate_similarity(query_seq, seq)
            results.append((symbol, sim_score))

        # Sort and Display
        results.sort(key=lambda x: x[1], reverse=True)
        print("\nResults (Sorted by Similarity):")
        for symbol, score in results:
            print(f"{symbol}: {score:.4f}")
    else:
        print("Error: Could not retrieve sequence for the target gene.")
