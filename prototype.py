import ssl
import os

# This tells the system to ignore SSL certificate verification
if (not os.environ.get('PYTHONHTTPSVERIFY', '') and getattr(ssl, '_create_unverified_context', None)):
    ssl._create_default_https_context = ssl._create_unverified_context
from Bio import SeqIO
from Bio import Entrez
from Bio.Align import PairwiseAligner, substitution_matrices
import os
from io import StringIO

# --- CONFIGURATION (Protein Mode) ---
Entrez.email = "hishamalsaadi06@gmail.com"
fasta_path = "proteins.fasta" # Rename your local DB file to reflect Proteins

# Setup the Protein Engine (BLOSUM62)
aligner = PairwiseAligner()
aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
aligner.mode = 'local'
aligner.open_gap_score = -10
aligner.extend_gap_score = -0.5

# --- FUNCTIONS ---
# --- FUNCTIONS ---
def fetch_protein_from_ncbi(gene_name):
    """Scouts NCBI for a human PROTEIN sequence."""
    print(f"\n--- Scouting NCBI Protein Database for: {gene_name} ---")
    try:
        # DB changed to 'protein'
        handle = Entrez.esearch(db="protein", term=f"{gene_name}[Gene Name] AND human[Organism] AND srcdb_refseq[PROP]", retmax=1)
        record = Entrez.read(handle)
        handle.close()

        if not record["IdList"]:
            print(f"No protein results found for {gene_name}")
            return None

        protein_id = record["IdList"][0]
        handle = Entrez.efetch(db="protein", id=protein_id, rettype="fasta", retmode="text")
        fasta_data = handle.read()
        handle.close()

        seq_record = SeqIO.read(StringIO(fasta_data), "fasta")
        print(f"Successfully retrieved Protein {gene_name} (Length: {len(seq_record.seq)} AA)")
        return str(seq_record.seq)
    except Exception as e:
        print(f"Failed to fetch Protein {gene_name}: {e}")
        return None

# --- EXECUTION ---
if not os.path.exists(fasta_path):
    print(f"Error: {fasta_path} not found. Ensure it's in the same folder.")
    exit()

# Load local database
gene_db = {record.id: str(record.seq) for record in SeqIO.parse(fasta_path, "fasta")}

# 1. Ask for input
target_gene = input("Enter a gene name to fetch from NCBI (e.g., TP53, BRCA2, KRAS): ")
query_seq = fetch_protein_from_ncbi(target_gene)

if query_seq:
    results = []
    print(f"\n--- Comparing {target_gene} against local database ---")
    
    # 2. Run Biological Alignment Loop
    for gene_id, local_seq in gene_db.items():
        score = aligner.score(query_seq, local_seq)
        # Normalize against perfect match (query vs itself)
        max_score = aligner.score(query_seq, query_seq)
        similarity = max(0, score / max_score)
        results.append((gene_id, similarity))

    # 3. Sort and Print Similarity Scores
    results.sort(key=lambda x: x[1], reverse=True)
    for gene, score in results:
        print(f"{gene}: {score:.4f}")

    # 4. Visualization of top match
    if len(results) > 0:
        top_id = results[0][0]
        top_similarity = results[0][1]
        print(f"\n--- Best Local Alignment: {target_gene} vs {top_id} ---")
        print(f"Similarity Score: {top_similarity:.4f}")
        
        alignments = aligner.align(query_seq, gene_db[top_id])
        print(alignments[0])