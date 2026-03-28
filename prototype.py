from Bio import SeqIO
from Bio.Align import PairwiseAligner, substitution_matrices
import os

# 1. LOAD DATA (Keep your FASTA loading logic)
fasta_path = "genes.fasta"
if not os.path.exists(fasta_path):
    print(f"Error: {fasta_path} not found. Ensure it's in the same folder.")
    exit()

# Using dictionary comprehension to load genes efficiently
genes = {record.id: str(record.seq) for record in SeqIO.parse(fasta_path, "fasta")}

# 2. THE UPGRADE: Setup Professional Scoring
# We use BLOSUM62 - the industry standard for biological similarity
matrix = substitution_matrices.load("BLOSUM62")

aligner = PairwiseAligner()
aligner.substitution_matrix = matrix
aligner.mode = 'local'  # Local alignment (Smith-Waterman) is better for different length genes

# Biological penalties: Nature hates gaps
aligner.open_gap_score = -10
aligner.extend_gap_score = -0.5

# 3. ANALYSIS LOGIC
query_id = "BRCA1"
if query_id not in genes:
    print(f"Error: {query_id} not found in {fasta_path}")
    exit()

query_seq = genes[query_id]
results = []

print(f"--- Running Tactical Bio-Alignment for {query_id} ---")

for target_id, target_seq in genes.items():
    # Calculate the raw biological score
    score = aligner.score(query_seq, target_seq)
    
    # Normalize against a perfect match (query vs itself)
    max_score = aligner.score(query_seq, query_seq)
    similarity = max(0, score / max_score)
    results.append((target_id, similarity))

# 4. OUTPUT (Sorted by highest similarity)
results.sort(key=lambda x: x[1], reverse=True)
for gene, score in results:
    print(f"{gene}: {score:.4f}")