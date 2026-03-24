from Bio import SeqIO
import Levenshtein

# Load genes
gene_db = {}
for record in SeqIO.parse("genes.fasta", "fasta"):
    gene_db[record.id] = str(record.seq)

print("Loaded genes:", gene_db.keys())


# Similarity function
def similarity_score(seq1, seq2):
    distance = Levenshtein.distance(seq1, seq2)
    return 1 - distance / max(len(seq1), len(seq2))


# Query gene (you can change this)
query_gene = gene_db["BRCA1"]


# Compare with all genes
results = []
for gene_id, seq in gene_db.items():
    score = similarity_score(query_gene, seq)
    results.append((gene_id, score))


# Sort results
results.sort(key=lambda x: x[1], reverse=True)


# Print top matches
print("\nTop similar genes to BRCA1:\n")
for gene, score in results:
    print(f"{gene}: {score:.2f}")