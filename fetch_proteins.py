from Bio import Entrez

Entrez.email = "your_email@example.com"

genes = [
    "TP53","BRCA1","KRAS","EGFR","HBB","GAPDH","MDM2","HRAS","MYC","AKT1",
    "CDK2","CDK4","MTOR","PTEN","RB1","VEGFA","INS","ALB","TNF","IL6"
]

all_sequences = []

for gene in genes:
    handle = Entrez.esearch(
        db="protein",
        term=f"{gene}[Gene] AND Homo sapiens[Organism] AND refseq[filter]",
        retmax=25   # 20 genes × 25 = ~500 proteins
    )
    record = Entrez.read(handle)

    for protein_id in record["IdList"]:
        fetch = Entrez.efetch(
            db="protein",
            id=protein_id,
            rettype="fasta",
            retmode="text"
        )
        seq = fetch.read()
        all_sequences.append(seq)

# Save to file
with open("proteins.fasta", "w") as f:
    f.write("\n".join(all_sequences))
