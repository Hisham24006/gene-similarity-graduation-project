import ssl
import os
import sys

if not os.environ.get('PYTHONHTTPSVERIFY', '') and getattr(ssl, '_create_unverified_context', None):
    ssl._create_default_https_context = ssl._create_unverified_context

# Add database/ folder to path so we can import GeneDatabase
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'database'))

from Bio import SeqIO, Entrez
from Bio.Align import PairwiseAligner, substitution_matrices
from io import StringIO
from database_manager import GeneDatabase

# --- Configuration ---
Entrez.email = "hishamalsaadi06@gmail.com"
DB_PATH = os.path.join(os.path.dirname(__file__), '..', 'database', 'gene_vault.db')

# --- Aligner setup (BLOSUM62) ---
aligner = PairwiseAligner()
aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
aligner.mode = 'local'
aligner.open_gap_score = -10
aligner.extend_gap_score = -0.5


def fetch_protein_from_ncbi(gene_name):
    """Fetches the top RefSeq protein for a gene from NCBI.
    Returns (isoform_id, sequence) or (None, None).
    """
    print(f"\n--- Fetching {gene_name} from NCBI ---")
    try:
        handle = Entrez.esearch(
            db="protein",
            term=f"{gene_name}[Gene Name] AND human[Organism] AND srcdb_refseq[PROP]",
            retmax=1
        )
        record = Entrez.read(handle)
        handle.close()

        if not record["IdList"]:
            print(f"No results found for {gene_name}")
            return None, None

        handle = Entrez.efetch(
            db="protein", id=record["IdList"][0],
            rettype="fasta", retmode="text"
        )
        seq_record = SeqIO.read(StringIO(handle.read()), "fasta")
        handle.close()

        print(f"Retrieved isoform : {seq_record.id}")
        print(f"Length            : {len(seq_record.seq)} AA")
        return seq_record.id, str(seq_record.seq)

    except Exception as e:
        print(f"Failed to fetch {gene_name}: {e}")
        return None, None


# --- Main ---
db = GeneDatabase(db_path=DB_PATH)

print("Genes in database:", ", ".join(db.list_genes()))

target_gene = input("\nEnter a gene name to search (e.g., TP53, BRCA2, KRAS): ").upper()
query_isoform_id, query_seq = fetch_protein_from_ncbi(target_gene)

if query_seq:
    # get_all_sequences now also returns isoform_id
    all_seqs = db.get_all_sequences_with_isoform("protein")
    results = []

    print(f"\n--- Comparing {target_gene} ({query_isoform_id}) against {len(all_seqs)} isoforms in database ---")

    for symbol, isoform_id, local_seq in all_seqs:
        score = aligner.score(query_seq, local_seq)
        max_score = aligner.score(query_seq, query_seq)
        similarity = max(0, score / max_score)
        results.append((symbol, isoform_id, similarity))

    results.sort(key=lambda x: x[2], reverse=True)

    print("\nTop 10 most similar isoforms:")
    print(f"  {'Gene':<10} {'Isoform':<25} {'Score'}")
    print(f"  {'-'*10} {'-'*25} {'-'*6}")
    for symbol, isoform_id, score in results[:10]:
        print(f"  {symbol:<10} {isoform_id:<25} {score:.4f}")

    top_symbol, top_isoform_id, top_sim = results[0]
    print(f"\n--- Best match ---")
    print(f"  Query   : {target_gene} | {query_isoform_id}")
    print(f"  Match   : {top_symbol} | {top_isoform_id}")
    print(f"  Score   : {top_sim:.4f}")
    print(f"\n--- Alignment ---")
    top_seq = db.get_isoform_sequence(top_symbol, "protein", top_isoform_id)
    print(aligner.align(query_seq, top_seq)[0])

db.close()