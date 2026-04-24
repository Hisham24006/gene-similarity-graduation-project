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
from metrics import kmer_similarity

# --- Configuration ---
Entrez.email = "hishamalsaadi06@gmail.com"
DB_PATH = os.path.join(os.path.dirname(__file__), '..', 'database', 'gene_vault.db')
KMER_K = 3  # k-mer size for proteins

# --- BLOSUM62 Aligner setup ---
aligner = PairwiseAligner()
aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
aligner.mode = 'local'
aligner.open_gap_score = -10
aligner.extend_gap_score = -0.5


def fetch_protein_from_ncbi(gene_name):
    """
    Fetches the top RefSeq protein for a gene from NCBI.
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


def blosum_similarity(query_seq, target_seq):
    """Normalized BLOSUM62 local alignment similarity score (0 to 1)."""
    score = aligner.score(query_seq, target_seq)
    max_score = aligner.score(query_seq, query_seq)
    if max_score == 0:
        return 0.0
    return max(0.0, score / max_score)


# --- Main ---
db = GeneDatabase(db_path=DB_PATH)

print("Genes in database:", ", ".join(db.list_genes()))

target_gene = input("\nEnter a gene name to search (e.g., TP53, BRCA2, KRAS): ").upper()
query_isoform_id, query_seq = fetch_protein_from_ncbi(target_gene)

if query_seq:
    all_seqs = db.get_all_sequences_with_isoform("protein")

    print(f"\n--- Running similarity metrics for {target_gene} ({query_isoform_id}) ---")
    print(f"    Database: {len(all_seqs)} protein isoforms")
    print(f"    Metrics : BLOSUM62 alignment, K-mer (k={KMER_K})\n")

    results = []
    for symbol, isoform_id, local_seq in all_seqs:
        blosum_score = blosum_similarity(query_seq, local_seq)
        kmer_score   = kmer_similarity(query_seq, local_seq, k=KMER_K)
        avg_score    = (blosum_score + kmer_score) / 2
        results.append((symbol, isoform_id, blosum_score, kmer_score, avg_score))

    # Sort by average score
    results.sort(key=lambda x: x[4], reverse=True)

    # Print combined results table
    print(f"{'Gene':<10} {'Isoform':<28} {'BLOSUM62':>10} {'K-mer':>8} {'Average':>9}")
    print(f"{'-'*10} {'-'*28} {'-'*10} {'-'*8} {'-'*9}")
    for symbol, isoform_id, blosum, kmer, avg in results[:15]:
        print(f"{symbol:<10} {isoform_id:<28} {blosum:>10.4f} {kmer:>8.4f} {avg:>9.4f}")

    # Best match section
    top_symbol, top_isoform_id, top_blosum, top_kmer, top_avg = results[0]
    print(f"\n--- Best match ---")
    print(f"  Query   : {target_gene} | {query_isoform_id}")
    print(f"  Match   : {top_symbol} | {top_isoform_id}")
    print(f"  BLOSUM62: {top_blosum:.4f}")
    print(f"  K-mer   : {top_kmer:.4f}")
    print(f"  Average : {top_avg:.4f}")

    # Show alignment for best non-self match
    non_self = [(s, i, b, k, a) for s, i, b, k, a in results if s != target_gene]
    if non_self:
        best_symbol, best_isoform, best_blosum, best_kmer, best_avg = non_self[0]
        print(f"\n--- Best match outside {target_gene} ---")
        print(f"  Gene    : {best_symbol} | {best_isoform}")
        print(f"  BLOSUM62: {best_blosum:.4f}")
        print(f"  K-mer   : {best_kmer:.4f}")
        print(f"\n--- Alignment ---")
        best_seq = db.get_isoform_sequence(best_symbol, "protein", best_isoform)
        print(aligner.align(query_seq, best_seq)[0])

db.close()
