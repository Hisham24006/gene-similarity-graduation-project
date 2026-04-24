import ssl
import os
import sys

if not os.environ.get('PYTHONHTTPSVERIFY', '') and getattr(ssl, '_create_unverified_context', None)):
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

# Valid amino acid characters
VALID_AA = set("ACDEFGHIKLMNPQRSTVWYacdefghiklmnpqrstvwy")


def is_protein_sequence(text):
    """
    Returns True if the input looks like a raw protein sequence
    (only amino acid letters, no spaces or special characters).
    """
    cleaned = text.replace("\n", "").replace(" ", "").replace("\r", "")
    return len(cleaned) > 10 and all(c in VALID_AA for c in cleaned)


def clean_sequence(text):
    """Strips whitespace and newlines from a pasted sequence."""
    return text.replace("\n", "").replace(" ", "").replace("\r", "").upper()


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


def get_query():
    """
    Asks the user how they want to input their query.
    Returns (label, sequence) where label is either the gene name or 'custom_sequence'.
    """
    print("\n========================================")
    print("  Gene Similarity Search")
    print("========================================")
    print("How would you like to input your query?")
    print("  [1] Search by gene name (fetched from NCBI)")
    print("  [2] Paste a raw protein sequence")
    print("========================================")

    choice = input("Enter 1 or 2: ").strip()

    if choice == "1":
        gene_name = input("Enter gene name (e.g., TP53, BRCA2, KRAS): ").upper().strip()
        isoform_id, seq = fetch_protein_from_ncbi(gene_name)
        if seq:
            return gene_name, isoform_id, seq
        else:
            print("Could not retrieve sequence. Exiting.")
            return None, None, None

    elif choice == "2":
        print("\nPaste your protein sequence below.")
        print("You can paste multiple lines — when done, enter a blank line:")
        lines = []
        while True:
            line = input()
            if line.strip() == "":
                break
            lines.append(line.strip())

        raw = "".join(lines)

        # Handle FASTA format (starts with >)
        if raw.startswith(">"):
            fasta_lines = raw.split("\n")
            header = fasta_lines[0]
            seq = clean_sequence("".join(fasta_lines[1:]))
            label = header[1:].split()[0]  # use first word after > as label
            print(f"\nDetected FASTA format.")
            print(f"Label  : {label}")
            print(f"Length : {len(seq)} AA")
            return label, label, seq

        # Plain sequence
        seq = clean_sequence(raw)
        if is_protein_sequence(seq):
            print(f"\nSequence accepted.")
            print(f"Length : {len(seq)} AA")
            return "custom_sequence", "custom_sequence", seq
        else:
            print("Invalid sequence — contains non-amino acid characters.")
            return None, None, None

    else:
        print("Invalid choice.")
        return None, None, None


# --- Main ---
db = GeneDatabase(db_path=DB_PATH)

print("Genes in database:", ", ".join(db.list_genes()))

target_gene, query_isoform_id, query_seq = get_query()

if query_seq:
    all_seqs = db.get_all_sequences_with_isoform("protein")

    print(f"\n--- Running similarity metrics for {target_gene} ({query_isoform_id}) ---")
    print(f"    Database : {len(all_seqs)} protein isoforms")
    print(f"    Metrics  : BLOSUM62 alignment, K-mer (k={KMER_K})\n")

    results = []
    for symbol, isoform_id, local_seq in all_seqs:
        blosum_score = blosum_similarity(query_seq, local_seq)
        kmer_score   = kmer_similarity(query_seq, local_seq, k=KMER_K)
        avg_score    = (blosum_score + kmer_score) / 2
        results.append((symbol, isoform_id, blosum_score, kmer_score, avg_score))

    results.sort(key=lambda x: x[4], reverse=True)

    # Print combined results table
    print(f"{'Gene':<10} {'Isoform':<28} {'BLOSUM62':>10} {'K-mer':>8} {'Average':>9}")
    print(f"{'-'*10} {'-'*28} {'-'*10} {'-'*8} {'-'*9}")
    for symbol, isoform_id, blosum, kmer, avg in results[:15]:
        print(f"{symbol:<10} {isoform_id:<28} {blosum:>10.4f} {kmer:>8.4f} {avg:>9.4f}")

    # Best match
    top_symbol, top_isoform_id, top_blosum, top_kmer, top_avg = results[0]
    print(f"\n--- Best match ---")
    print(f"  Query   : {target_gene} | {query_isoform_id}")
    print(f"  Match   : {top_symbol} | {top_isoform_id}")
    print(f"  BLOSUM62: {top_blosum:.4f}")
    print(f"  K-mer   : {top_kmer:.4f}")
    print(f"  Average : {top_avg:.4f}")

    # Best non-self match + alignment
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