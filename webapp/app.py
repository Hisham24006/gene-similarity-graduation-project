"""
app.py — Gene Similarity Search Flask App
"""

import os
import sys
import ssl
import glob
import json
import base64
import time
import random

import numpy as np
from flask import Flask, render_template, request, jsonify

# ── Paths ──
BASE_DIR    = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SIM_DIR     = os.path.join(BASE_DIR, 'similarity')
DB_DIR      = os.path.join(BASE_DIR, 'database')
RESULTS_DIR     = os.path.join(BASE_DIR, 'results')
RECENT_SEARCHES = os.path.join(BASE_DIR, 'data', 'recent_searches.json')
sys.path.insert(0, SIM_DIR)
sys.path.insert(0, DB_DIR)

if not os.environ.get('PYTHONHTTPSVERIFY', '') and getattr(ssl, '_create_unverified_context', None):
    ssl._create_default_https_context = ssl._create_unverified_context

from database_manager import GeneDatabase
from metrics import kmer_similarity, edit_distance_similarity
from motif_similarity import load_jaspar_motifs, fetch_promoter, cached_scan_motifs, jaccard_motif_similarity
from visualizer import plot_bar_chart, plot_heatmap
from Bio import SeqIO, Entrez
from Bio.Align import PairwiseAligner, substitution_matrices
from io import StringIO

# ── Config ──
Entrez.email = "hishamalsaadi06@gmail.com"
DB_PATH      = os.path.join(DB_DIR, 'gene_vault.db')
KMER_K       = 3

app = Flask(__name__)

# ── Gene descriptions (general, not isoform-specific) ──
GENE_DESCRIPTIONS = {
    "TP53":   "Tumor protein p53 — master regulator of the cell cycle and tumor suppressor. Activates DNA repair and apoptosis in response to stress.",
    "BRCA1":  "Breast cancer type 1 susceptibility protein — involved in DNA repair and maintenance of genomic stability. Mutations greatly increase breast/ovarian cancer risk.",
    "BRCA2":  "Breast cancer type 2 susceptibility protein — essential for homologous recombination DNA repair. Mutations increase susceptibility to breast, ovarian, and other cancers.",
    "KRAS":   "KRAS proto-oncogene — GTPase that regulates cell proliferation and differentiation. One of the most commonly mutated oncogenes in human cancers.",
    "EGFR":   "Epidermal growth factor receptor — receptor tyrosine kinase that drives cell growth. Frequently overexpressed or mutated in lung, colorectal, and other cancers.",
    "HBB":    "Hemoglobin subunit beta — forms the oxygen-carrying hemoglobin protein in red blood cells. Mutations cause sickle cell disease and beta-thalassemia.",
    "GAPDH":  "Glyceraldehyde-3-phosphate dehydrogenase — key enzyme in glycolysis. Widely used as a housekeeping gene and loading control in molecular biology.",
    "MDM2":   "MDM2 proto-oncogene — E3 ubiquitin ligase that negatively regulates TP53 by targeting it for degradation. Amplified in many cancers.",
    "HRAS":   "Harvey rat sarcoma viral proto-oncogene — GTPase involved in signal transduction regulating cell growth and differentiation.",
    "MYC":    "MYC proto-oncogene — transcription factor that regulates expression of many genes controlling cell growth, proliferation, and apoptosis. A major oncogene.",
    "AKT1":   "AKT serine/threonine kinase 1 — key mediator of the PI3K signaling pathway, regulating cell survival, growth, and metabolism.",
    "CDK2":   "Cyclin-dependent kinase 2 — regulates cell cycle progression through S phase and G2/M transition.",
    "CDK4":   "Cyclin-dependent kinase 4 — drives G1/S cell cycle transition by phosphorylating RB1. Frequently dysregulated in cancer.",
    "MTOR":   "Mechanistic target of rapamycin — serine/threonine kinase that integrates nutrient and growth factor signals to regulate cell growth and metabolism.",
    "PTEN":   "Phosphatase and tensin homolog — tumor suppressor that negatively regulates the PI3K/AKT pathway. One of the most frequently mutated genes in cancer.",
    "RB1":    "Retinoblastoma protein — tumor suppressor that controls cell cycle entry. Loss of RB1 is a founding event in retinoblastoma and many other cancers.",
    "VEGFA":  "Vascular endothelial growth factor A — major regulator of angiogenesis (new blood vessel formation). Overexpressed in many tumors.",
    "INS":    "Insulin — pancreatic hormone that regulates blood glucose by promoting glucose uptake. Deficiency leads to diabetes mellitus.",
    "ALB":    "Serum albumin — most abundant plasma protein, responsible for maintaining oncotic pressure and transporting hormones, fatty acids, and drugs.",
    "TNF":    "Tumor necrosis factor — pro-inflammatory cytokine involved in systemic inflammation, immune response, and apoptosis.",
    "IL6":    "Interleukin 6 — multifunctional cytokine involved in immune response, inflammation, and hematopoiesis.",
    "BRCA2":  "Breast cancer type 2 susceptibility protein — essential for homologous recombination DNA repair.",
    "ATM":    "ATM serine/threonine kinase — master regulator of DNA double-strand break repair and cell cycle checkpoint signaling.",
    "BRAF":   "BRAF proto-oncogene — serine/threonine kinase in the RAS/MAPK signaling pathway. V600E mutation is common in melanoma and other cancers.",
    "PIK3CA": "Phosphatidylinositol-4,5-bisphosphate 3-kinase catalytic subunit alpha — activates the PI3K/AKT/mTOR pathway. Frequently mutated in cancer.",
    "ERBB2":  "Erb-b2 receptor tyrosine kinase 2 (HER2) — receptor that drives cell growth. Amplified in ~25% of breast cancers; target of trastuzumab.",
    "NRAS":   "NRAS proto-oncogene — RAS family GTPase regulating cell proliferation. Mutations common in melanoma, leukemia, and thyroid cancer.",
    "ABL1":   "ABL proto-oncogene 1 — non-receptor tyrosine kinase involved in cell differentiation and division. BCR-ABL fusion drives chronic myeloid leukemia.",
    "BCL2":   "BCL2 apoptosis regulator — anti-apoptotic protein that promotes cell survival by preventing mitochondrial cytochrome c release.",
    "CDKN2A": "Cyclin-dependent kinase inhibitor 2A — tumor suppressor encoding p16(INK4a) and p14(ARF), both inhibiting cell cycle progression.",
    "FGFR1":  "Fibroblast growth factor receptor 1 — receptor tyrosine kinase involved in cell proliferation, differentiation, and angiogenesis.",
    "SRC":    "SRC proto-oncogene — non-receptor tyrosine kinase involved in multiple signaling pathways controlling cell growth and survival.",
    "JAK2":   "Janus kinase 2 — tyrosine kinase critical for cytokine signaling. JAK2 V617F mutation is the hallmark of polycythemia vera.",
    "STAT3":  "Signal transducer and activator of transcription 3 — transcription factor activated by cytokines, frequently constitutively active in cancers.",
    "RAF1":   "RAF1 proto-oncogene — serine/threonine kinase in the RAS/MAPK signaling cascade regulating cell proliferation and survival.",
    "MAP2K1": "Mitogen-activated protein kinase kinase 1 (MEK1) — dual-specificity kinase in the MAPK/ERK signaling pathway.",
    "MAP2K2": "Mitogen-activated protein kinase kinase 2 (MEK2) — dual-specificity kinase that phosphorylates and activates ERK1/2.",
    "MAPK1":  "Mitogen-activated protein kinase 1 (ERK2) — serine/threonine kinase mediating cellular responses to growth signals.",
    "MAPK3":  "Mitogen-activated protein kinase 3 (ERK1) — serine/threonine kinase involved in cell growth and differentiation.",
    "AKT2":   "AKT serine/threonine kinase 2 — PI3K effector important for insulin signaling and glucose metabolism.",
    "AKT3":   "AKT serine/threonine kinase 3 — PI3K effector involved in neuronal survival and brain development.",
    "MYOD1":  "Myogenic differentiation 1 — master transcription factor for skeletal muscle differentiation.",
    "E2F1":   "E2F transcription factor 1 — regulates genes needed for cell cycle progression from G1 to S phase.",
    "SP1":    "Specificity protein 1 — constitutive transcription factor that activates or represses expression of many genes.",
    "FOXP3":  "Forkhead box P3 — transcription factor essential for regulatory T cell development and immune tolerance.",
    "GATA1":  "GATA binding protein 1 — transcription factor essential for erythroid and megakaryocyte development.",
    "TP63":   "Tumor protein p63 — TP53 family member critical for epithelial development and stem cell maintenance.",
    "TP73":   "Tumor protein p73 — TP53 family member involved in neuronal development and apoptosis.",
    "PCNA":   "Proliferating cell nuclear antigen — DNA clamp protein essential for DNA replication and repair.",
    "CCND1":  "Cyclin D1 — regulatory subunit of CDK4/6 that drives G1/S cell cycle transition. Frequently overexpressed in cancer.",
    "CCNE1":  "Cyclin E1 — regulatory subunit of CDK2, required for G1/S transition. Amplified in ovarian and other cancers.",
    "MLH1":   "MutL homolog 1 — DNA mismatch repair protein. Mutations cause Lynch syndrome and microsatellite instability in colorectal cancer.",
    "MSH2":   "MutS homolog 2 — DNA mismatch repair protein. Mutations cause Lynch syndrome predisposing to colorectal and other cancers.",
    "RAD51":  "RAD51 recombinase — essential for homologous recombination repair of DNA double-strand breaks.",
    "XRCC1":  "X-ray repair cross complementing 1 — scaffold protein coordinating base excision repair and single-strand break repair.",
    "ERCC1":  "Excision repair cross-complementation group 1 — essential for nucleotide excision repair of bulky DNA lesions.",
    "CHEK1":  "Checkpoint kinase 1 — serine/threonine kinase that mediates cell cycle arrest in response to DNA damage.",
    "CHEK2":  "Checkpoint kinase 2 — tumor suppressor kinase activated by ATM in response to DNA double-strand breaks.",
    "FANCA":  "Fanconi anemia complementation group A — involved in the Fanconi anemia DNA repair pathway. Mutations cause Fanconi anemia.",
    "PALB2":  "Partner and localizer of BRCA2 — bridges BRCA1 and BRCA2 in homologous recombination repair. Mutations increase breast cancer risk.",
    "NBN":    "Nibrin (NBS1) — component of the MRN complex involved in DNA double-strand break sensing and repair.",
}

# ── Startup ──
print("Loading JASPAR motifs...")
BIO_MOTIFS = load_jaspar_motifs()
print(f"Loaded {len(BIO_MOTIFS)} motif models.")

aligner = PairwiseAligner()
aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
aligner.mode = 'local'
aligner.open_gap_score  = -10
aligner.extend_gap_score = -0.5

dna_aligner = PairwiseAligner()
dna_aligner.mode = 'local'
dna_aligner.match_score    = 2
dna_aligner.mismatch_score = -3
dna_aligner.open_gap_score  = -5
dna_aligner.extend_gap_score = -2

VALID_AA  = set("ACDEFGHIKLMNPQRSTVWYacdefghiklmnpqrstvwy")
VALID_DNA = set("ATGCatgcNn")


# ── Helpers ──

def blosum_similarity(seq_a, seq_b):
    score = aligner.score(seq_a, seq_b)
    mx    = aligner.score(seq_a, seq_a)
    return max(0.0, score / mx) if mx > 0 else 0.0


def dna_similarity(seq_a, seq_b):
    score = dna_aligner.score(seq_a, seq_b)
    mx    = dna_aligner.score(seq_a, seq_a)
    return max(0.0, score / mx) if mx > 0 else 0.0


def fetch_from_ncbi(gene_name, seq_type='protein'):
    """Fetch protein or DNA sequence from NCBI. Returns (isoform_id, seq, error)."""
    try:
        if seq_type == 'protein':
            db_name = "protein"
            term    = f"{gene_name}[Gene Name] AND human[Organism] AND srcdb_refseq[PROP]"
        else:
            db_name = "nucleotide"
            term    = f"{gene_name}[Gene Name] AND Homo sapiens[Organism] AND mRNA[Filter] AND srcdb_refseq[PROP]"

        handle = Entrez.esearch(db=db_name, term=term, retmax=1)
        record = Entrez.read(handle)
        handle.close()

        if not record["IdList"]:
            return None, None, f"No {seq_type} results found on NCBI for {gene_name}"

        fmt    = "fasta"
        handle = Entrez.efetch(db=db_name, id=record["IdList"][0], rettype=fmt, retmode="text")
        seq_record = SeqIO.read(StringIO(handle.read()), "fasta")
        handle.close()
        return seq_record.id, str(seq_record.seq), None

    except Exception as e:
        return None, None, str(e)


def get_promoter_tfs(symbol):
    promoter = fetch_promoter(symbol)
    return cached_scan_motifs(symbol, promoter, BIO_MOTIFS) if promoter else set()


def image_to_base64(path):
    if not os.path.exists(path):
        return None
    with open(path, 'rb') as f:
        return base64.b64encode(f.read()).decode('utf-8')


def clean_gene_description(raw_desc, symbol):
    """Return a clean general description for a gene, ignoring isoform-specific text."""
    if symbol.upper() in GENE_DESCRIPTIONS:
        return GENE_DESCRIPTIONS[symbol.upper()]
    if not raw_desc or raw_desc == 'N/A':
        return None
    # Strip isoform suffixes like "isoform a", "isoform 1", "transcript variant X"
    import re
    cleaned = re.sub(r',?\s*(isoform\s*\w+|transcript variant\s*\w+|mRNA\s*\w*)', '', raw_desc, flags=re.IGNORECASE)
    cleaned = re.sub(r'\[.*?\]', '', cleaned).strip().strip(',').strip()
    return cleaned or None


def save_recent_search(query_label, seq_mode, metrics, result_count):
    """Save a search to the recent searches log (max 20 entries)."""
    try:
        searches = []
        if os.path.exists(RECENT_SEARCHES):
            with open(RECENT_SEARCHES, 'r') as f:
                searches = json.load(f)

        entry = {
            'query':    query_label,
            'mode':     seq_mode,
            'metrics':  metrics,
            'results':  result_count,
            'time':     time.strftime('%Y-%m-%d %H:%M'),
        }

        # Remove duplicate query if exists
        searches = [s for s in searches if s['query'] != query_label]
        searches.insert(0, entry)
        searches = searches[:20]  # keep last 20

        os.makedirs(os.path.dirname(RECENT_SEARCHES), exist_ok=True)
        with open(RECENT_SEARCHES, 'w') as f:
            json.dump(searches, f)
    except Exception:
        pass


def load_recent_searches():
    """Load recent searches from file."""
    try:
        if os.path.exists(RECENT_SEARCHES):
            with open(RECENT_SEARCHES, 'r') as f:
                return json.load(f)
    except Exception:
        pass
    return []


def generate_charts(results, query_label, all_seqs, promoter_tfs):
    """Generate bar chart and heatmaps after a search."""
    try:
        os.makedirs(RESULTS_DIR, exist_ok=True)
        plot_bar_chart(results, query_label, output_dir=RESULTS_DIR)
        plot_heatmap(all_seqs, lambda a, b: kmer_similarity(a, b, k=KMER_K),
                     "K-mer (k=3)", query_label, output_dir=RESULTS_DIR)
        plot_heatmap(all_seqs, edit_distance_similarity,
                     "Edit distance", query_label, output_dir=RESULTS_DIR)
        if promoter_tfs:
            gene_tf_map = promoter_tfs
            plot_heatmap(
                all_seqs,
                lambda a, b: jaccard_motif_similarity(
                    gene_tf_map.get(next((s for s, i, sq in all_seqs if sq == a), ''), set()),
                    gene_tf_map.get(next((s for s, i, sq in all_seqs if sq == b), ''), set())
                ),
                "Motif (JASPAR)", query_label, output_dir=RESULTS_DIR
            )
    except Exception as e:
        print(f"Chart generation error: {e}")


# ── Page routes ──

@app.route('/')
def index():
    return render_template('search.html')


@app.route('/database')
def database():
    return render_template('database.html')


@app.route('/visualize')
def visualize():
    charts = {}
    patterns = {
        'bar':       'bar_*.png',
        'kmer':      'heatmap_k-mer_k=3.png',
        'edit':      'heatmap_edit_distance.png',
        'motif':     'heatmap_motif_jaspar.png',
        'speed':     'benchmark_speed.png',
        'agreement': 'benchmark_agreement.png',
    }
    os.makedirs(RESULTS_DIR, exist_ok=True)
    for key, pattern in patterns.items():
        matches = glob.glob(os.path.join(RESULTS_DIR, pattern))
        if matches:
            charts[key] = image_to_base64(sorted(matches)[-1])

    return render_template('visualize.html', charts=charts)


@app.route('/benchmark')
def benchmark():
    return render_template('benchmark.html')


# ── API routes ──

@app.route('/api/search', methods=['POST'])
def api_search():
    data       = request.json
    input_type = data.get('input_type', 'gene_name')
    fetch_type = data.get('fetch_type', 'protein')   # 'protein' or 'dna'
    seq_mode   = data.get('seq_mode', 'protein')      # what to search against
    metrics    = data.get('metrics', ['blosum', 'kmer', 'edit', 'motif'])
    top_n      = int(data.get('top_n', 10))

    # ── Get query sequence ──
    if input_type == 'gene_name':
        gene_name   = data.get('gene_name', '').upper().strip()
        isoform_id, query_seq, error = fetch_from_ncbi(gene_name, fetch_type)
        if error:
            return jsonify({'error': error}), 400
        query_label = gene_name
        actual_mode = fetch_type

    else:
        raw   = data.get('sequence', '').strip()
        lines = [l.strip() for l in raw.splitlines() if l.strip()]
        if lines and lines[0].startswith('>'):
            isoform_id  = lines[0][1:].split()[0]
            query_seq   = ''.join(lines[1:]).upper().replace(' ', '')
            query_label = isoform_id
        else:
            query_seq   = ''.join(lines).upper().replace(' ', '')
            isoform_id  = 'custom'
            query_label = 'custom_sequence'

        # Auto-detect sequence type
        dna_chars = set(query_seq) - set('ATGCN')
        if dna_chars:
            actual_mode = 'protein'
            if not all(c in VALID_AA for c in query_seq):
                return jsonify({'error': 'Invalid sequence — unrecognized characters'}), 400
        else:
            actual_mode = 'dna'

    # ── Load DB sequences ──
    db       = GeneDatabase(db_path=DB_PATH)
    all_seqs = db.get_all_sequences_with_isoform(actual_mode)
    db.close()

    if not all_seqs:
        return jsonify({'error': f'No {actual_mode} sequences found in database.'}), 400

    # ── Motif TF sets ──
    promoter_tfs = {}
    query_tfs    = set()
    if 'motif' in metrics:
        query_tfs    = get_promoter_tfs(query_label) if query_label != 'custom_sequence' else set()
        gene_symbols = list(set(s for s, i, _ in all_seqs))
        for symbol in gene_symbols:
            promoter_tfs[symbol] = get_promoter_tfs(symbol)

    # ── Run metrics ──
    results = []
    for symbol, iso_id, local_seq in all_seqs:
        row   = {'gene': symbol, 'isoform': iso_id, 'scores': {}, 'average': 0.0}
        total = 0; count = 0

        if actual_mode == 'protein':
            if 'blosum' in metrics:
                s = blosum_similarity(query_seq, local_seq)
                row['scores']['blosum'] = round(s, 4); total += s; count += 1
        else:
            if 'blosum' in metrics:
                s = dna_similarity(query_seq, local_seq)
                row['scores']['blosum'] = round(s, 4); total += s; count += 1

        if 'kmer' in metrics:
            s = kmer_similarity(query_seq, local_seq, k=KMER_K)
            row['scores']['kmer'] = round(s, 4); total += s; count += 1

        if 'edit' in metrics:
            s = edit_distance_similarity(query_seq, local_seq)
            row['scores']['edit'] = round(s, 4); total += s; count += 1

        if 'motif' in metrics:
            s = jaccard_motif_similarity(query_tfs, promoter_tfs.get(symbol, set()))
            row['scores']['motif'] = round(s, 4); total += s; count += 1

        row['average'] = round(total / count, 4) if count > 0 else 0.0
        results.append(row)

    results.sort(key=lambda x: x['average'], reverse=True)

    # ── Save recent search ──
    save_recent_search(query_label, actual_mode, metrics, len(results))

    # ── Generate charts in background ──
    import threading
    threading.Thread(
        target=generate_charts,
        args=(results[:15], query_label, all_seqs, promoter_tfs),
        daemon=True
    ).start()

    # ── Best alignment ──
    alignment_text = None
    non_self = [r for r in results if r['gene'] != query_label]
    if non_self and ('blosum' in metrics):
        db      = GeneDatabase(db_path=DB_PATH)
        top     = non_self[0]
        top_seq = db.get_isoform_sequence(top['gene'], actual_mode, top['isoform'])
        db.close()
        if top_seq:
            aln = (aligner if actual_mode == 'protein' else dna_aligner).align(query_seq, top_seq)
            alignment_text = str(aln[0])

    return jsonify({
        'query_label':    query_label,
        'isoform_id':     isoform_id,
        'query_length':   len(query_seq),
        'seq_mode':       actual_mode,
        'total_compared': len(all_seqs),
        'results':        results[:top_n],
        'alignment':      alignment_text,
        'best_non_self':  non_self[0] if non_self else None
    })


@app.route('/api/recent_searches')
def api_recent_searches():
    return jsonify(load_recent_searches())


@app.route('/api/clear_recent', methods=['POST'])
def api_clear_recent():
    try:
        if os.path.exists(RECENT_SEARCHES):
            os.remove(RECENT_SEARCHES)
    except Exception:
        pass
    return jsonify({'ok': True})


@app.route('/api/genes')
def api_genes():
    db     = GeneDatabase(db_path=DB_PATH)
    cursor = db.conn.cursor()
    cursor.execute('''
        SELECT g.symbol, g.organism, g.description,
               COUNT(CASE WHEN s.type = 'protein' THEN 1 END),
               COUNT(CASE WHEN s.type = 'dna'     THEN 1 END)
        FROM genes g
        LEFT JOIN sequences s ON s.gene_id = g.id
        GROUP BY g.id ORDER BY g.symbol
    ''')
    rows = cursor.fetchall()
    db.close()

    return jsonify([{
        'symbol':        r[0],
        'organism':      r[1],
        'description':   clean_gene_description(r[2], r[0]),
        'protein_count': r[3],
        'dna_count':     r[4],
    } for r in rows])


@app.route('/api/gene/<symbol>')
def api_gene(symbol):
    sym    = symbol.upper()
    db     = GeneDatabase(db_path=DB_PATH)
    info   = db.get_gene_info(sym)
    isoforms = db.get_isoforms(sym)
    db.close()

    if not info:
        return jsonify({'error': 'Gene not found'}), 404

    return jsonify({
        'symbol':      info[0],
        'organism':    info[1],
        'description': clean_gene_description(info[2], sym),
        'isoforms':    isoforms
    })


@app.route('/api/gene/<symbol>/sequence')
def api_gene_sequence(symbol):
    """Return the raw sequence for a specific isoform."""
    sym       = symbol.upper()
    iso_id    = request.args.get('isoform_id')
    seq_type  = request.args.get('type', 'protein')

    db  = GeneDatabase(db_path=DB_PATH)
    seq = db.get_isoform_sequence(sym, seq_type, iso_id)
    db.close()

    if not seq:
        return jsonify({'error': 'Sequence not found'}), 404

    return jsonify({'symbol': sym, 'isoform_id': iso_id, 'type': seq_type, 'sequence': seq, 'length': len(seq)})


@app.route('/api/benchmark', methods=['POST'])
def api_benchmark():
    n_pairs = int(request.json.get('n_pairs', 20))

    db       = GeneDatabase(db_path=DB_PATH)
    all_seqs = db.get_all_sequences_with_isoform("protein")
    db.close()

    seen = {}
    for symbol, _, seq in all_seqs:
        if symbol not in seen:
            seen[symbol] = seq

    genes     = list(seen.keys())
    all_pairs = [(a, b) for i, a in enumerate(genes) for b in genes[i+1:]]
    random.seed(42)
    pairs = random.sample(all_pairs, min(n_pairs, len(all_pairs)))

    tf_map = {s: get_promoter_tfs(s) for s in genes}

    metric_data = {
        'BLOSUM62':       {'times': [], 'scores': []},
        'K-mer (k=3)':    {'times': [], 'scores': []},
        'Edit distance':  {'times': [], 'scores': []},
        'Motif (JASPAR)': {'times': [], 'scores': []},
    }

    for sym_a, sym_b in pairs:
        sa, sb = seen[sym_a], seen[sym_b]

        t0 = time.perf_counter(); s = blosum_similarity(sa, sb)
        metric_data['BLOSUM62']['times'].append((time.perf_counter()-t0)*1000)
        metric_data['BLOSUM62']['scores'].append(s)

        t0 = time.perf_counter(); s = kmer_similarity(sa, sb, k=KMER_K)
        metric_data['K-mer (k=3)']['times'].append((time.perf_counter()-t0)*1000)
        metric_data['K-mer (k=3)']['scores'].append(s)

        t0 = time.perf_counter(); s = edit_distance_similarity(sa, sb)
        metric_data['Edit distance']['times'].append((time.perf_counter()-t0)*1000)
        metric_data['Edit distance']['scores'].append(s)

        t0 = time.perf_counter(); s = jaccard_motif_similarity(tf_map.get(sym_a,set()), tf_map.get(sym_b,set()))
        metric_data['Motif (JASPAR)']['times'].append((time.perf_counter()-t0)*1000)
        metric_data['Motif (JASPAR)']['scores'].append(s)

    blosum_scores = np.array(metric_data['BLOSUM62']['scores'])
    summary = []
    for name, d in metric_data.items():
        times  = np.array(d['times'])
        scores = np.array(d['scores'])
        corr   = float(np.corrcoef(blosum_scores, scores)[0,1]) if name != 'BLOSUM62' else 1.0
        summary.append({
            'metric': name,
            'avg_ms': round(float(np.mean(times)), 3),
            'min_ms': round(float(np.min(times)),  3),
            'max_ms': round(float(np.max(times)),  3),
            'corr':   round(corr, 4),
        })

    return jsonify({'pairs_tested': len(pairs), 'results': summary})


if __name__ == '__main__':
    os.makedirs(RESULTS_DIR, exist_ok=True)
    app.run(debug=True, port=5000)