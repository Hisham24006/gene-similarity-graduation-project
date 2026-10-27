"""
app.py — Gene Similarity Search Flask App
"""

import os
import sys
import glob
import json
import base64
import time
import random
import threading

import numpy as np
from flask import Flask, render_template, request, jsonify

# ── Dynamic Path Setup ──
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SIM_DIR = os.path.join(BASE_DIR, 'similarity')
DB_DIR = os.path.join(BASE_DIR, 'database')
sys.path.insert(0, SIM_DIR)
sys.path.insert(0, DB_DIR)

from config import (
    DB_PATH,
    RESULTS_DIR,
    RECENT_SEARCHES_PATH as RECENT_SEARCHES,
    KMER_K,
    VALID_AA,
    VALID_DNA
)
from database_manager import GeneDatabase
from metrics import (
    kmer_similarity,
    edit_distance_similarity,
    blosum_similarity,
    dna_similarity,
    protein_aligner,
    dna_aligner
)
from motif_similarity import (
    load_jaspar_motifs,
    fetch_promoter,
    cached_scan_motifs,
    jaccard_motif_similarity
)
from fetcher import (
    fetch_from_ncbi,
    fetch_from_uniprot,
    fetch_from_ensembl
)
from visualizer import plot_bar_chart, plot_heatmap

app = Flask(__name__)

# ── Startup ──
print("Loading JASPAR motifs...")
BIO_MOTIFS = load_jaspar_motifs()
print(f"Loaded {len(BIO_MOTIFS)} motif models.")


# ── Helpers ──

def get_promoter_tfs(symbol):
    promoter = fetch_promoter(symbol)
    return cached_scan_motifs(symbol, promoter, BIO_MOTIFS) if promoter else set()


def image_to_base64(path):
    if not os.path.exists(path):
        return None
    with open(path, 'rb') as f:
        return base64.b64encode(f.read()).decode('utf-8')


def save_recent_search(query_label, seq_mode, metrics, result_count):
    """Save a search to the recent searches log (max 20 entries)."""
    try:
        searches = []
        if os.path.exists(RECENT_SEARCHES):
            with open(RECENT_SEARCHES, 'r') as f:
                searches = json.load(f)

        entry = {
            'query':   query_label,
            'mode':    seq_mode,
            'metrics': metrics,
            'results': result_count,
            'time':    time.strftime('%Y-%m-%d %H:%M'),
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


@app.route('/fetch')
def fetch_page():
    return render_template('fetch.html')


# ── API routes ──

@app.route('/api/fetch', methods=['POST'])
def api_fetch():
    data = request.json or {}
    gene_name = data.get('gene_name', '').strip().upper()
    source = data.get('source', 'ncbi').strip().lower()
    seq_type = data.get('seq_type', 'protein').strip().lower()

    if not gene_name:
        return jsonify({'error': 'Please enter a gene symbol (e.g. TP53, BRCA1, EGFR).'}), 400

    if source == 'ncbi':
        iso_id, seq, desc, error = fetch_from_ncbi(gene_name, seq_type)
        source_name = 'NCBI RefSeq'
    elif source == 'uniprot':
        iso_id, seq, desc, error = fetch_from_uniprot(gene_name, seq_type)
        source_name = 'UniProt'
    elif source == 'ensembl':
        iso_id, seq, desc, error = fetch_from_ensembl(gene_name, seq_type)
        source_name = 'Ensembl'
    else:
        return jsonify({'error': f"Unknown website '{source}'. Please choose NCBI, UniProt, or Ensembl."}), 400

    if error:
        return jsonify({'error': error}), 400

    fasta_header = f">{iso_id} {desc or (gene_name + ' ' + seq_type)}"
    fasta_str = f"{fasta_header}\n{seq}"

    source_url = None
    if source == 'ncbi':
        if seq_type == 'protein':
            source_url = f"https://www.ncbi.nlm.nih.gov/protein/{iso_id}"
        else:
            source_url = f"https://www.ncbi.nlm.nih.gov/nuccore/{iso_id}"
    elif source == 'uniprot':
        uniprot_acc = iso_id.split('|')[1] if '|' in iso_id else iso_id
        source_url = f"https://www.uniprot.org/uniprotkb/{uniprot_acc}/entry"
    elif source == 'ensembl':
        source_url = f"https://www.ensembl.org/id/{iso_id}"

    return jsonify({
        'symbol': gene_name,
        'source': source_name,
        'source_key': source,
        'seq_type': seq_type,
        'isoform_id': iso_id,
        'sequence': seq,
        'length': len(seq),
        'description': desc or f"{gene_name} {seq_type} sequence from {source_name}",
        'fasta': fasta_str,
        'source_url': source_url
    })


@app.route('/api/search', methods=['POST'])
def api_search():
    data = request.json
    input_type = data.get('input_type', 'gene_name')
    fetch_type = data.get('fetch_type', 'protein')   # 'protein' or 'dna'
    seq_mode = data.get('seq_mode', 'protein')      # what to search against
    metrics = data.get('metrics', ['blosum', 'kmer', 'edit', 'motif'])
    top_n = int(data.get('top_n', 10))

    # ── Get query sequence ──
    if input_type == 'gene_name':
        gene_name = data.get('gene_name', '').upper().strip()
        isoform_id, query_seq, _, error = fetch_from_ncbi(gene_name, fetch_type)
        if error:
            return jsonify({'error': error}), 400
        query_label = gene_name
        actual_mode = fetch_type

    else:
        raw = data.get('sequence', '').strip()
        lines = [l.strip() for l in raw.splitlines() if l.strip()]
        if lines and lines[0].startswith('>'):
            isoform_id = lines[0][1:].split()[0]
            query_seq = ''.join(lines[1:]).upper().replace(' ', '')
            query_label = isoform_id
        else:
            query_seq = ''.join(lines).upper().replace(' ', '')
            isoform_id = 'custom'
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
    db = GeneDatabase(db_path=DB_PATH)
    all_seqs = db.get_all_sequences_with_isoform(actual_mode)
    db.close()

    if not all_seqs:
        return jsonify({'error': f'No {actual_mode} sequences found in database.'}), 400

    # ── Motif TF sets ──
    promoter_tfs = {}
    query_tfs = set()
    if 'motif' in metrics:
        query_tfs = get_promoter_tfs(query_label) if query_label != 'custom_sequence' else set()
        gene_symbols = list(set(s for s, i, _ in all_seqs))
        for symbol in gene_symbols:
            promoter_tfs[symbol] = get_promoter_tfs(symbol)

    # ── Run metrics ──
    results = []
    for symbol, iso_id, local_seq in all_seqs:
        row = {'gene': symbol, 'isoform': iso_id, 'scores': {}, 'average': 0.0}
        total = 0
        count = 0

        if actual_mode == 'protein':
            if 'blosum' in metrics:
                s = blosum_similarity(query_seq, local_seq)
                row['scores']['blosum'] = round(s, 4)
                total += s
                count += 1
        else:
            if 'blosum' in metrics:
                s = dna_similarity(query_seq, local_seq)
                row['scores']['blosum'] = round(s, 4)
                total += s
                count += 1

        if 'kmer' in metrics:
            s = kmer_similarity(query_seq, local_seq, k=KMER_K)
            row['scores']['kmer'] = round(s, 4)
            total += s
            count += 1

        if 'edit' in metrics:
            s = edit_distance_similarity(query_seq, local_seq)
            row['scores']['edit'] = round(s, 4)
            total += s
            count += 1

        if 'motif' in metrics:
            s = jaccard_motif_similarity(query_tfs, promoter_tfs.get(symbol, set()))
            row['scores']['motif'] = round(s, 4)
            total += s
            count += 1

        row['average'] = round(total / count, 4) if count > 0 else 0.0
        results.append(row)

    results.sort(key=lambda x: x['average'], reverse=True)

    # ── Save recent search ──
    save_recent_search(query_label, actual_mode, metrics, len(results))

    # ── Generate charts in background ──
    threading.Thread(
        target=generate_charts,
        args=(results[:15], query_label, all_seqs, promoter_tfs),
        daemon=True
    ).start()

    # ── Best alignment ──
    alignment_text = None
    non_self = [r for r in results if r['gene'] != query_label]
    if non_self and ('blosum' in metrics):
        db = GeneDatabase(db_path=DB_PATH)
        top = non_self[0]
        top_seq = db.get_isoform_sequence(top['gene'], actual_mode, top['isoform'])
        db.close()
        if top_seq:
            aln = (protein_aligner if actual_mode == 'protein' else dna_aligner).align(query_seq, top_seq)
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
    db = GeneDatabase(db_path=DB_PATH)
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
        'description':   r[2],
        'protein_count': r[3],
        'dna_count':     r[4],
    } for r in rows])


@app.route('/api/gene/<symbol>')
def api_gene(symbol):
    sym = symbol.upper()
    db = GeneDatabase(db_path=DB_PATH)
    info = db.get_gene_info(sym)
    isoforms = db.get_isoforms(sym)
    db.close()

    if not info:
        return jsonify({'error': 'Gene not found'}), 404

    return jsonify({
        'symbol':      info[0],
        'organism':    info[1],
        'description': info[2],
        'isoforms':    isoforms
    })


@app.route('/api/gene/<symbol>/sequence')
def api_gene_sequence(symbol):
    """Return the raw sequence for a specific isoform."""
    sym = symbol.upper()
    iso_id = request.args.get('isoform_id')
    seq_type = request.args.get('type', 'protein')

    db = GeneDatabase(db_path=DB_PATH)
    seq = db.get_isoform_sequence(sym, seq_type, iso_id)
    db.close()

    if not seq:
        return jsonify({'error': 'Sequence not found'}), 404

    return jsonify({'symbol': sym, 'isoform_id': iso_id, 'type': seq_type, 'sequence': seq, 'length': len(seq)})


@app.route('/api/benchmark', methods=['POST'])
def api_benchmark():
    n_pairs = int(request.json.get('n_pairs', 20))

    db = GeneDatabase(db_path=DB_PATH)
    all_seqs = db.get_all_sequences_with_isoform("protein")
    db.close()

    seen = {}
    for symbol, _, seq in all_seqs:
        if symbol not in seen:
            seen[symbol] = seq

    genes = list(seen.keys())
    all_pairs = [(a, b) for i, a in enumerate(genes) for b in genes[i + 1:]]
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

        t0 = time.perf_counter()
        s = blosum_similarity(sa, sb)
        metric_data['BLOSUM62']['times'].append((time.perf_counter() - t0) * 1000)
        metric_data['BLOSUM62']['scores'].append(s)

        t0 = time.perf_counter()
        s = kmer_similarity(sa, sb, k=KMER_K)
        metric_data['K-mer (k=3)']['times'].append((time.perf_counter() - t0) * 1000)
        metric_data['K-mer (k=3)']['scores'].append(s)

        t0 = time.perf_counter()
        s = edit_distance_similarity(sa, sb)
        metric_data['Edit distance']['times'].append((time.perf_counter() - t0) * 1000)
        metric_data['Edit distance']['scores'].append(s)

        t0 = time.perf_counter()
        s = jaccard_motif_similarity(tf_map.get(sym_a, set()), tf_map.get(sym_b, set()))
        metric_data['Motif (JASPAR)']['times'].append((time.perf_counter() - t0) * 1000)
        metric_data['Motif (JASPAR)']['scores'].append(s)

    blosum_scores = np.array(metric_data['BLOSUM62']['scores'])
    summary = []
    for name, d in metric_data.items():
        times = np.array(d['times'])
        scores = np.array(d['scores'])
        corr = float(np.corrcoef(blosum_scores, scores)[0, 1]) if name != 'BLOSUM62' else 1.0
        summary.append({
            'metric': name,
            'avg_ms': round(float(np.mean(times)), 3),
            'min_ms': round(float(np.min(times)), 3),
            'max_ms': round(float(np.max(times)), 3),
            'corr':   round(corr, 4),
        })

    return jsonify({'pairs_tested': len(pairs), 'results': summary})


if __name__ == '__main__':
    os.makedirs(RESULTS_DIR, exist_ok=True)
    app.run(debug=True, port=5000)
