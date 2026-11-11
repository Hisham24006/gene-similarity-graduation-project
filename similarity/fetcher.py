"""
fetcher.py
----------
Unified client for fetching gene and protein/DNA sequences from
NCBI RefSeq, UniProt, and Ensembl REST APIs.

Two ways to fetch:
    fetch_from_ncbi / fetch_from_uniprot / fetch_from_ensembl(gene_symbol, seq_type)
        Look a gene up by symbol on the chosen database.
    fetch_by_accession(accession, seq_type)
        Fetch one exact record by its accession number; the database and
        sequence type follow from the accession format (see detect_accession).
"""

import re
import requests
from io import StringIO
from Bio import SeqIO, Entrez
from config import ENTREZ_EMAIL

Entrez.email = ENTREZ_EMAIL

ENSEMBL_REST = "https://rest.ensembl.org"
ENSEMBL_HEADERS = {'Content-Type': 'application/json', 'User-Agent': 'GeneSim/1.0'}
ENSEMBL_TIMEOUT = 30  # Ensembl REST can take 10-40 s at busy times


# -----------------------------------------------------------------------
# Accession numbers
# -----------------------------------------------------------------------

# UniProt accession format (from uniprot.org/help/accession_numbers),
# optionally with an isoform suffix such as P04637-2.
_UNIPROT = r'(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})(?:-\d+)?'

# (regex, database, sequence type, label). A sequence type of None means
# the record is a gene, so the protein/DNA choice comes from the user.
ACCESSION_PATTERNS = [
    (re.compile(r'(?:NP|XP)_\d+(?:\.\d+)?'),       'ncbi',    'protein', 'RefSeq protein'),
    (re.compile(r'(?:NM|XM|NR|XR)_\d+(?:\.\d+)?'), 'ncbi',    'dna',     'RefSeq RNA'),
    (re.compile(r'ENSP\d{11}(?:\.\d+)?'),          'ensembl', 'protein', 'Ensembl protein'),
    (re.compile(r'ENST\d{11}(?:\.\d+)?'),          'ensembl', 'dna',     'Ensembl transcript'),
    (re.compile(r'ENSG\d{11}(?:\.\d+)?'),          'ensembl', None,      'Ensembl gene'),
    (re.compile(_UNIPROT),                         'uniprot', 'protein', 'UniProt'),
]


def detect_accession(text):
    """
    Returns {'accession', 'source', 'seq_type', 'kind'} if text is an
    accession number in a supported format, otherwise None (e.g. a symbol).
    """
    acc = (text or '').strip().upper()
    for pattern, source, seq_type, kind in ACCESSION_PATTERNS:
        if pattern.fullmatch(acc):
            return {'accession': acc, 'source': source, 'seq_type': seq_type, 'kind': kind}
    return None


def _result(isoform_id=None, sequence=None, description=None, symbol=None,
            seq_type=None, source=None, kind=None, organism=None, error=None):
    return {
        'isoform_id': isoform_id, 'sequence': sequence, 'description': description,
        'symbol': symbol, 'seq_type': seq_type, 'source': source, 'kind': kind,
        'organism': organism, 'error': error,
    }


def _ensembl_get(path):
    """GET an Ensembl REST path; returns parsed JSON or None if not found."""
    res = requests.get(f"{ENSEMBL_REST}{path}", headers=ENSEMBL_HEADERS, timeout=ENSEMBL_TIMEOUT)
    return res.json() if res.status_code == 200 else None


def _symbol_from_transcript_name(name):
    """Ensembl transcript names are '<gene symbol>-<number>', e.g. TP53-201."""
    return name.rsplit('-', 1)[0] if name and '-' in name else name


def _fetch_ncbi_accession(acc, seq_type, kind):
    db, rettype = ('protein', 'gp') if seq_type == 'protein' else ('nuccore', 'gb')
    try:
        handle = Entrez.efetch(db=db, id=acc, rettype=rettype, retmode='text')
        record = SeqIO.read(StringIO(handle.read()), 'genbank')
        handle.close()
        sequence = str(record.seq)
    except Exception:
        return _result(error=f"{kind} accession '{acc}' was not found on NCBI.")

    symbol = None
    for feature in record.features:
        if feature.type in ('gene', 'CDS', 'Protein') and 'gene' in feature.qualifiers:
            symbol = feature.qualifiers['gene'][0]
            break
    return _result(record.id, sequence, record.description, symbol, seq_type, 'ncbi', kind,
                   record.annotations.get('organism'))


def _fetch_uniprot_accession(acc, kind):
    try:
        res = requests.get(f"https://rest.uniprot.org/uniprotkb/{acc}.fasta", timeout=15)
    except requests.exceptions.RequestException as e:
        return _result(error=f"UniProt error: {e}")
    if res.status_code != 200 or not res.text.startswith('>'):
        return _result(error=f"UniProt accession '{acc}' was not found.")

    record = SeqIO.read(StringIO(res.text), 'fasta')
    gene = re.search(r'\bGN=(\S+)', record.description)
    organism = re.search(r'\bOS=(.+?)\s+OX=', record.description)
    return _result(record.id, str(record.seq), record.description,
                   gene.group(1) if gene else None, 'protein', 'uniprot', kind,
                   organism.group(1) if organism else None)


def _fetch_ensembl_accession(acc, seq_type, kind):
    stable_id = acc.split('.')[0]  # Ensembl REST does not accept version suffixes
    try:
        info = _ensembl_get(f"/lookup/id/{stable_id}")
        if not info:
            return _result(error=f"{kind} ID '{stable_id}' was not found on Ensembl.")
        species = (info.get('species') or '').replace('_', ' ').capitalize() or None

        if info.get('object_type') == 'Gene':
            # A gene ID: fetch its canonical transcript like a symbol lookup
            symbol = info.get('display_name') or stable_id
            iso_id, seq, desc, error = fetch_from_ensembl(symbol, seq_type or 'protein')
            if error:
                return _result(error=error)
            return _result(iso_id, seq, desc, symbol, seq_type or 'protein', 'ensembl', kind, species)

        if info.get('object_type') == 'Translation':
            transcript = _ensembl_get(f"/lookup/id/{info.get('Parent')}") or {}
            symbol = _symbol_from_transcript_name(transcript.get('display_name'))
            seq_kind = 'protein'
        else:
            symbol = _symbol_from_transcript_name(info.get('display_name'))
            seq_kind = 'cdna'  # spliced mRNA; the default 'genomic' includes introns

        seq_data = _ensembl_get(f"/sequence/id/{stable_id}?type={seq_kind}")
        if not seq_data or not seq_data.get('seq'):
            return _result(error=f"Could not retrieve the sequence for '{stable_id}' from Ensembl.")
        desc = f"{symbol or stable_id} | {stable_id} | {species or 'Ensembl'}"
        return _result(stable_id, seq_data['seq'], desc, symbol, seq_type, 'ensembl', kind, species)

    except requests.exceptions.Timeout:
        return _result(error="Ensembl REST API timed out. Please try again, or use an NCBI or UniProt accession.")
    except requests.exceptions.RequestException as e:
        return _result(error=f"Ensembl error: {e}")


def fetch_by_accession(accession, seq_type='protein'):
    """
    Fetches one record by accession number. Supported formats:
        RefSeq   NP_/XP_ (protein), NM_/XM_/NR_/XR_ (RNA)        -> NCBI
        UniProt  e.g. P04637, P04637-2 (isoform)                  -> UniProt
        Ensembl  ENSP (protein), ENST (transcript, as cDNA),
                 ENSG (gene: canonical transcript of seq_type)    -> Ensembl

    seq_type only matters for gene IDs (ENSG); other accessions determine
    their own type. Returns a dict with isoform_id, sequence, description,
    symbol (gene symbol, if the source reports one), seq_type, source,
    kind, organism and error.
    """
    info = detect_accession(accession)
    if not info:
        return _result(error=f"'{accession}' is not a recognised accession number.")
    acc, source, kind = info['accession'], info['source'], info['kind']

    if source == 'ncbi':
        return _fetch_ncbi_accession(acc, info['seq_type'], kind)
    if source == 'uniprot':
        return _fetch_uniprot_accession(acc, kind)
    return _fetch_ensembl_accession(acc, info['seq_type'] or seq_type, kind)


def fetch_from_ncbi(gene_name, seq_type='protein'):
    """
    Fetch protein or DNA sequence from NCBI RefSeq.
    Returns (isoform_id, sequence, description, error_message).
    """
    try:
        if seq_type == 'protein':
            db_name = "protein"
            term = f"{gene_name}[Gene Name] AND human[Organism] AND srcdb_refseq[PROP]"
        else:
            db_name = "nucleotide"
            term = f"{gene_name}[Gene Name] AND Homo sapiens[Organism] AND mRNA[Filter] AND srcdb_refseq[PROP]"

        handle = Entrez.esearch(db=db_name, term=term, retmax=1)
        record = Entrez.read(handle)
        handle.close()

        if not record.get("IdList"):
            return None, None, None, f"No {seq_type} results found on NCBI for '{gene_name}'"

        handle = Entrez.efetch(db=db_name, id=record["IdList"][0], rettype="fasta", retmode="text")
        seq_record = SeqIO.read(StringIO(handle.read()), "fasta")
        handle.close()
        return seq_record.id, str(seq_record.seq), seq_record.description, None

    except Exception as e:
        return None, None, None, str(e)


def fetch_from_uniprot(gene_name, seq_type='protein'):
    """
    Fetch protein sequence from UniProt.
    Returns (isoform_id, sequence, description, error_message).
    """
    if seq_type == 'dna':
        return (
            None,
            None,
            None,
            "UniProt is a curated protein sequence database and does not store nucleotide/DNA sequences. Please choose NCBI or Ensembl to fetch DNA."
        )

    try:
        # 1. Try reviewed (Swiss-Prot)
        url = (
            f"https://rest.uniprot.org/uniprotkb/search"
            f"?query=gene_exact:{gene_name}+AND+organism_id:9606+AND+reviewed:true"
            f"&format=fasta&size=1"
        )
        res = requests.get(url, timeout=12)
        if res.status_code == 200 and res.text.strip():
            seq_record = SeqIO.read(StringIO(res.text), "fasta")
            return seq_record.id, str(seq_record.seq), seq_record.description, None

        # 2. Fallback to unreviewed if not found
        url2 = (
            f"https://rest.uniprot.org/uniprotkb/search"
            f"?query=gene_exact:{gene_name}+AND+organism_id:9606"
            f"&format=fasta&size=1"
        )
        res2 = requests.get(url2, timeout=12)
        if res2.status_code == 200 and res2.text.strip():
            seq_record = SeqIO.read(StringIO(res2.text), "fasta")
            return seq_record.id, str(seq_record.seq), seq_record.description, None

        return None, None, None, f"No protein sequence found on UniProt for '{gene_name}'"

    except Exception as e:
        return None, None, None, f"UniProt error: {str(e)}"


def fetch_from_ensembl(gene_name, seq_type='protein'):
    """
    Fetch sequence from Ensembl REST API.
    Returns (isoform_id, sequence, description, error_message).
    """
    headers = {'Content-Type': 'application/json', 'User-Agent': 'GeneSim/1.0'}
    try:
        url_lookup = f"http://rest.ensembl.org/lookup/symbol/homo_sapiens/{gene_name}?expand=1"
        try:
            res = requests.get(url_lookup, headers=headers, timeout=ENSEMBL_TIMEOUT)
        except Exception:
            url_lookup = f"https://rest.ensembl.org/lookup/symbol/homo_sapiens/{gene_name}?expand=1"
            res = requests.get(url_lookup, headers=headers, timeout=ENSEMBL_TIMEOUT)

        if res.status_code != 200:
            return None, None, None, f"Gene '{gene_name}' not found on Ensembl"

        data = res.json()
        gene_id = data.get("id")
        display_name = data.get("display_name", gene_name)
        desc = data.get("description", "")

        transcripts = data.get("Transcript", [])
        target_id = None
        for tr in transcripts:
            if tr.get("is_canonical") == 1:
                if seq_type == "protein" and tr.get("Translation"):
                    target_id = tr["Translation"].get("id")
                else:
                    target_id = tr.get("id")
                break

        if not target_id and transcripts:
            tr0 = transcripts[0]
            target_id = tr0["Translation"].get("id") if (seq_type == "protein" and tr0.get("Translation")) else tr0.get("id")

        if not target_id:
            target_id = gene_id

        # Transcripts default to the genomic sequence (introns included);
        # request the spliced cDNA so DNA fetches match mRNA records.
        seq_param = "?type=cdna" if target_id.startswith("ENST") else ""
        url_seq = f"http://rest.ensembl.org/sequence/id/{target_id}{seq_param}"
        try:
            res_seq = requests.get(url_seq, headers=headers, timeout=ENSEMBL_TIMEOUT)
        except Exception:
            url_seq = f"https://rest.ensembl.org/sequence/id/{target_id}{seq_param}"
            res_seq = requests.get(url_seq, headers=headers, timeout=ENSEMBL_TIMEOUT)

        if res_seq.status_code == 200:
            seq = res_seq.json().get("seq")
            if seq:
                desc_str = f"{display_name} | {target_id} | {desc}" if desc else f"{display_name} | {target_id} | Homo sapiens"
                return target_id, str(seq), desc_str, None

        return None, None, None, f"Could not retrieve sequence from Ensembl for '{gene_name}'"

    except requests.exceptions.Timeout:
        return None, None, None, "Ensembl REST API timed out. Ensembl servers may be experiencing high latency; please try NCBI or UniProt."
    except Exception as e:
        return None, None, None, f"Ensembl error: {str(e)}"

