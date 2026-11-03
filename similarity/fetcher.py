"""
fetcher.py
----------
Unified client for fetching gene and protein/DNA sequences from
NCBI RefSeq, UniProt, and Ensembl REST APIs.
"""

import requests
from io import StringIO
from Bio import SeqIO, Entrez
from config import ENTREZ_EMAIL

Entrez.email = ENTREZ_EMAIL


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
            res = requests.get(url_lookup, headers=headers, timeout=8)
        except Exception:
            url_lookup = f"https://rest.ensembl.org/lookup/symbol/homo_sapiens/{gene_name}?expand=1"
            res = requests.get(url_lookup, headers=headers, timeout=8)

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

        url_seq = f"http://rest.ensembl.org/sequence/id/{target_id}"
        try:
            res_seq = requests.get(url_seq, headers=headers, timeout=8)
        except Exception:
            url_seq = f"https://rest.ensembl.org/sequence/id/{target_id}"
            res_seq = requests.get(url_seq, headers=headers, timeout=8)

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

