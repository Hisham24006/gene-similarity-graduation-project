import requests
from Bio import Entrez, SeqIO
from io import StringIO

class BioFetcher:
    def __init__(self, email):
        Entrez.email = email

    def fetch_ncbi(self, symbol, db_type="protein"):
        """Fetches DNA or Protein from NCBI."""
        term = f"{symbol}[Gene Name] AND human[Organism]"
        handle = Entrez.esearch(db=db_type, term=term, retmax=1)
        record = Entrez.read(handle)
        if record["IdList"]:
            id = record["IdList"][0]
            fetch = Entrez.efetch(db=db_type, id=id, rettype="fasta", retmode="text")
            return str(SeqIO.read(StringIO(fetch.read()), "fasta").seq)
        return None

    def fetch_uniprot(self, symbol):
        """Fetches Protein sequence from UniProt REST API."""
        url = f"https://rest.uniprot.org/uniprotkb/search?query=gene:{symbol}+AND+organism_id:9606&format=fasta"
        response = requests.get(url)
        if response.status_code == 200 and response.text:
            return str(SeqIO.read(StringIO(response.text), "fasta").seq)
        return None
