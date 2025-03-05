# Python script to convert Ensembl IDs to FASTA sequences using NCBI's genomic coordinates
# Requires: pip install biopython

"""
This script handles:
1. Ensembl ID version stripping
2. Coordinate conversion from 1-based (Ensembl) to 0-based (NCBI)
3. Chromosome mapping for GRCh38 assembly
4. FASTA sequence retrieval from NCBI

Key features:
- Uses Ensembl REST API for gene location data[6][7]
- Converts chromosome numbers to NCBI RefSeq format[3][8]
- Handles X/Y chromosome mapping[8]
- Uses Biopython's Entrez module for NCBI access[3]

For the example ID `ENSG00000092847.14`, this code retrieves the FASTA sequence between positions 35,869,761 and 35,930,532 on chromosome 1 (NC_000001.11)[1][8].

Note: Replace "your.email@example.com" with your actual email address for NCBI compliance.

Citations:
[1] https://www.ncbi.nlm.nih.gov/nuccore/NC_000001.11?report=fasta&from=35869761&to=35930532
[2] https://support.bioconductor.org/p/9145757/
[3] https://guides.lib.berkeley.edu/ncbi/sequences
[4] https://www.biostars.org/p/428605/
[5] https://support.bioconductor.org/p/9155010/
[6] https://www.youtube.com/watch?v=OIAuBoUPOHw
[7] https://www.ensembl.org/Homo_sapiens/Gene/Summary?g=ENSG00000141510
[8] https://www.ncbi.nlm.nih.gov/gene?Db=gene&Cmd=DetailsSearch&Term=26523
[9] https://www.uniprot.org/uniprot/A0A8V8TPF2-1
[10] https://www.uniprot.org/uniprotkb/Q9UL18/entry
[11] https://david.ncifcrf.gov/conversion.jsp
[12] https://www.ensembl.org/Homo_sapiens/Gene/Summary?g=ENSG00000114771
[13] https://www.biorxiv.org/content/10.1101/2024.06.12.598692v1.full.pdf
[14] https://www.ncbi.nlm.nih.gov/gene/

"""

import requests
from Bio import Entrez

def get_genomic_coordinates(ensembl_id):
    """Retrieve chromosome, start, and end positions from Ensembl REST API"""
    base_id = ensembl_id.split('.')[0]  # Remove version number
    url = f"https://rest.ensembl.org/lookup/id/{base_id}?expand=0"
    
    response = requests.get(url, headers={"Content-Type": "application/json"})
    response.raise_for_status()
    
    data = response.json()
    return {
        'chrom': data['seq_region_name'],
        'start': data['start'] - 1,  # Adjust for 0-based indexing
        'end': data['end']
    }

def get_ncbi_refseq_id(chrom):
    """Construct NCBI RefSeq ID from chromosome number"""
    chrom_map = {'X': '23', 'Y': '24'}
    chrom_num = chrom_map.get(chrom, chrom.zfill(2))
    return f"NC_{chrom_num.rjust(6, '0')}.11"  # GRCh38 assembly format

def fetch_fasta(ensembl_id, email = "lars.ericson@catskillsresearch.com"):
    coordinates = get_genomic_coordinates(ensembl_id)
    refseq_id = get_ncbi_refseq_id(coordinates['chrom'])
    
    Entrez.email = email # Required by NCBI
    handle = Entrez.efetch(
        db="nuccore",
        id=refseq_id,
        rettype="fasta",
        retmode="text",
        seq_start=coordinates['start'],
        seq_stop=coordinates['end']
    )
    return handle.read()

if __name__=="__main__":
    fasta_sequence = fetch_fasta("ENSG00000092847.14")
    print(fasta_sequence)
