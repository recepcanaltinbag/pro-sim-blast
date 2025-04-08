import os
import argparse
import requests
from Bio.Blast import NCBIXML
import time
from Bio import Entrez

Entrez.email = "your_email@example.com"

def fetch_genome_from_protein(protein_id):
    # Step 1: Link Protein ID to Nucleotide ID
    handle = Entrez.elink(dbfrom="protein", db="nuccore", id=protein_id)
    record = Entrez.read(handle)
    handle.close()

    # Get the first nucleotide ID linked to the protein
    try:
        nucleotide_id = record[0]["LinkSetDb"][0]["Link"][0]["Id"]
        print(f"Nucleotide ID: {nucleotide_id}")
    except (IndexError, KeyError):
        print("No nucleotide ID found for this protein.")
        return None

    # Step 2: Fetch the Nucleotide Sequence
    handle = Entrez.efetch(db="nuccore", id=nucleotide_id, rettype="fasta", retmode="text")
    genome_sequence = handle.read()
    handle.close()

    return genome_sequence



def parse_args():
    parser = argparse.ArgumentParser(description="Perform NCBI BLAST search for nucleotide or protein FASTA files.")
    parser.add_argument(
        "-i", "--input", required=True, help="Path to the input FASTA file."
    )
    parser.add_argument(
        "-o", "--output", required=True, help="Path to save the results in CSV format."
    )
    parser.add_argument(
        "-e", "--evalue", type=float, default=1e-10, help="E-value threshold for BLAST hits. Default: 1e-10."
    )
    parser.add_argument(
        "-t", "--type", choices=["nucleotide", "protein"], required=True,
        help="Type of query sequence: 'nucleotide' or 'protein'."
    )
    return parser.parse_args()

def perform_blast(fasta_file, query_type, evalue_threshold):
    # Read the FASTA file content
    with open(fasta_file, "r") as f:
        fasta_content = f.read()

    # Set up the BLAST API request
    url = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"
    program = "blastn" if query_type == "nucleotide" else "blastp"
    database = "nt" if query_type == "nucleotide" else "nr"
    params = {
        "CMD": "Put",
        "DATABASE": database,
        "PROGRAM": program,
        "MEGABLAST": "on" if query_type == "nucleotide" else "off",
        "FILTER": "L",
        "EXPECT": evalue_threshold,
        "HITLIST_SIZE": 100,  # Retrieve top 50 hits
        "FORMAT_TYPE": "XML"
    }
    files = {"QUERY": fasta_content}

    # Submit the BLAST search
    print("Submitting BLAST search...")
    response = requests.post(url, data=params, files=files)
    response.raise_for_status()
    rid = response.text.split("RID = ")[1].split("\n")[0]
    print(f"BLAST job submitted. RID: {rid}")

    # Check job status and retrieve results
    print("Checking BLAST job status...")
    params = {"CMD": "Get", "RID": rid, "FORMAT_TYPE": "XML"}
    while True:
        result = requests.get(url, params=params)
        if "Status=WAITING" not in result.text:
            break
        print("Job still running, waiting...")
        time.sleep(10)

    return result.text



def parse_blast_results(blast_xml, output_file, evalue_threshold):
    with open(output_file, "w") as out:
        out.write("Query_ID,Subject_ID,Accession,E-value,Alignment_Length\n")

        # Parse the BLAST XML results
        blast_records = NCBIXML.parse(open(blast_xml))
        for record in blast_records:
            for alignment in record.alignments:
                for hsp in alignment.hsps:
                    if hsp.expect <= evalue_threshold:
                        out.write(
                            f"{record.query_id},{alignment.hit_id},{alignment.accession},{hsp.expect},{hsp.align_length}\n"
                        )

if __name__ == "__main__":
    args = parse_args()
    blast_results = perform_blast(args.input, args.type, 1e-10)
    with open("temp_blast_results.xml", "w") as temp:
        temp.write(blast_results)
    parse_blast_results("temp_blast_results.xml", args.output, 0.0)
    os.remove("temp_blast_results.xml")
    print(f"Results saved to {args.output}")

