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
def create_folder(folder_name):
    """Create a folder if it does not exist."""
    if not os.path.exists(folder_name):
        os.makedirs(folder_name)

def parse_filename(filename):
    """Parse the filename and return GroupID_UniqueID_abbrevation."""
    parts = filename.split('_')
    if len(parts) >= 4:
        return '_'.join(parts[:3])
    raise ValueError(f"Filename format is incorrect: {filename}")



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


from Bio import Entrez, SeqIO

def fetch_protein_sequence(protein_id, email="your_email@example.com"):
    """
    Fetch the protein sequence for a given protein ID from NCBI.
    
    Args:
        protein_id (str): The protein ID to fetch.
        email (str): Your email address for NCBI Entrez.

    Returns:
        str: The protein amino acid sequence in FASTA format.
    
    Raises:
        ValueError: If no sequence is found for the given protein ID.
    """
    Entrez.email = email  # NCBI requires an email for usage
    
    try:
        # Fetch the protein record from NCBI
        with Entrez.efetch(db="protein", id=protein_id, rettype="fasta", retmode="text") as handle:
            fasta_data = handle.read()
            if not fasta_data:
                raise ValueError(f"No sequence found for protein ID: {protein_id}")
            return fasta_data
    except Exception as e:
        raise ValueError(f"Error fetching sequence for protein ID {protein_id}: {str(e)}")



# Step 1: Create required folders
create_folder("out_sequences")
out_protein_folder = "out_proteins"
create_folder(out_protein_folder)

blast_result_folder = "blast_result"
create_folder(blast_result_folder)

# Step 2: List all files in the current folder
current_dir = "Proteins_ROs_alphaSubUnits"
files = [f for f in os.listdir(current_dir) if os.path.isfile(f)]

# Step 3: Process each file
for file in files:
    try:
        group_id = parse_filename(file)
        blast_file = os.path.join(blast_result_folder, group_id)
        blast_results = perform_blast(file, "nucleotide", 1e-10)
            with open("temp_blast_results.xml", "w") as temp:
                temp.write(blast_results)
        parse_blast_results("temp_blast_results.xml", blast_file, 0.0)
        os.remove("temp_blast_results.xml")
        print(f"Results saved to {blast_file}")

        # Parse filename
        
        incremental_id = 1

        
        # Read BLAST results
        with open(blast_file, 'r') as blast_results:
            next(blast_results)  # Skip header
            max_line = 10
            if max_line == 0:
                continue
            for line in blast_results:
                query_id, hit_id, accession, e_value, align_length = line.strip().split(',')
                if e_value != 0.0:
                    continue
                max_line = max_line - 1
                # Fetch genome sequence
                genome_sequence = fetch_genome_from_protein(accession)
                protein_sequence = fetch_protein_sequence(accession)

                # Save genome sequence to a file
                output_file = os.path.join("out_sequences", f"{group_id}_{incremental_id}.fasta")
                with open(output_file, 'w') as out_seq:
                    out_seq.write(genome_sequence)
                
                output_file = os.path.join(out_protein_folder, f"{group_id}_{incremental_id}.fasta")
                with open(output_file, 'w') as out_seq:
                    out_seq.write(protein_sequence)               


                incremental_id += 1
    except ValueError as e:
        print(f"Skipping file due to error: {e}")
