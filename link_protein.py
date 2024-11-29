from Bio import Entrez

# Set your email for NCBI Entrez
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

# Example: Fetch genome for a Protein ID
protein_id = "WP_068587002"  # Replace with your Protein ID
genome_sequence = fetch_genome_from_protein(protein_id)
if genome_sequence:
    print(f"Genome Sequence:\n{genome_sequence}")

