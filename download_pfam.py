#!/usr/bin/env python3

# standard library modules
import sys, json, ssl
import os
from urllib import request
from urllib.error import HTTPError
from time import sleep

BASE_URL = "https://www.ebi.ac.uk:443/interpro/api/protein/UniProt/entry/pfam/PF00355/?page_size=200"
UNIPROT_BASE_URL = "https://www.uniprot.org/uniprot/"

def get_sequence(accession):
    # Get the protein sequence from UniProt
    context = ssl._create_unverified_context()
    url = f"{UNIPROT_BASE_URL}{accession}.fasta"
    try:
        req = request.Request(url)
        res = request.urlopen(req, context=context)
        if res.status == 200:
            sequence_data = res.read().decode()
            return sequence_data
        else:
            print(f"Failed to get sequence for {accession}, status code: {res.status}")
            return None
    except HTTPError as e:
        print(f"Error fetching sequence for {accession}: {e}")
        return None

def output_list(output_dir="output_pfam"):
    # Create output directory if it doesn't exist
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Disable SSL verification to avoid configuration issues
    context = ssl._create_unverified_context()

    next = BASE_URL
    last_page = False

    attempts = 0
    while next:
        try:
            req = request.Request(next, headers={"Accept": "application/json"})
            res = request.urlopen(req, context=context)
            # If the API times out due to a long running query
            if res.status == 408:
                # wait just over a minute
                sleep(61)
                # continue this loop with the same URL
                continue
            elif res.status == 204:
                # No data, so leave loop
                break
            payload = json.loads(res.read().decode())
            next = payload["next"]
            attempts = 0
            if not next:
                last_page = True
        except HTTPError as e:
            if e.code == 408:
                sleep(61)
                continue
            else:
                # Retry up to 3 times before failing on other HTTP errors
                if attempts < 3:
                    attempts += 1
                    sleep(61)
                    continue
                else:
                    sys.stderr.write("LAST URL: " + next)
                    raise e

        for i, item in enumerate(payload["results"]):
            # Extract information for the header
            accession = item["metadata"].get("accession", "unknown")
            name = item["metadata"].get("name", "unknown")

            # Define file path for each protein
            file_path = os.path.join(output_dir, f"{accession}.fasta")

            # Check if the file already exists
            if os.path.exists(file_path):
                print(f"File {accession}.fasta already exists, skipping...")
                continue

            # Get the protein sequence
            sequence = get_sequence(accession)
            if sequence:
                # Write header and sequence to the file
                print('.',end='',flush=True)
                with open(file_path, "w") as f:
                    
                    f.write(sequence)
            else:
                print(f"Skipping {accession} due to sequence retrieval failure.")
            
            if last_page and i + 1 == len(payload["results"]):
                # End of output
                pass

        # Don't overload the server, give it time before asking for more
        if next:
            sleep(1)

if __name__ == "__main__":
    output_list()
