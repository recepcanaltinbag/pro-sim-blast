import requests
import pprint

def fetch_uniprot_metadata(uniprot_id):
    print(uniprot_id)
    url = f"https://www.uniprot.org/uniprotkb/{uniprot_id}.json"
    response = requests.get(url)
    
    if response.status_code == 200:
        data = response.json()
        
        # Extracting relevant metadata
        metadata = {}
        # RefSeq bilgilerini çekme
        # RefSeq bilgilerini çekme (Liste olarak)
        refseq_entries = [
            entry for entry in data.get('uniProtKBCrossReferences', []) if entry.get('database') == 'RefSeq'
        ]

        metadata['refseq_ids'] = [entry.get('id', 'N/A') for entry in refseq_entries]
        metadata['nucleotide_sequence_ids'] = [
            prop.get('value')
            for entry in refseq_entries
            for prop in entry.get('properties', [])
            if prop.get('key') == 'NucleotideSequenceId'
        ]
        metadata['protein_name'] = data.get('proteinDescription', {}).get('submissionNames', [{}])[0].get('fullName', {}).get('value', 'N/A')
        if metadata['protein_name'] == 'N/A':
            metadata['protein_name'] = data.get('proteinDescription', {}).get('recommendedName', {}).get('fullName', {}).get('value', 'N/A')
        metadata['refs'] = [ref.get("citation", {}).get("title", "N/A") for ref in data.get("references", [])]
        metadata['ec_number'] = data.get('proteinDescription', {}).get('submissionNames', [{}])[0].get('ecNumbers', [{}])[0].get('value', 'N/A')
        metadata['gene_names'] = data.get('genes', [{}])[0].get('geneName', {}).get('value', 'N/A')
        metadata['organism'] = data.get('organism', {}).get('scientificName', 'N/A')
        metadata['sequence'] = data.get('sequence', {}).get('value', 'N/A')
        metadata['go_terms'] = [go['term']['name'] for go in data.get('goTerms', [])]  # Gene Ontology terms
        metadata['function'] = data.get('comments', [{}])[0].get('text', 'N/A')  # Function description
        
        # Return the extracted metadata
        return metadata
    else:
        return None


def uniprot_id_to_info_dict(uniprot_id = 'A0A009H5P1'):
    metadata = fetch_uniprot_metadata(uniprot_id)

    if metadata:
        print(f"Metadata for UniProt ID {uniprot_id}:")
        for key, value in metadata.items():
            print(f"{key}: {value}")
        return metadata
    else:
        print(f"Failed to fetch metadata for {uniprot_id}")
        return None
