from pymongo.mongo_client import MongoClient
from dotenv import load_dotenv
import os
from pprint import pprint
from f_uniprot_metadata_retrieve import uniprot_id_to_info_dict
from filter_hmm_out import parse_filter_hmm
import subprocess
from Bio import SeqIO

def get_seq_ids(fasta_path):
    seq_ids = []
    # FASTA dosyasını okuyarak her bir sekansın ID'sini al
    for record in SeqIO.parse(fasta_path, "fasta"):
        seq_ids.append(record.id)  # seq_id'yi alır
    return seq_ids

load_dotenv()
MONGODB_URI = os.environ['MONGODB_URI']


def fasta_one_by_one(input_dir, hmm_file, e_threshold, Rieske_Collection):
        # Klasördeki tüm FASTA dosyalarını bul
    for file1 in os.listdir(input_dir):
        if file1.endswith(".fasta"):
            fasta_path = os.path.join(input_dir, file1)
            tblout_path = "pfamPart.out"
            print(file1)
            # hmmsearch komutunu çalıştır
            command = [
                "hmmsearch",
                "--tblout", tblout_path,
                hmm_file,
                fasta_path
            ]
            subprocess.run(command, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            print('End of Hmmer')
            # Çıktıyı işle, e_threshold = 1e-60
            processed_result = parse_filter_hmm(tblout_path, e_threshold)
            if processed_result == None:
                print('No result for this e threshold')
                the_id = get_seq_ids(fasta_path)[0]
                no_res = {'_id': the_id, 'PredictedCluster': 'N/A'}
                
                try:
                    Rieske_Collection.insert_one(no_res)
                    print(f"Inserted: {the_id}")
                except Exception as e:
                    print(f"Skipped duplicate entry with _id: {the_id}",e)
            else:
                extra_info_dict = uniprot_id_to_info_dict(processed_result["Accession Name"].split('|')[1])
                
                print(processed_result)
                print('EXTRA')
                print(extra_info_dict)
                processed_result['_id'] = processed_result['UniProt Accession']
                combined_data = {**processed_result, **extra_info_dict}  # İki sözlüğü birleştir
                the_id = processed_result['UniProt Accession']
                try:
                    Rieske_Collection.insert_one(combined_data)
                    print(f"Inserted: {the_id}")
                except Exception as e:
                    print(f"Skipped duplicate entry with _id: {the_id}",e)


e_threshold = 1e-10
hmm_file = "RieskeDB.hmm"
protein_folder = "output_pfam"


# Create a new client and connect to the server
client = MongoClient(MONGODB_URI)
# Send a ping to confirm a successful connection
#try:
#    client.admin.command('ping')
#    print("Pinged your deployment. You successfully connected to MongoDB!")
#except Exception as e:
#    print(e)

ROAR_DB = client['ROAR-DB']
Rieske_Collection = ROAR_DB['Rieske']

fasta_one_by_one(protein_folder, hmm_file, e_threshold, Rieske_Collection)
