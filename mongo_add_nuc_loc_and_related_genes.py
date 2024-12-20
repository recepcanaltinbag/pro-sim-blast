from pymongo.mongo_client import MongoClient
from dotenv import load_dotenv
import os
from pprint import pprint
from f_uniprot_metadata_retrieve import uniprot_id_to_info_dict
from filter_hmm_out import parse_filter_hmm
import subprocess
from Bio import Entrez, SeqIO, SearchIO
import subprocess
import sys
import time

# tblastn komutunu çalıştırmak için bir fonksiyon
def run_tblastn(query_file, db_path, output_file,  num_threads=4):
    try:
        # tblastn komutunu hazırlayın
        command = [
            "tblastn",
            "-query", query_file,            # Protein sorgu dosyası
            "-db", db_path,                  # Nükleotid veri tabanı
            "-out", output_file,             # Çıkış dosyası          # E-value eşiği
            "-outfmt", "5",                  # Tab-delimited format
            "-num_threads", str(num_threads) # Paralel iş parçacığı sayısı
        ]

        # Komutu çalıştırın
        subprocess.run(command, check=True, stdout=subprocess.DEVNULL)
        #print(f"tblastn analizi tamamlandı. Çıkış dosyası: {output_file}")
    except subprocess.CalledProcessError as e:
        print(f"tblastn çalıştırılırken bir hata oluştu: {e}")
    except FileNotFoundError:
        print("tblastn programı bulunamadı. PATH ayarlarınızı kontrol edin.")

def run_blastdb(in_file, output_file):
    try:
        # tblastn komutunu hazırlayın
        command = [
            "makeblastdb",
            "-in", in_file,            # Protein sorgu dosyası
            "-dbtype", "nucl",                  # Nükleotid veri tabanı
            "-out", output_file,             # Çıkış dosyası          # E-value eşiği
        ]

        # Komutu çalıştırın
        subprocess.run(command, check=True, stdout=subprocess.DEVNULL)
        #print(f"tblastn analizi tamamlandı. Çıkış dosyası: {output_file}")
    except subprocess.CalledProcessError as e:
        print(f"tblastn çalıştırılırken bir hata oluştu: {e}")
    except FileNotFoundError:
        print("tblastn programı bulunamadı. PATH ayarlarınızı kontrol edin.")


# NCBI'den nükleotid dizisi indirme
def fetch_nucleotide_sequences(nucleotide_id, output_dir):
    Entrez.email = "your_email@example.com"  # Entrez için e-posta adresinizi girin
    os.makedirs(output_dir, exist_ok=True)
    handle = Entrez.efetch(db="nucleotide", id=nucleotide_id, rettype="fasta", retmode="text")
    file_path = os.path.join(output_dir, f"S.fasta")
    with open(file_path, "w") as file:
        file.write(handle.read())
    #print(f"Nucleotide ID {nucleotide_id} için FASTA dosyası kaydedildi.")
    return file_path


def perform_blast(protein_sequence, file_path, blast_db, output_file):
    os.makedirs(blast_db, exist_ok=True)
    blast_db_out = os.path.join(blast_db, "db")

    run_blastdb(file_path, blast_db_out)
    #os.system(f"makeblastdb -in {file_path} -dbtype nucl -out {blast_db_out}")
    #print(blast_db_out)

    # BLAST çalıştır
    
    protein_file = os.path.join(blast_db, "query_protein.fasta")
    with open(protein_file, "w") as file:
        file.write(f">query_protein\n{protein_sequence}")
    
    run_tblastn(protein_file, blast_db_out, output_file)

    #print("BLAST tamamlandı. Sonuçlar kaydedildi.")

# MongoDB'de güncelleme
def update_mongodb_with_blast_results(my_id, nuc_id, collection, blast_results_file):
    
    # BLAST sonuçlarını analiz et
    blast_records = SearchIO.parse(blast_results_file, "blast-xml")
    for record in blast_records:
        query_id = record.id
        for hit in record.hits:
            subject_id = hit.id
            hsp = hit.hsps[0]
            start, end = hsp.hit_start + 1, hsp.hit_end  # BLAST sonuçları 0-indexlidir
            #print(start, end)
            # MongoDB'de eşleşen belgeyi güncelle
            result = collection.update_one(
                {"_id": my_id},
                {"$set": {"location_in_seq": {"start": start, "end": end, "nucleotide_id": nuc_id}}},
                upsert=True
            )
            #print(f"Protein ID {query_id} için güncelleme tamamlandı.")



load_dotenv()
MONGODB_URI = os.environ['MONGODB_URI']

client = MongoClient(MONGODB_URI)
# Send a ping to confirm a successful connection
#try:
#    client.admin.command('ping')
#    print("Pinged your deployment. You successfully connected to MongoDB!")
#except Exception as e:
#    print(e)

ROAR_DB = client['ROAR-DB']
Rieske_Collection = ROAR_DB['Rieske']
total_documents = Rieske_Collection.count_documents({})
processed = 0


#group_query = {"Group and ID": "505"}




batch_size = 100
for i in range(0, total_documents, batch_size):
    cursor = Rieske_Collection.find().skip(i).limit(batch_size)
    for document in cursor:
#for document in Rieske_Collection.find():
        #print(document['_id'])
        processed += 1
        progress = int((processed / total_documents) * 50)  # 50 birimlik çubuk
        sys.stdout.write(f"\rProcessing: [{'#' * progress}{'.' * (50 - progress)}] {processed}/{total_documents}")
        sys.stdout.flush()

        #print(document['sequence'])
        if 'nucleotide_sequence_ids' not in document:  # Eğer location_in_seq alanı varsa
            #print('Dont exist any Group')
            continue  # Bu dökümanı atla ve bir sonrakine geç
        nuc_ids = document['nucleotide_sequence_ids']

        for nuc_id in nuc_ids:
            if "location_in_seq" in document:  # Eğer location_in_seq alanı varsa
                #print('Already exist')
                continue  # Bu dökümanı atla ve bir sonrakine geç
            #print(nuc_id)
            fasta_file = fetch_nucleotide_sequences(nuc_id, "blast_results_temp")
            perform_blast(document['sequence'], fasta_file, "blast_results_temp", "blast_results_temp/blast_results.xml")
            update_mongodb_with_blast_results(document['_id'], nuc_id, Rieske_Collection, "blast_results_temp/blast_results.xml")






