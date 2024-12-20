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



def is_matching_product(product_name):
    """
    Check if the product name matches any keyword in the filter_products list
    (case-insensitive and partial matching).
    
    Args:
        product_name (str): The product name from GenBank annotation.
        filter_products (list): List of keywords to filter by (e.g., ["transposon", "transposase"]).
    
    Returns:
        bool: True if a match is found, False otherwise.
    """
    filter_products = ["transposon","transposase","regulator","transcriptional regulator","TetR", "transporter",
    "insertion sequence", "reductase", "ferredoxin", "subunit beta", "ferredoxin subunit", "binding protein", 
    "subunit alpha", "small subunit", "electron transfer component"]


    product_name = product_name.lower()  # Convert to lowercase for case-insensitive comparison
    for term in filter_products:
        if term.lower() in product_name:  # Partial matching (substring)
            return True
    return False


def get_sequence_type(genbank_file):
    """
    Extract whether the sequence is from a plasmid or genome based on GenBank annotations.

    Args:
        genbank_file (str): Path to the GenBank file.

    Returns:
        str: A message indicating whether the sequence is a plasmid, genomic DNA, or unclassified.
    """
    sequence_type = "Unclassified"
    
    try:
        with open(genbank_file, "r") as file:
            for record in SeqIO.parse(file, "genbank"):
                # Check 'mol_type' in 'source' feature for genomic DNA or plasmid info
                for feature in record.features:
                    if feature.type == "source":
                        mol_type = feature.qualifiers.get("mol_type", [""])[0].lower()
                        # Check if mol_type contains 'genomic' or 'plasmid'
                        if "genomic" in mol_type:
                            sequence_type = "Genome"
                        elif "plasmid" in mol_type:
                            sequence_type = "Plasmid"
                
                # Alternatively, check the 'definition' line for sequence info
                if 'plasmid' in record.description.lower():
                    sequence_type = "Plasmid"
                elif 'genome' in record.description.lower():
                    sequence_type = "Genome"
    
    except Exception as e:
        print(f"Error reading GenBank file: {e}")
    
    return sequence_type

def get_total_sequence_length(genbank_file):
    """
    Calculate the total length of all sequences in a GenBank file.

    Args:
        genbank_file (str): Path to the GenBank file.

    Returns:
        int: The total length of all sequences.
    """
    total_length = 0
    
    try:
        with open(genbank_file, "r") as file:
            for record in SeqIO.parse(file, "genbank"):
                total_length += len(record.seq)
    except Exception as e:
        print(f"Error reading GenBank file: {e}")
    
    return total_length

# Belirli bir sekans aralığını indirirken uzunluğu kontrol etme
def download_genbank_range_with_check(accession_id, output_file):


    # GenBank aralığını indir
    try:
        with Entrez.efetch(
            db="nucleotide",
            id=accession_id,
            rettype="gb",
            retmode="text",
        ) as handle:
            record = SeqIO.read(handle, "genbank")
            with open(output_file, "w") as output_handle:
                SeqIO.write(record, output_handle, "genbank")
            print(f"{accession_id} aralığı '{output_file}' olarak kaydedildi.")
    except Exception as e:
        print(f"Hata oluştu: {e}")

def calculate_overlap_percentage(range1, range2):
    """
    İki sayı aralığı arasındaki yüzdelik örtüşmeyi hesaplar.
    
    Args:
        range1 (tuple): İlk aralık (başlangıç, bitiş)
        range2 (tuple): İkinci aralık (başlangıç, bitiş)
        
    Returns:
        float: Örtüşme yüzdesi (0-100 arası)
    """
    start1, end1 = range1
    start2, end2 = range2
    
    # Ortak aralığın başlangıç ve bitiş noktalarını bul
    overlap_start = max(start1, start2)
    overlap_end = min(end1, end2)
    
    # Ortak aralık uzunluğunu hesapla
    overlap_length = max(0, overlap_end - overlap_start)
    
    # Birleşik aralık uzunluğunu hesapla
    union_start = min(start1, start2)
    union_end = max(end1, end2)
    union_length = union_end - union_start
    
    # Örtüşme yüzdesini hesapla
    overlap_percentage = (overlap_length / union_length) * 100 if union_length > 0 else 0
    
    return overlap_percentage

def gbk_analysis(genbank_file, target_start, target_end, extention):

    limit_end = target_end + extention
    limit_start = target_start - extention
    if limit_start < 0:
        limit_start = 0

    seq_length = get_total_sequence_length(genbank_file)
    if seq_length is None:
        print("Sekans uzunluğu alınamadı, işlem iptal edildi.")
        return

    # Eğer end değeri toplam uzunluğu aşarsa, onu maksimum uzunluğa ayarla
    if limit_end > seq_length:
        print(f"Uyarı: 'end' değeri ({limit_end}) sekans uzunluğunu ({seq_length}) aşıyor. Maksimum uzunluk alınacak.")
        limit_end = seq_length


    # Sonuçları saklamak için iki liste
    flank = []

    # GenBank dosyasını oku
    with open(genbank_file, "r") as file:
        for record in SeqIO.parse(file, "genbank"):
            # Tüm özellikleri (features) sırayla işle
            for feature in record.features:
                # Özellik bir "gene" veya "CDS" ise işleme al
                if feature.type in ["gene", "CDS"]:
                    # Konum bilgilerini al
                    start = int(feature.location.start)
                    end = int(feature.location.end)

                    if start < limit_start:
                        continue
                    if end > limit_end:
                        continue

                    # Ürün adı varsa al
                    product_name = feature.qualifiers.get("product", ["Unknown"])[0]
                    if is_matching_product(product_name):
                        # Flank kontrolü

                        if calculate_overlap_percentage((target_start, target_end),(start, end)) > 80:
                            continue

                        if start < target_start:  # Sol flanka girenler
                            flank.append({
                                "product": product_name,
                                "start": start,
                                "end": end,
                                "relative": target_start - end,
                                "type": "left"
                            })
                        elif start >= target_start:  # Sağ flanka girenler
                            flank.append({
                                "product": product_name,
                                "start": start,
                                "end": end,
                                "relative": start - target_end,
                                "type": "right"
                            })

    # Sonuçları yazdır
    '''
    print(target_start, target_end)
    print("Left Flank (Start < Target Start):")
    for annotation in left_flank:
        print(f"Product: {annotation['product']}, Start: {annotation['start']}, End: {annotation['end']}, Relative: {annotation['relative']}")

    print("\nRight Flank (End > Target End):")
    for annotation in right_flank:
        print(f"Product: {annotation['product']}, Start: {annotation['start']}, End: {annotation['end']}, Relative: {annotation['relative']}")
    '''
    return flank






load_dotenv()
MONGODB_URI = os.environ['MONGODB_URI']
client = MongoClient(MONGODB_URI)
ROAR_DB = client['ROAR-DB']
Rieske_Collection = ROAR_DB['Rieske']
Relative_Collection = ROAR_DB['RelativeGenes']

total_documents = Rieske_Collection.count_documents({})
processed = 0
os.makedirs('gbk_files', exist_ok=True)

extention = 10000

batch_size = 100
for i in range(0, total_documents, batch_size):
    cursor = Rieske_Collection.find().skip(i).limit(batch_size)
    for document in cursor:
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
                start = int(document['location_in_seq']['start'])
                end = int(document['location_in_seq']['end'])
                gbk_file = f"gbk_files/{nuc_id}.gbk"
                if not os.path.isfile(gbk_file):
                    download_genbank_range_with_check(nuc_id, gbk_file)
                

                if not Relative_Collection.find_one({'_id': document['_id']}):
                    flank = gbk_analysis(gbk_file, start, end, extention)
                    type_seq = get_sequence_type(gbk_file)
                    result = Relative_Collection.insert_one(
                        {"_id": document['_id'], "nuc_id": nuc_id, "flank": flank, "seq_type": type_seq}
                    )   

