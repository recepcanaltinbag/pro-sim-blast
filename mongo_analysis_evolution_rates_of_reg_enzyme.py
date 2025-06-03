from rapidfuzz import process
import pandas as pd
from pymongo.mongo_client import MongoClient
from dotenv import load_dotenv
import os
from Bio import SeqIO
from pathlib import Path



load_dotenv()
MONGODB_URI = os.environ['MONGODB_URI']
client = MongoClient(MONGODB_URI)
ROAR_DB = client['ROAR-DB']
Rieske_Collection = ROAR_DB['Rieske']
Relative_Collection = ROAR_DB['RelativeGenes']


# 1. Rieske'den filtrele
rieske_filtered = list(Rieske_Collection.find({
    "Group and ID": "102",
    "E-value": {"$lte": 1.3e-14},
    "location_in_seq": {"$ne": None},
    "sequence": {"$ne": None}
}))

print(f"{len(rieske_filtered)} kayıt bulundu.\n")

# 2. Dizinler
gbk_dir = Path("gbk_files")
seq_cache = {}

output_path = Path("proteins_output.fasta")
with open(output_path, "w") as fasta_out:
    for entry in rieske_filtered:
        loc = entry.get("location_in_seq", {})
        nucleotide_id = loc.get("nucleotide_id")
        main_start = loc.get("start")
        main_end = loc.get("end")
        aa_seq = entry.get("sequence")

        if not nucleotide_id or not aa_seq or main_start is None:
            continue

        # Relative'lerden ilgili nucleotide_id'yi al
        relative_matches = list(Relative_Collection.find({
            "nuc_id": nucleotide_id
        }))
        if not relative_matches:
            continue

        # gbk yükle
        if nucleotide_id not in seq_cache:
            gbk_path = gbk_dir / f"{nucleotide_id}.gbk"
            if not gbk_path.exists():
                print(f"GBK bulunamadı: {gbk_path}")
                continue
            try:
                seq_cache[nucleotide_id] = SeqIO.read(str(gbk_path), "genbank")
            except Exception as e:
                print(f"GBK hatası: {nucleotide_id} -> {e}")
                continue

        record = seq_cache[nucleotide_id]

        # En yakın regulator'ı bul
        closest_regulator = None
        closest_translation = ""
        min_distance = float('inf')

        for match in relative_matches:
            flanks = match.get("flank", [])
            for flank in flanks:
                product = flank.get("product", "").lower()
                if "regulator" in product:
                    flank_start = flank.get("start")
                    flank_end = flank.get("end")
                    if flank_start is None or flank_end is None:
                        continue

                    # GBK içinden translation'ı bul
                    for feature in record.features:
                        if feature.type == "CDS":
                            f_start = int(feature.location.start)
                            f_end = int(feature.location.end)
                            if f_start == flank_start and f_end == flank_end:
                                raw_product = feature.qualifiers.get("product", ["unknown"])[0]
                                translation = feature.qualifiers.get("translation", [""])[0]
                                if translation:
                                    distance = abs(main_start - flank_start)
                                    if distance < min_distance:
                                        min_distance = distance
                                        closest_regulator = {
                                            "start": flank_start,
                                            "end": flank_end,
                                            "product": raw_product
                                        }
                                        closest_translation = translation
                                break

        if not closest_regulator or not closest_translation:
            continue  # regulator yoksa main'i de yazma

        # FASTA yaz: önce main sonra regulator
        header_main = f">MAIN_GENE | {entry.get('Protein ID', 'unknown')} | {main_start}-{main_end} | {nucleotide_id}"
        fasta_out.write(f"{header_main}\n{aa_seq}\n")

        header_reg = f">REGULATOR | {closest_regulator['product']} | {closest_regulator['start']}-{closest_regulator['end']} | {nucleotide_id}"
        fasta_out.write(f"{header_reg}\n{closest_translation}\n")

print(f"İşlem tamamlandı. Çıktı dosyası: {output_path}")