import os

# Klasör yolunu belirt
input_folder = "out_proteins"  # Burayı değiştir
output_folder = "out_pro_combined"

# Çıktı klasörünü oluştur
os.makedirs(output_folder, exist_ok=True)


# Her bir dosyayı işle
for file_name in os.listdir(input_folder):
    if file_name.endswith(".fasta"):  # Sadece FASTA dosyalarını işle
        file_path = os.path.join(input_folder, file_name)
        seq_id_base = "_".join(file_name.split("_")[:3])  # İlk iki `_` ile ayrılmış kısmı al
        number_s = file_name.split("_")[3].split('.')[0]
        # Seq_ID'ye göre çıktı dosyası oluştur
        output_file = os.path.join(output_folder, f"{seq_id_base}.fasta")
        
        with open(file_path, "r") as infile, open(output_file, "a") as outfile:
            for line in infile:
                if line.startswith(">"):  # Header satırı
                    original_seq_id = line.split(" ")[0][1:]  # Orijinal ID'nin ilk elemanı
                    original_seq_id = original_seq_id.split("|")[0]  # Orijinal ID'nin ilk elemanı
                    new_seq_id = f"{seq_id_base}_{number_s}_{original_seq_id}"  # Yeni ID
                    outfile.write(f">{new_seq_id}\n")  # Yeni header yaz
                else:
                    # Boş satırları kontrol et ve yazma
                    if line.strip():  # line.strip() boşlukları temizler
                        outfile.write(line)  # Sekans satırını ekle

print("Çıktılar başarıyla oluşturuldu.")
