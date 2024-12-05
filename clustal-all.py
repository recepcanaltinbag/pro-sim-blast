import os
import subprocess


input_folder = "out_pro_combined"  # all proteins families in a single fasta file
output_folder = "aligned"  # alignment outs

if not os.path.exists(output_folder):
    os.makedirs(output_folder)

# Klasördeki .fasta dosyalarını işle
for file_name in os.listdir(input_folder):
    if file_name.endswith(".fasta"):
        file_path = os.path.join(input_folder, file_name)
        base_name = os.path.splitext(file_name)[0]  # Dosya adını uzantıdan ayır
        output_path = os.path.join(output_folder, f"{base_name}.sto")

        # Clustalo komutunu çalıştır
        subprocess.run(["clustalo", "-i", file_path, "-o", output_path, "--outfmt=st"])

print("Hizalamalar tamamlandı.")


# Hizalamaların bulunduğu klasör
alignment_folder = output_folder  # Stockholm formatlı hizalamalar buradan alınacak
output_file = "Rieske_DB_For_HMMER.sto"  # Çıktı dosyası

# Çıktı dosyasını oluştur
with open(output_file, "w") as outfile:
    # Dosyanın başına sadece bir kez Stockholm başlığını ekle
    

    for file_name in os.listdir(alignment_folder):
        if file_name.endswith(".sto"):  # Sadece .sto dosyalarını işle
            file_path = os.path.join(alignment_folder, file_name)
            alignment_name = os.path.splitext(file_name)[0]  # Dosya adını al
            with open(file_path, "r") as infile:
                lines = infile.readlines()

                # İlk satırdaki # STOCKHOLM 1.0 başlığını atla
                if lines[0].startswith("# STOCKHOLM 1.0"):
                    lines = lines[1:]

                # Name annotation ekle
                outfile.write("# STOCKHOLM 1.0\n")
                outfile.write(f"#=GF ID {alignment_name}")
                outfile.writelines(lines)




print(f"Birleştirme tamamlandı. Çıktı: {output_file}")
