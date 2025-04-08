
import os
from Bio.Align.Applications import ClustalwCommandline
from Bio import Phylo
from Bio.SeqIO import parse, write
import matplotlib.pyplot as plt
from ete3 import Tree, TreeStyle, NodeStyle
from Bio import AlignIO
# ALN formatındaki hizalamayı oku
input_aln = "trees/combined.aln"
output_aln = "trees/combinedv2.aln"

alignment = AlignIO.read(input_aln, "clustal")

# Dizilerdeki 'N' karakterlerini '-' ile değiştir
for record in alignment:
    record.seq = record.seq.tomutable()  # Mutable hale getir
    record.seq = record.seq.replace('N', '-')  # 'N' karakterlerini '-' ile değiştir

# Değiştirilmiş hizalamayı yeni dosyaya kaydet
with open(output_aln, "w") as output_handle:
    AlignIO.write(alignment, output_handle, "clustal")

print(f"Değiştirilmiş hizalama {output_aln} dosyasına kaydedildi.")


# Orijinal veriyi hizalanmış FASTA formatına dönüştürmek için
output_aln = "trees/combinedv2.aln"
output_fasta = "trees/rastInv3.fasta"  # Yukarıdaki veri dosyasını kaydedin

# Aln dosyasını FASTA formatına dönüştür
from Bio import AlignIO
AlignIO.convert(output_aln, "clustal", output_fasta, "phylip")
print(f"Hizalama FASTA formatında kaydedildi: {output_fasta}")

input()





# Klasör ve ClustalW yolunu tanımlayın
input_folder = "./Proteins_ROs_alphaSubUnits"  # FASTA dosyalarının olduğu klasör
clustalw_path = "clustalw"  # ClustalW'nin sistemdeki yolu
output_folder = "./trees"  # Çıktıların kaydedileceği klasör
os.makedirs(output_folder, exist_ok=True)

# Birleştirilmiş FASTA dosyasının yolu
combined_fasta = os.path.join(output_folder, "combined.fasta")
'''
# Tüm FASTA dosyalarını birleştir
with open(combined_fasta, "w") as outfile:
    for fasta_file in os.listdir(input_folder):
        if fasta_file.endswith(".fasta"):
            fasta_path = os.path.join(input_folder, fasta_file)
            # Dosya adını etikete dönüştür ve dizileri yaz
            for record in parse(fasta_path, "fasta"):
                record.id =  '_'.join(os.path.basename(fasta_file).split('_')[:3]) 
                record.description = ""
                write(record, outfile, "fasta")
'''
print(f"Tüm diziler {combined_fasta} dosyasına birleştirildi.")

# ClustalW ile hizalama yap
aligned_fasta = os.path.join(output_folder, "combined.fasta")
clustalw_cline = ClustalwCommandline(
    clustalw_path,
    infile=combined_fasta,
    bootstrap=1000,  # Bootstrap tekrar sayısı
    output="PHYLIP"  # Ağaç oluştururken PHYLIP formatını kullanır
)
stdout, stderr = clustalw_cline()


# Filogenetik ağaç dosyasını BioPython ile oku ve dönüştür
tree_file = os.path.join(output_folder, "combined.dnd")
# Ağaç dosyasını oku
tree = Phylo.read(tree_file, "newick")

# Ağacı çiz (örneğin, yuvarlak bir ağaç)
Phylo.draw(tree, branch_labels=lambda c: f"{c.confidence:.1f}" if c.confidence else "")