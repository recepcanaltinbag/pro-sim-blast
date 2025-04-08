import random
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio.Phylo import draw_ascii
from Bio.Phylo import write

# 1. FASTA dosyasını okuyun
fasta_file = "trees/combinedv2.aln"  # FASTA dosyasının adı
alignment = AlignIO.read(fasta_file, "clustal")

excluded_species = ["TutE", "2_205_IsoMO", "1_113_CdnD"]  # Çıkarmak istediğiniz tür adı
for excluded in excluded_species:
    alignment = MultipleSeqAlignment(
        [record for record in alignment if record.id != excluded]
    )

# 2. Bootstrap için yeniden örnekleme işlevi
def resample_alignment(alignment):
    num_positions = alignment.get_alignment_length()  # Alignment sütun sayısı
    indices = [random.randint(0, num_positions - 1) for _ in range(num_positions)]
    resampled_seqs = [
        SeqRecord("".join(record.seq[i] for i in indices), id=record.id) for record in alignment
    ]
    return MultipleSeqAlignment(resampled_seqs)

# 3. Bootstrap ağacı oluştur ve destek değerlerini hesapla
bootstrap_replicates = 100  # Bootstrap tekrarı sayısı
calculator = DistanceCalculator("identity")
constructor = DistanceTreeConstructor()

bootstrap_trees = []
for _ in range(bootstrap_replicates):
    resampled_alignment = resample_alignment(alignment)
    distance_matrix = calculator.get_distance(resampled_alignment)
    tree = constructor.nj(distance_matrix)
    bootstrap_trees.append(tree)

# 4. Destek değerlerini hesaplama
def calculate_support(nj_tree, bootstrap_trees):
    for clade in nj_tree.get_nonterminals():  # İç kladları dolaş
        clade_bootstrap_support = sum(
            1 for tree in bootstrap_trees if clade in tree.find_clades()
        )
        clade.confidence = clade_bootstrap_support / bootstrap_replicates * 100  # Destek yüzdesi
    return nj_tree

# Orijinal ağacı oluştur
original_distance_matrix = calculator.get_distance(alignment)
nj_tree = constructor.nj(original_distance_matrix)

# Destek değerlerini orijinal ağaca ekle
nj_tree_with_support = calculate_support(nj_tree, bootstrap_trees)

# 5. Outgroup belirleyin (opsiyonel)
#outgroup_name = "1_101_OxoO"  # Outgroup olarak kullanmak istediğiniz tür adı
#if outgroup_name in [term.name for term in nj_tree_with_support.get_terminals()]:
#    nj_tree_with_support.root_with_outgroup(outgroup_name)

# 6. Bootstrap destek değerlerini içeren ağacı Newick formatına kaydedin
tree_file = "trees/nj_tree_with_support.newick"
with open(tree_file, "w") as output_file:
    write(nj_tree_with_support, output_file, format="newick")

print(f"Ağaç başarıyla {tree_file} dosyasına yazdırıldı.")

# 7. ASCII formatında ağacı da yazdırın
print("\nNeighbor-Joining Tree with Bootstrap Support:")
draw_ascii(nj_tree_with_support)
