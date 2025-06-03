from Bio import SeqIO
from Bio.Align.Applications import MuscleCommandline
from Bio import AlignIO
import tempfile
import os
import numpy as np

def parse_fasta_header(header):
    parts = header.split('|')
    gene_type = parts[0].strip('>').strip()
    start_end = parts[-2].strip()
    contig = parts[-1].strip()
    start, end = map(int, start_end.split('-'))
    return gene_type, start, end, contig

def calculate_pdistance(seq1, seq2):
    mismatches = sum(a != b for a, b in zip(seq1, seq2))
    length = min(len(seq1), len(seq2))
    return mismatches / length if length > 0 else 0

def run_muscle(input_fasta, output_fasta):
    muscle_cline = MuscleCommandline(input=input_fasta, out=output_fasta, clwstrict=True)
    stdout, stderr = muscle_cline()

def average_pairwise_distance(alignment_file):
    alignment = AlignIO.read(alignment_file, "clustal")
    n = len(alignment)
    distances = []
    for i in range(n):
        for j in range(i+1, n):
            seq1 = str(alignment[i].seq)
            seq2 = str(alignment[j].seq)
            pdist = calculate_pdistance(seq1, seq2)
            distances.append(pdist)
    return np.mean(distances) if distances else 0

def main(distance_threshold=1000):
    fasta_path = "proteins_output.fasta"

    main_seqs = []
    regulator_seqs = []

    # Fasta dosyasını oku
    for record in SeqIO.parse(fasta_path, "fasta"):
        gene_type, start, end, contig = parse_fasta_header(record.description)
        seq = str(record.seq)

        entry = {
            "id": record.id,
            "type": gene_type,
            "start": start,
            "end": end,
            "contig": contig,
            "seq": seq
        }

        if gene_type == "MAIN_GENE":
            main_seqs.append(entry)
        elif gene_type == "REGULATOR":
            regulator_seqs.append(entry)

    filtered_main = []
    filtered_regulator = []

    for main in main_seqs:
        contig = main['contig']
        main_mid = (main['start'] + main['end']) // 2
        regulators_same_contig = [r for r in regulator_seqs if r['contig'] == contig]

        min_dist = None
        closest_reg = None
        for reg in regulators_same_contig:
            reg_mid = (reg['start'] + reg['end']) // 2
            dist = abs(main_mid - reg_mid)
            if min_dist is None or dist < min_dist:
                min_dist = dist
                closest_reg = reg

        if closest_reg and min_dist <= distance_threshold:
            filtered_main.append(main)
            filtered_regulator.append(closest_reg)

    count_pairs = len(filtered_main)

    if count_pairs == 0:
        print("Yeterli çift bulunamadı, işlem sonlandırılıyor.")
        return

    print(f"Filtrelenmiş çift sayısı (mesafe ≤ {distance_threshold} bp): {count_pairs}")

    with open("filtered_main.fasta", "w") as f_main, open("filtered_regulator.fasta", "w") as f_reg:
        for i, main in enumerate(filtered_main):
            f_main.write(f">main_{i}\n{main['seq']}\n")
        for i, reg in enumerate(filtered_regulator):
            f_reg.write(f">reg_{i}\n{reg['seq']}\n")

    with tempfile.TemporaryDirectory() as tmpdir:
        main_aln = os.path.join(tmpdir, "main_aln.clw")
        reg_aln = os.path.join(tmpdir, "reg_aln.clw")

        run_muscle("filtered_main.fasta", main_aln)
        run_muscle("filtered_regulator.fasta", reg_aln)

        main_diversity = average_pairwise_distance(main_aln)
        reg_diversity = average_pairwise_distance(reg_aln)

    print(f"Main gen ortalama çeşitlilik (p-distance): {main_diversity:.4f}")
    print(f"Regülatör ortalama çeşitlilik (p-distance): {reg_diversity:.4f}")

    if main_diversity > reg_diversity:
        print("Main genler regulatorlerden daha çeşitli.")
    elif reg_diversity > main_diversity:
        print("Regülatörler main genlerden daha çeşitli.")
    else:
        print("Her iki grup da benzer çeşitlilik gösteriyor.")

if __name__ == "__main__":
    # İstersen threshold'u buradan değiştirebilirsin
    main(distance_threshold=1000)
