import os
import csv
import tempfile
import numpy as np
from Bio import SeqIO, AlignIO
from Bio.Align.Applications import MuscleCommandline
import random
import subprocess

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
    

    cmd = [
        "muscle",  # muscle5 olmalı
        "-threads", "24",
        "-align", input_fasta,
        "-output", output_fasta,

    ]
    subprocess.run(cmd, check=True)

    #muscle_cline = MuscleCommandline(input=input_fasta, out=output_fasta, clwstrict=True)
    #stdout, stderr = muscle_cline()
    #print(stdout)

def average_pairwise_distance(alignment_file):
    alignment = AlignIO.read(alignment_file, "fasta")
    n = len(alignment)
    distances = []
    for i in range(n):
        for j in range(i+1, n):
            seq1 = str(alignment[i].seq)
            seq2 = str(alignment[j].seq)
            pdist = calculate_pdistance(seq1, seq2)

            distances.append(pdist)
    return np.mean(distances) if distances else 0, distances

def main_analysis_function(fasta_path, threshold=1500):
    main_seqs = []
    regulator_seqs = []

    try:
        records = list(SeqIO.parse(fasta_path, "fasta"))
        i = 0
        while i < len(records):
            # MAIN_GENE ve REGULATOR çiftlerini al
            if i + 1 >= len(records):
                break  # REGULATOR eksikse döngüyü bitir

            main_record = records[i]
            reg_record = records[i + 1]

            main_seq = str(main_record.seq)
            reg_seq = str(reg_record.seq)

            # Eğer main gene 300 bp'den kısa ise, ikisini de atla
            if len(main_seq) < 280:
                i += 2
                continue

            # MAIN_GENE'i ekle
            gene_type, start, end, contig = parse_fasta_header(main_record.description)
            main_seqs.append({
                "id": main_record.id,
                "type": gene_type,
                "start": start,
                "end": end,
                "contig": contig,
                "seq": main_seq
            })

            # REGULATOR'u ekle
            gene_type, start, end, contig = parse_fasta_header(reg_record.description)
            regulator_seqs.append({
                "id": reg_record.id,
                "type": gene_type,
                "start": start,
                "end": end,
                "contig": contig,
                "seq": reg_seq
            })

            i += 2

    except Exception as e:
        print(f"{fasta_path} okunurken hata oluştu: {e}")
        return None

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

        if closest_reg and min_dist <= threshold:
            filtered_main.append(main)
            filtered_regulator.append(closest_reg)

    count_pairs = len(filtered_main)
    print(count_pairs)
    if count_pairs == 0:
        return None
    if count_pairs == 1:
        return None
    if count_pairs > 50:
        indices = random.sample(range(count_pairs), 50)
        filtered_main = [filtered_main[i] for i in indices]
        filtered_regulator = [filtered_regulator[i] for i in indices]
        print(f"Çok fazla eşleşme vardı ({count_pairs}), rastgele 500 tanesi alındı.")
    else:
        print(f"Eşleşme sayısı: {count_pairs}")
    
    with tempfile.TemporaryDirectory() as tmpdir:
        main_fasta = os.path.join(tmpdir, "main.fasta")
        reg_fasta = os.path.join(tmpdir, "reg.fasta")
        main_aln = os.path.join(tmpdir, "main.clw")
        reg_aln = os.path.join(tmpdir, "reg.clw")

        with open(main_fasta, "w") as f_main, open(reg_fasta, "w") as f_reg:
            for i, main in enumerate(filtered_main):
                f_main.write(f">main_{i}\n{main['seq']}\n")
            for i, reg in enumerate(filtered_regulator):
                f_reg.write(f">reg_{i}\n{reg['seq']}\n")
        try:
            run_muscle(main_fasta, main_aln)
            run_muscle(reg_fasta, reg_aln)

            main_div, main_all_div = average_pairwise_distance(main_aln)
            reg_div, reg_all_div = average_pairwise_distance(reg_aln)
        except:
            print('No')
    return {
        "pair_count": count_pairs,
        "main_diversity": main_div,
        "reg_diversity": reg_div,
        "all_main_diversity": main_all_div,
        "all_reg_diversity": reg_all_div
    }

def analyze_filtered_fastas(base_dir="Filtered_Proteins", threshold=1500):
    summary_file = "Enzyme_Regulator_summary.csv"
    overall_stats_file = "Enzyme_Regulator_overall_summary.txt"
    per_file_diversity_file = "Per_File_Diversity_Values.csv"

    results = []
    all_main_diversities = []
    all_reg_diversities = []

    fasta_files = [f for f in os.listdir(base_dir) if f.endswith(".fasta")]

    for fasta_name in fasta_files:
        fasta_path = os.path.join(base_dir, fasta_name)
        print(f"Analiz başlatılıyor: {fasta_name}")
        result = main_analysis_function(fasta_path, threshold=threshold)
        if result:
            result["filename"] = fasta_name
            results.append(result)
            all_main_diversities.append(result["main_diversity"])
            all_reg_diversities.append(result["reg_diversity"])

    # Toplu özet CSV
    with open(summary_file, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["filename", "pair_count", "main_diversity", "reg_diversity", "all_main_diversity", "all_reg_diversity"])
        writer.writeheader()
        for row in results:
            writer.writerow(row)

    # Her dosya için çeşitlilik değerlerini ayrı CSV dosyasına yaz
    with open(per_file_diversity_file, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["filename", "main_diversity", "reg_diversity"])
        for r in results:
            the_main_values = r["all_main_diversity"]
            the_reg_values = r["all_reg_diversity"]
            filename = r["filename"]
            for main_val, reg_val in zip(the_main_values, the_reg_values):
                writer.writerow([filename, main_val, reg_val])

    # Genel istatistikler
    with open(overall_stats_file, "w") as f:
        f.write(f"Analiz edilen dosya sayısı: {len(results)}\n")
        f.write(f"Toplam eşleşen çift sayısı: {sum(r['pair_count'] for r in results)}\n")
        if all_main_diversities:
            avg_main = sum(all_main_diversities) / len(all_main_diversities)
            avg_reg = sum(all_reg_diversities) / len(all_reg_diversities)
            f.write(f"Main ortalama çeşitlilik: {avg_main:.4f}\n")
            f.write(f"Reg ortalama çeşitlilik: {avg_reg:.4f}\n")
            if avg_main > avg_reg:
                f.write("Genel olarak main genler daha çeşitli.\n")
            elif avg_main < avg_reg:
                f.write("Genel olarak regülatör genler daha çeşitli.\n")
            else:
                f.write("Her iki grup da benzer çeşitlilikte.\n")
        else:
            f.write("Yeterli veri yok.\n")


if __name__ == "__main__":
    analyze_filtered_fastas()
