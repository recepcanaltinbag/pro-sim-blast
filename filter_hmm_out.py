from Bio import SearchIO
import pandas as pd
import re
import os

def parse_filter_hmm(tblout_file = 'pfamPart.out', e_threshold = 1e-60, output_file="parsed_results.csv"):
    # HMMER3 tblout dosyasını aç
    

    # Verileri saklamak için liste
    data = []

    #  hmmsearch --tblout pfamPart.out RieskeDB.hmm part_pfam.fasta
    # PARSING THE HMMER TBL_OUT OUTPUT 
    with open(tblout_file, "r") as file:
        for line in file:
            # Yorum satırlarını atla
            if line.startswith("#"):
                continue
            
            # Satırdaki bilgileri düzenli ifadeyle ayır
            match = re.match(r"^(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)", line)
            if match:
                target_name = match.group(1)
                accession = match.group(2)
                query_name = match.group(3)
                query_accession = match.group(4)
                e_value = float(match.group(5))
                score = float(match.group(6))

                # Verileri listeye ekle
                data.append([target_name, accession, query_name, query_accession, e_value, score])

    # DataFrame oluştur
    columns = ["target_name", "accession", "query_name", "query_accession", "e_value", "score"]
    df = pd.DataFrame(data, columns=columns)
    df = df.sort_values(by="e_value", ascending=True)
    # Her query_name ve target_name çifti için ilk sonucu seç
    #result_df = df.sort_values(by="e_value").drop_duplicates(subset=["target_name", "query_name"], keep="first")

    # E-value filtrele (e-10'dan daha küçük olanları al)
    result_df = df[df["e_value"] < e_threshold]

    if result_df.empty:
        print('There is no similar results can be a new one')
        return None
    else:
        first_acc_name = df.iloc[0]["target_name"]
        first_query_Rieske_name = df.iloc[0]["query_name"]
        first_e_value = df.iloc[0]["e_value"]
        first_score = df.iloc[0]["score"]
        group = first_query_Rieske_name.split('_')[0]
        group_and_ID = first_query_Rieske_name.split('_')[1]
        gene_abr = first_query_Rieske_name.split('_')[2]
        uniProtAcc = first_acc_name.split('|')[1]
        print(f"First Accession Name: {first_acc_name}, First Query Rieske Name: {first_query_Rieske_name}, "
            f"First E-value: {first_e_value}, First Score: {first_score}, Group: {group}, "
            f"Group and ID: {group_and_ID}, Gene Abbreviation: {gene_abr}, UniProt Accession: {uniProtAcc}")
        result_dict = {
            "Accession Name": first_acc_name,
            "PredictedCluster": first_query_Rieske_name,
            "E-value": first_e_value,
            "Score": first_score,
            "Group": group,
            "Group and ID": group_and_ID,
            "Gene Abbreviation": gene_abr,
            "UniProt Accession": uniProtAcc
        }
        # Sonuçları CSV dosyasına kaydet
        if not os.path.isfile(output_file):  # Dosya yoksa, başlıkla yaz
            result_df.to_csv(output_file, index=False)
        else:  # Dosya varsa, sadece veri ekle
            result_df.to_csv(output_file, mode="a", index=False, header=False)
        # Çıktıyı görüntüle
        #print(result_df)
        return result_dict






















