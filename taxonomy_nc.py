import pandas as pd
from Bio import Entrez

# NCBI Taxonomy API için e-posta ayarı
Entrez.email = "your_email@example.com"

def fetch_taxonomy_info(tax_id):
    """NCBI'den Taxonomy bilgilerini sorgular."""
    try:
        handle = Entrez.efetch(db="taxonomy", id=str(tax_id), retmode="xml")
        records = Entrez.read(handle)
        if records:
            record = records[0]
            lineage = record.get("Lineage", "")
            lineage_list = lineage.split("; ")
            return {
                "Domain": lineage_list[0] if len(lineage_list) > 0 else "",
                "Phylum": lineage_list[1] if len(lineage_list) > 1 else "",
                "Class": lineage_list[2] if len(lineage_list) > 2 else ""
            }
    except Exception as e:
        print(f"Error fetching Tax ID {tax_id}: {e}")
        return {"Domain": "", "Phylum": "", "Class": ""}

# 1. Girdi dosyasını oku
input_file = "protein-matching-PF00355.tsv"
output_file = "protein_matchin_with_taxonomy.tsv"

df = pd.read_csv(input_file, sep="\t")

# 2. Her bir Tax ID için taxonomy bilgilerini topla
taxonomy_data = {}
for tax_id in df["Tax ID"].unique():
    taxonomy_data[tax_id] = fetch_taxonomy_info(tax_id)

# 3. Taxonomy bilgilerini DataFrame'e ekle
df["Domain"] = df["Tax ID"].map(lambda x: taxonomy_data[x]["Domain"])
df["Phylum"] = df["Tax ID"].map(lambda x: taxonomy_data[x]["Phylum"])
df["Class"] = df["Tax ID"].map(lambda x: taxonomy_data[x]["Class"])

# 4. Yeni dosyayı kaydet
df.to_csv(output_file, sep="\t", index=False)

print(f"Yeni dosya oluşturuldu: {output_file}")

