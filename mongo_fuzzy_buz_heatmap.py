import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

# Klasördeki tüm CSV dosyalarını toplamak
folder_path = "out_Relative_Info"  # Çalışma dizinini belirtin
csv_files = [f for f in os.listdir(folder_path) if f.endswith(".csv")]


# Grupları saklamak için dictionary oluştur
grouped_data = {}

for file in csv_files:
    # Dosya ismini `_`'ye göre split edip grup ismini al
    parts = file.split("_")
    if len(parts) < 3 or not parts[2].endswith(".csv"):
        continue  # Dosya formatı uygun değilse atla
    group_name = parts[0]
    a_value = parts[2].replace(".csv", "")
    
    # CSV dosyasını oku
    file_path = os.path.join(folder_path, file)
    df = pd.read_csv(file_path)
    
    # "G_" ile başlayan grupları filtrele
    filtered_df = df[df["Group"].str.startswith("G_")]
    
    # Normalize et (Min-Max normalization)
    if not filtered_df.empty:
        filtered_df["Score"] = np.log10(filtered_df["Score"] + 1e-6)  # Küçük bir epsilon eklenir
        min_score = filtered_df["Score"].min()
        max_score = filtered_df["Score"].max()
        filtered_df["Score"] = (filtered_df["Score"] - min_score) / (max_score - min_score)
    
    # Grup adı için dictionary oluştur
    if group_name not in grouped_data:
        grouped_data[group_name] = {}
    
    # a_value ve Group-Score çiftlerini sakla
    for _, row in filtered_df.iterrows():
        group = row["Group"]
        score = row["Score"]
        if group not in grouped_data[group_name]:
            grouped_data[group_name][group] = {}
        grouped_data[group_name][group][a_value] = score

# Her grup için ayrı ayrı ısı haritası oluştur ve kaydet
output_folder = "heatmaps"
os.makedirs(output_folder, exist_ok=True)

for group_name, group_data in grouped_data.items():
    # DataFrame'e dönüştür
    heatmap_data = pd.DataFrame(group_data).fillna(0)  # Eksik verileri 0 ile doldur
    heatmap_data = heatmap_data[sorted(heatmap_data.columns)]
    # Isı haritası oluştur
    plt.figure(figsize=(10, 8))
    sns.heatmap(heatmap_data, annot=True, fmt=".2f", cmap="vlag")
    plt.title(f"Heatmap for Group {group_name}")
    plt.xlabel("Group")
    plt.ylabel("ROs")
    plt.yticks(rotation=0)   # Y ekseni etiketlerini yere paralel yap
    plt.tight_layout()
    
    # Dosyayı kaydet
    output_file = os.path.join(output_folder, f"heatmap_{group_name}.png")
    plt.savefig(output_file)
    plt.close()

print(f"Heatmaps saved in {output_folder}")

for group_name, group_data in grouped_data.items():
    # DataFrame'e dönüştür
    heatmap_data = pd.DataFrame(group_data).fillna(0)  # Eksik verileri 0 ile doldur
    heatmap_data = heatmap_data[sorted(heatmap_data.columns)]
    
    # Isı haritası oluştur (clustermap ile)
    cluster = sns.clustermap(
        heatmap_data,
        cmap="vlag",
        annot=True,
        fmt=".2f",
        figsize=(10, 8),
        row_cluster=True,  # Satırları clusterla
        col_cluster=False,  # Sütunları clusterlama
        dendrogram_ratio=(.1, .1),  # Dendrogram boyutlarını ayarla
        cbar_pos=(0.02, 0.8, 0.03, 0.18),  # Colorbar pozisyonu
    )
    
    # Başlık ve eksen ayarları
    cluster.ax_heatmap.set_title(f"Clustered Heatmap for Group {group_name}")
    cluster.ax_heatmap.set_xlabel("Group")
    cluster.ax_heatmap.set_ylabel("ROs")
    
    # Y eksen etiketlerini yere paralel yap
    cluster.ax_heatmap.yaxis.set_tick_params(labelrotation=0)
    cluster.fig.tight_layout()
    # Dosyayı kaydet
    output_file = os.path.join(output_folder, f"clustered_heatmap_{group_name}.png")
    plt.savefig(output_file)
    plt.close()








