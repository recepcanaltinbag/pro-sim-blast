import pandas as pd
import plotly.graph_objects as go
import os
from plotly.colors import n_colors

# CSV dosyasını oku
csv_file = "trees/regulator_having_Rieske.csv"  # CSV dosyanızın adını yazın
df = pd.read_csv(csv_file, delimiter=",")  # Dosyanızın ayracı tab ise `delimiter="\t"` kullanın

# Cluster'ları grup numarasına göre ayır
df["Group"] = df["Cluster"].str.extract(r"^(\d+)_")  # Cluster isimlerinin başındaki sayıyı al

# "Family" sütununu cümlenin ilk kelimesini alarak oluştur
df["Family"] = df["Product"].str.split(r'[-/ ]').str[0]  # Cümlenin ilk kelimesi veya önceki kısmı (family) olarak al

# Çıktıları kaydedeceğiniz dizini oluştur
output_dir = "sankey_outputs"
os.makedirs(output_dir, exist_ok=True)

# Renk paletini tanımla
color_palette = [
    "rgba(31, 119, 180, 0.4)",
    "rgba(255, 127, 14, 0.4)",
    "rgba(44, 160, 44, 0.4)",
    "rgba(214, 39, 40, 0.4)",
    "rgba(148, 103, 189, 0.4)",
    "rgba(140, 86, 75, 0.4)",
    "rgba(227, 119, 194, 0.4)",
    "rgba(127, 127, 127, 0.4)",
    "rgba(188, 189, 34, 0.4)",
    "rgba(23, 190, 207, 0.4)",
    "rgba(31, 119, 180, 0.4)",
    "rgba(255, 127, 14, 0.4)",
    "rgba(44, 160, 44, 0.4)",
    "rgba(214, 39, 40, 0.4)",
    "rgba(148, 103, 189, 0.4)",
    "rgba(140, 86, 75, 0.4)",
    "rgba(227, 119, 194, 0.4)",
    "rgba(127, 127, 127, 0.4)",
    "rgba(188, 189, 34, 0.4)",
    "rgba(23, 190, 207, 0.4)",
    "rgba(31, 119, 180, 0.4)",
    "rgba(255, 127, 14, 0.4)",
    "rgba(44, 160, 44, 0.4)",
    "rgba(214, 39, 40, 0.4)",
    "rgba(148, 103, 189, 0.4)",
    "magenta",
    "rgba(227, 119, 194, 0.4)",
    "rgba(127, 127, 127, 0.4)",
    "rgba(188, 189, 34, 0.4)",
    "rgba(23, 190, 207, 0.4)"
]

# Her grup için Sankey diyagramı oluştur
unique_groups = df["Group"].unique()
threshold_percentage = 0.02  # %5'ten az olanları "Other" olarak grupla
for group in unique_groups:
    group_df = df[df["Group"] == group]  # İlgili grubu seç

    # Sankey için veri hazırlama: Cluster -> Family
    df_counts = group_df.groupby(["Cluster", "Family"]).size().reset_index(name="Value")

    # Ağırlıklı toplam değeri hesapla
    total_value = df_counts["Value"].sum()

    # "Family" kategorilerini %5'ten az olanları "Other" olarak grupla
    df_counts["Family"] = df_counts.apply(
        lambda row: row["Family"] if row["Value"] / total_value >= threshold_percentage else "Other", axis=1
    )

    # Yeniden Family ve Cluster'a göre gruplama
    df_counts = df_counts.groupby(["Cluster", "Family"]).sum().reset_index()

    # Benzersiz kaynaklar, hedefler olarak Cluster ve Family'yi listele
    sources = list(df_counts["Cluster"].unique())
    targets = list(df_counts["Family"].unique())

    # Tüm düğümleri birleştir
    all_nodes = sources + targets

    # İndeksler
    source_indices = df_counts["Cluster"].apply(lambda x: all_nodes.index(x))
    target_indices = df_counts["Family"].apply(lambda x: all_nodes.index(x))
    cluster_colors = {cluster: color_palette[i % len(color_palette)] for i, cluster in enumerate(sources)}

    # Linkler: Cluster -> Family
    links = pd.DataFrame({
        "source": source_indices,
        "target": target_indices,
        "value": df_counts["Value"],
        "color": [cluster_colors.get(cluster, "lightblue") for cluster in df_counts["Cluster"]]  # Color based on source's cluster
    })

    # Her grup için benzersiz bir renk paleti oluştur
    # Cluster renk paleti oluştur
    cluster_colors = {cluster: color_palette[i % len(color_palette)] for i, cluster in enumerate(sources)}
    # Family renk paleti oluştur
    family_colors = {family: color_palette[i % len(color_palette)] for i, family in enumerate(targets)}

    # Node renklerini ve link renklerini ayarla
    node_colors = [
        family_colors.get(node, cluster_colors.get(node, "lightblue"))
        if node in targets else cluster_colors.get(node, "lightblue")
        for node in all_nodes
    ]

    # Link renklerini doğrudan rgba formatında ayarlayın
    link_colors = [node_colors[src].replace("rgb", "rgba") for src in links["source"]]

    # Sankey diyagramını çizdir
    fig = go.Figure(data=[go.Sankey(
        valueformat = ".0f",
        valuesuffix = "TWh",
        # Define nodes
        node=dict(
            pad=15,
            thickness=15,
            line=dict(color="black", width=0.5),
            label=all_nodes,
            color=node_colors
        ),
        # Add links
        link=dict(
            source=links["source"],
            target=links["target"],
            value=links["value"],
            label=links["value"].astype(str),
            color=links["color"]
        )
    )])

    # Grafik başlığı
    fig.update_layout(
        title_text=f"Sankey Diagramı: {group}_ Grubu için Cluster -> Family",
        font_size=10
    )

    # Diyagramı PDF olarak kaydet
    pdf_file = os.path.join(output_dir, f"sankey_{group}_Regulator_diagram.pdf")
    fig.write_image(pdf_file, format="pdf")
    print(f"{group}_ grubu için Sankey diyagramı {pdf_file} olarak kaydedildi.")
