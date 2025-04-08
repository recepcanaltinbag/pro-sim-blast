import pandas as pd
import plotly.graph_objects as go
import os
# CSV dosyasını oku
csv_file = "trees/regulator_having_Rieske.csv"
df = pd.read_csv(csv_file, delimiter=",")

# Grup isimlerini düzenle
group_mapping = {
    "1_": "Group I",
    "2_": "Group II",
    "3_": "Group III",
    "4_": "Group IV",
    "5_": "Group V",
}
df["Group"] = df["Cluster"].str.extract(r"^(\d+_)").replace(group_mapping)

# Family sütununu oluştur
df["Family"] = df["Product"].str.split(r"[-/ ]").str[0]

# Sankey için veri hazırlama
# Group -> Family
group_family = df.groupby(["Group", "Family"]).size().reset_index(name="Value")

# Total değeri hesapla ve oranlara göre filtre uygula
total_value = group_family["Value"].sum()
threshold = 0.025  # %5 eşik değeri
group_family["Family"] = group_family.apply(
    lambda row: row["Family"] if row["Value"] / total_value >= threshold else "Other",
    axis=1
)

# "Other" için yeniden gruplama
group_family = group_family.groupby(["Group", "Family"]).sum().reset_index()

# Düğümleri oluştur
groups = group_family["Group"].unique()
families = group_family["Family"].unique()

all_nodes = list(groups) + list(families)

# İndeks haritalama
group_indices = {group: i for i, group in enumerate(groups)}
family_indices = {family: i + len(groups) for i, family in enumerate(families)}

# Linkleri oluştur
links = pd.DataFrame({
    "source": group_family["Group"].map(group_indices),
    "target": group_family["Family"].map(family_indices),
    "value": group_family["Value"]
})

# Renkleri ayarla
color_palette = [
    "rgba(31, 119, 180, 0.4)",  # Group I
    "rgba(255, 127, 14, 0.4)",  # Group II
    "rgba(44, 160, 44, 0.4)",   # Group III
    "rgba(214, 39, 40, 0.4)",   # Group IV
    "rgba(148, 103, 189, 0.4)"  # Group V
]
node_colors = (
    [color_palette[i % len(color_palette)] for i in range(len(groups))] +  # Groups için renkler
    ["rgba(200, 200, 200, 0.4)" if family == "Other" else "rgba(100, 100, 255, 0.4)" for family in families]  # Families için renkler
)

# Link renklerini grupların renklerine göre ayarla
link_colors = links["source"].map(lambda x: color_palette[x])

# Sankey diyagramı oluştur
fig = go.Figure(data=[go.Sankey(
    node=dict(
        pad=15,
        thickness=15,
        line=dict(color="black", width=0.5),
        label=all_nodes,
        color=node_colors
    ),
    link=dict(
        source=links["source"],
        target=links["target"],
        value=links["value"],
        color=link_colors  # Link renkleri
    )
)])

# Grafik başlığı
fig.update_layout(
    title_text="Group -> Family Sankey Diagram (Threshold Applied)",
    font_size=10
)

# Diyagramı göster
output_dir = "sankey_outputs"

fig.show()
pdf_file = os.path.join(output_dir, f"sankey_combined_R_diagram.pdf")
fig.write_image(pdf_file, format="pdf")