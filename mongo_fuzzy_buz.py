from rapidfuzz import process
import pandas as pd
from pymongo.mongo_client import MongoClient
from dotenv import load_dotenv
import os


# Örnek veri
def data_fuzz_analysis(data):
    '''
    data = [
        {"Product": "MFS transporter", "Start": 255537, "End": 256743, "Relative": 7203},
        {"Product": "assimilatory nitrate reductase electron transfer subunit NasB", "Start": 256922, "End": 259238, "Relative": 4708},
        {"Product": "assimilatory nitrate reductase catalytic subunit NasC", "Start": 259244, "End": 261377, "Relative": 2569},
        {"Product": "NADPH-nitrite reductase", "Start": 261497, "End": 263915, "Relative": 31},
        {"Product": "MFS type transporter", "Start": 264000, "End": 265000, "Relative": 1500},
        {"Product": "nitrate reductase NasD", "Start": 265100, "End": 267000, "Relative": 3400},
    ]
    '''

    # Veriyi DataFrame'e dönüştürme
    df = pd.DataFrame(data)
    sample_size = len(data)
    # Benzerlik eşik değeri
    threshold = 90

    '''
    # Grupları oluşturmak için bir liste
    groups = []
    group_map = {}  # Ürün ve grup eşlemesi

    # Gruplama işlemi
    for product in df['product']:
        # Mevcut gruplarla benzerlik kontrolü
        match = process.extractOne(product, groups, score_cutoff=threshold)
        if match:
            # Benzer grup bulunduysa aynı gruba ekle
            group_map[product] = match[0]
        else:
            # Yeni bir grup oluştur
            groups.append(product)
            group_map[product] = product

    # Grup bilgilerini DataFrame'e ekle
    df['Group'] = df['product'].map(group_map)
    '''

    # Ayrı gruplarda olması gereken anahtar kelimeler ve grupları
    exclusive_keywords = {
        "transcriptional regulator": "G_Regulator",
        "transporter": "G_Transporter",
        "transposase": "G_Transposase",
        "regulator": "G_Regulator",
        "insertion sequence": "G_Transposase",
        "reductase": "G_Reductase",
        "ferredoxin": "G_Ferredoxin",
        "oxygenase small subunit": "G_OxygenaseSubunit",
        "oxygenase subunit beta": "G_OxygenaseSubunit"
    }

    # Grup bilgilerini takip edecek yapılar
    groups = []
    group_map = {}

    for product in df['product']:
        # Önce exclusive_keywords içinde kontrol et
        assigned = False
        for keyword, group_name in exclusive_keywords.items():
            if keyword in product.lower():  # Büyük/küçük harf duyarlılığını kaldırmak için
                group_map[product] = group_name
                assigned = True
                break
        
        if not assigned:
            # Mevcut gruplarla benzerlik kontrolü
            match = process.extractOne(product, groups, score_cutoff=threshold)
            if match:
                # Benzer grup bulunduysa aynı gruba ekle
                group_map[product] = match[0]
            else:
                # Yeni bir grup oluştur
                groups.append(product)
                group_map[product] = product

    # Grup bilgilerini DataFrame'e ekle
    df['Group'] = df['product'].map(group_map)


    # Her grup için puan hesaplama

    '''
    grouped = df.groupby('Group').apply(
        lambda g: pd.Series({
            "Score": (1 / g['relative'].apply(lambda r: ( abs(r)) / sample_size).sum() ) * 1000,
            "Products": list(g['product'])  # Gruplara ait ürünler
        })
    ).sort_values(by='Score', ascending=False)
    '''
    EPSILON = 1e-6  # Sıfıra bölünmeyi önlemek için küçük bir sabit

    grouped = df.groupby('Group').apply(
    lambda g: pd.Series({
    #"Score": (1 / ((g['relative'].sum()+EPSILON)/len(g))) * ((len(g)**2)/sample_size) * 1000,
    "Score": (1/(abs(g['relative'])+1)).sum() * (len(g)**2/sample_size) * 100,
    "Average R": g['relative'].sum()/len(g),
    "Total R": g['relative'].sum(),
    "Total Gene": len(g),
    "Products": list(g['product']),  # Gruplara ait ürünler
    })
    ).sort_values(by='Score', ascending=False)

    # Sonuçları inceleme
    print('Groups were determined')
    #print("Gruplar ve Puanları:")
    #print(grouped)
    return grouped


load_dotenv()
MONGODB_URI = os.environ['MONGODB_URI']
client = MongoClient(MONGODB_URI)
ROAR_DB = client['ROAR-DB']
Rieske_Collection = ROAR_DB['Rieske']
Relative_Collection = ROAR_DB['RelativeGenes']



distinct_values = Rieske_Collection.distinct("PredictedCluster")


with open(f"out_Relative_Infov2/GeneralInfo.txt", "w") as file1:

    file1.write(f"-------------------------------------\n")

print('File writed')

the_value = "ALL"
for the_value in distinct_values:
    documents = Rieske_Collection.find({"PredictedCluster": the_value}, {"_id": 1})



    #documents = Rieske_Collection.find()
    # _id'leri listeye almak
    ids = [doc["_id"] for doc in documents]


    flank_array = []
    seq_type_array = []

    # _id'ye göre dökümanları sorgulayıp işlemleri yapma

    documents_Relative = Relative_Collection.find({"_id": {"$in": ids}}, {"flank": 1, "seq_type": 1})

    # Sorgudan dönen her bir dökümanı işleme
    for document in documents_Relative:
        if "flank" in document and isinstance(document["flank"], list):
            flank_array.extend(document["flank"])  # Array'i mevcut array'e ekle
        if "seq_type" in document:
            seq_type_array.append(document["seq_type"])


    total_count = len(seq_type_array)

    # "genome" ve "plasmid" sayısını say
    genome_count = seq_type_array.count("Genome")
    plasmid_count = seq_type_array.count("Plasmid")
    if total_count == 0:
        genome_percentage = 0
        plasmid_percentage = 0
    else:
    # Yüzdelik hesaplamalar
        genome_percentage = (genome_count / total_count) * 100
        plasmid_percentage = (plasmid_count / total_count) * 100

    # Sonuçları yazdırma
    print(the_value)
    print(len(ids))
    print(f"Genome Yüzdesi: {genome_percentage:.2f}%")
    print(f"Plasmid Yüzdesi: {plasmid_percentage:.2f}%")
    print(len(flank_array))
    if len(flank_array) == 0:
        print('no array info, skipping')
        continue
    groups_dataframe = data_fuzz_analysis(flank_array)

        # Dosyayı açma (yazma modunda)

    groups_dataframe.to_csv(f"out_Relative_Infov2/{the_value}.csv", float_format="%.2f")  # `index=False` ile index'i yazmadan kaydediyoruz.

    with open(f"out_Relative_Infov2/GeneralInfo.txt", "a") as file1:
        # Verilerinizi hesapladıktan sonra
        file1.write(f"The Value: {the_value}\n")
        file1.write(f"Total IDs: {len(ids)}\n")
        file1.write(f"Relative List: {len(flank_array)}\n")
        file1.write(f"Genome Yüzdesi: {genome_percentage:.2f}%\n")
        file1.write(f"Plasmid Yüzdesi: {plasmid_percentage:.2f}%\n")
        file1.write(f"-------------------------------------\n")

    print('File writed')