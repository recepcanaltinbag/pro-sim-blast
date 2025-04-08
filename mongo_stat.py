import pymongo
import re
from collections import Counter
from pymongo.mongo_client import MongoClient
from dotenv import load_dotenv
import os
from wordcloud import WordCloud
import matplotlib.pyplot as plt

import pandas as pd 



load_dotenv()
MONGODB_URI = os.environ['MONGODB_URI']

client = MongoClient(MONGODB_URI)
# Send a ping to confirm a successful connection
#try:
#    client.admin.command('ping')
#    print("Pinged your deployment. You successfully connected to MongoDB!")
#except Exception as e:
#    print(e)

ROAR_DB = client['ROAR-DB']
RelativeGenes_Collection = ROAR_DB['RelativeGenes']
Rieske_Collection = ROAR_DB['Rieske']

total_documents = RelativeGenes_Collection.count_documents({})
processed = 0

# RelativeGenes koleksiyonundan verileri sorgulayın
query = {}  # İhtiyaç duyduğunuz filtreleme işlemleri burada yapılabilir
cursor = RelativeGenes_Collection.find(query)

# 'regulator' kelimesini içeren product keylerini toplayalım
regulator_list = []
transposon_list = []
transposon_having_Rieske = []
regulator_having_Rieske = []
for_each_count = 0
not_empty_ones = 0
transposon_count = 0

for entry in cursor:
    is_regulator = False
    is_transposon = False
    not_empty = False
    flank = entry.get('flank', [])  # flank, bir dizi olduğu için varsayılan olarak boş liste alıyoruz

    # flank içindeki her bir objeyi kontrol edelim
    for item in flank:
        product = item.get('product', '').lower()  # Küçük harfe çevirerek case insensitive yapıyoruz
        not_empty = True

        # Eğer 'regulator' kelimesi product içerisinde geçiyorsa
        if 'regulator' in product:
            regulator_list.append(product)
            regulator_having_Rieske.append((entry['_id'], product))

            is_regulator = True
        
        if 'insertion sequence' in product or 'transposase' in product or 'transposon' in product:
            transposon_list.append(product)
            transposon_having_Rieske.append((entry['_id'], product))

            is_transposon = True

    if is_regulator:
        for_each_count = for_each_count + 1
    if not_empty:
        not_empty_ones = not_empty_ones + 1
    if is_transposon:
        transposon_count = transposon_count + 1
print('finished')
ids = [element[0] for element in transposon_having_Rieske]
# Tek bir sorguda tüm ilgili belgeleri al
documents = Rieske_Collection.find({"_id": {"$in": ids}}, {"_id": 1, "PredictedCluster": 1})
# Gelen belgelerden PredictedCluster değerlerini ekle
tr_rieske_list = [doc["PredictedCluster"] for doc in documents]
documents = Rieske_Collection.find({"_id": {"$in": ids}}, {"_id": 1, "PredictedCluster": 1})
print(documents[0])
id_to_predicted_cluster = {doc["_id"]: doc["PredictedCluster"] for doc in documents}

print(id_to_predicted_cluster)

# PredictedCluster ve element[1]'i birleştir
cluster_product_list = [
    (element[0], id_to_predicted_cluster[element[0]], element[1]) for element in transposon_having_Rieske
]

idsR = [element[0] for element in regulator_having_Rieske]
documentsR = Rieske_Collection.find({"_id": {"$in": idsR}}, {"_id": 1, "PredictedCluster": 1})
id_to_predicted_clusterR = {doc["_id"]: doc["PredictedCluster"] for doc in documentsR}

regulator_product_list = [
    (element[0], id_to_predicted_clusterR[element[0]], element[1]) for element in regulator_having_Rieske
]

df = pd.DataFrame(cluster_product_list, columns=['id','Cluster', 'Product'])
df.to_csv('trees/transposon_having_Rieske.csv', index=False, encoding='utf-8')
print('file writed')


df = pd.DataFrame(regulator_product_list, columns=['id','Cluster', 'Product'])
df.to_csv('trees/regulator_having_Rieske.csv', index=False, encoding='utf-8')
print('file writed')
input()
# Regulator terimlerinin istatistiğini çıkaralım
regulator_counts = Counter(regulator_list)

# İstatistiği sıralayalım (en çok geçen terim en üstte olacak şekilde)
sorted_regulator_counts = regulator_counts.most_common()

# Sonuçları yazdıralım
for term, count in sorted_regulator_counts:
    print(f"{term}: {count}")

print('-----TRANPSONS')
transposon_counts = Counter(transposon_list)

# İstatistiği sıralayalım (en çok geçen terim en üstte olacak şekilde)
sorted_transposon_counts = transposon_counts.most_common()

# Sonuçları yazdıralım
for term, count in sorted_transposon_counts:
    print(f"{term}: {count}")


print('-----Rieske')
transposonRieske_counts = Counter(tr_rieske_list)

# İstatistiği sıralayalım (en çok geçen terim en üstte olacak şekilde)
sorted_transposonRieske_counts = transposonRieske_counts.most_common()

# Sonuçları yazdıralım
for term, count in sorted_transposonRieske_counts:
    print(f"{term}: {count}")


print('len:', len(regulator_list))
print('len for each: ', for_each_count)
print('not empty ones:', not_empty_ones)
print('transposon:', len(transposon_list))
print('len for each: ', transposon_count)

# Her kelimeyi ayırıp tekrar sayma
words = ' '.join(regulator_list).split()
word_counts = Counter(words)

# Word Cloud oluşturma
wordcloud = WordCloud(width=2000, height=1000, background_color='white').generate_from_frequencies(word_counts)

'''
# Görselleştirme
plt.figure(figsize=(10, 5))
plt.imshow(wordcloud, interpolation='bilinear')
plt.axis('off')  # Eksenleri kaldır

# PDF olarak kaydetme
plt.savefig("trees/wordcloud_individual.pdf", format='pdf', dpi=300)

# Görselleştirmeyi gösterme (isteğe bağlı)
plt.show()
'''


# Her kelimeyi ayırıp tekrar sayma
words = ' '.join(transposon_list).split()
word_counts = Counter(words)

# Word Cloud oluşturma
wordcloud = WordCloud(width=2000, height=1000, background_color='white').generate_from_frequencies(word_counts)


# Görselleştirme
plt.figure(figsize=(10, 5))
plt.imshow(wordcloud, interpolation='bilinear')
plt.axis('off')  # Eksenleri kaldır

# PDF olarak kaydetme
plt.savefig("trees/transposon_individual.pdf", format='pdf', dpi=300)

# Görselleştirmeyi gösterme (isteğe bağlı)
plt.show()


'''
tetr/acrr family transcriptional regulator: 3294
lysr family transcriptional regulator: 2734
iclr family transcriptional regulator: 1830
marr family winged helix-turn-helix transcriptional regulator: 1479
gntr family transcriptional regulator: 1175
arac family transcriptional regulator: 899
response regulator transcription factor: 731
helix-turn-helix transcriptional regulator: 642
metalloregulator arsr/smtb family transcription factor: 601
laci family dna-binding transcriptional regulator: 574
antar domain-containing response regulator: 425
glxa family transcriptional regulator: 407
lrp/asnc family transcriptional regulator: 343
response regulator: 331
iclr family transcriptional regulator c-terminal domain-containing protein: 247
merr family transcriptional regulator: 236
dna-binding transcriptional regulator hcar: 227



tetr/acrr, 3294
lysr, 2734
iclr, 1830
marr, 1479
gntr, 1175
arac, 899
response, 731
helix-turn-helix, 642
metalloregulator arsr/smtb, 601
laci, 574
antar, 425
glxa, 407
lrp/asnc, 343
response 
iclr family 
merr family 
dna-binding 
'''