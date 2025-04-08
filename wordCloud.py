
from wordcloud import WordCloud
import matplotlib.pyplot as plt
from collections import Counter

# Verilen kelimeler ve tekrar sayıları
data = """
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
"""

# Veriyi işleyelim
word_counts = {}
for line in data.split('\n'):
    if line.strip():
        word, count = line.split(',')
        word_counts[word.strip()] = int(count.strip())

# Word Cloud oluşturma
wordcloud = WordCloud(width=800, height=400, background_color='white').generate_from_frequencies(word_counts)

# Görselleştirme
plt.figure(figsize=(10, 5))
plt.imshow(wordcloud, interpolation='bilinear')
plt.axis('off')  # Eksenleri kaldır
plt.show()
