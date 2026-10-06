"""
Veriden ogrenme -- "iyi bir veri setimiz var, daha iyi korelasyonlar bulabiliriz"
sorusunu DENETLENEBILIR bir denetimli ogrenme problemine cevirir.

NEDEN BU MODUL VAR
  Veritabaninda 11.422 dogrulanmis alpha alt birimi var ve her biri bir tipe,
  her tip de kuratorlu bir reaksiyon sinifina ve substrat sinifina bagli. Bu
  haliyle akla gelen ilk sey "dizi -> reaksiyon sinifi" siniflandiricisi
  kurmaktir. Ama bu veri setinde o is, dogru yapilmazsa, sonucu BASTAN BELLI
  bir sayi uretir.

SIZINTI (LEAKAGE) PROBLEMI -- bu modulun asil konusu
  Girisler BAGIMSIZ DEGIL. Iki ayri sebepten:
    1. Fazlalik: 11.422 giris yalnizca 10.019 tekil dizi (bkz. redundancy.json).
       Ayni protein birbirine cok yakin suslardan tekrar tekrar geliyor.
    2. Etiket tipin OZELLIGI: reaksiyon sinifi ve substrat sinifi tek tek
       girislere degil `ro.ro_cluster` tipine atanmis. Yani bir tipin bir uyesi
       egitimde baska bir uyesi testte oldugunda model dizi-reaksiyon iliskisini
       DEGIL, "bu dizi hangi tipe benziyor" sorusunu cozuyor ve cevabi egitim
       kumesinden okuyor.
  Sonuc: rastgele bolunmus bir train/test ayrimi yuksek ve ANLAMSIZ bir dogruluk
  verir. Bu sayi olculdu ve raporlaniyor -- ama SIZINTILI diye etiketlenerek,
  cunku onu gizlemek yerine gostermek ogreticidir.

NASIL ONLENDI
  * Butun basliklar GRUPLANMIS capraz dogrulama ile olculuyor; grup = enzim TIPI
    (`ro.ro_cluster`). Ayni tip hem egitimde hem testte asla bulunmuyor. Iki
    sema kosuluyor: GroupKFold (tip bazli, giris sayisina gore dengeli 10 kat)
    ve her tipi sirayla disarida birakma (leave-one-type-out, 61 kat).
  * Ayni model ve ayni hiperparametreler rastgele bolunmus semayla da kosuluyor;
    iki sayi arasindaki UCURUM bir sonuc olarak veriliyor. Ucurum model
    KAPASITESIYLE buyuyor; bu, ogrenmenin degil ezberlemenin imzasidir.
  * Dort taban cizgisi: (a) global cogunluk sinifi orani, (b) tip duzeyi
    cogunluk orani, (c) kat basina egitim cogunlugu, (d) KLAD KAHINI -- test
    girisinin HMM grubunu (`ro.ro_group`, dizi benzerliginden turer, yani yeni
    bir protein icin de bilinir) kullanip o grubun egitim tiplerindeki cogunluk
    etiketini soyleyen kural. (d) en onemlisi: hizalama kolonlarinin derin klad
    uyeliginin UZERINE bir sey katip katmadigini olcer.
  * Etiket permutasyon null'i: etiketler TIP duzeyinde karistirilip ayni
    gruplanmis sema yeniden kosulur, yuzlerce tekrar, ve null dagiliminin
    ortalamasi, standart sapmasi, %95'lik dilimi ve azamisi yazilir. Permutasyon
    tip duzeyinde yapilir, cunku etiket tipin ozelligi; giris duzeyinde
    karistirmak null'i yapay olarak dusururdu. Tipin ozelligi OLMAYAN tek
    baslikta (beta alt birimi) karistirma giris duzeyindedir.
  * Model ACIK: derinligi sinirli karar agaci (bolme kurallari JSON'a yaziliyor)
    ve L1 cezali cok sinifli lojistik regresyon. Kapali bir topluluk (rastgele
    orman) YALNIZCA degisken onemi icin kuruluyor ve ciktida boyle isaretli.
  * Agac derinligi test dogrulugundan secilmedi: IC ICE (nested) gruplanmis CV
    ile her dis katin kendi egitim verisi uzerinde secildi, secilen derinlikler
    de JSON'a yazildi. Butun derinlik taramasi ayrica veriliyor ama "en iyi
    derinlik" o taramadan okunmadi.
  * Her girisin agirligi 1/(tipinin buyuklugu): 1.345 uyeli bir tip ile 2 uyeli
    bir tip egitimde esit soz sahibi olsun diye. Agirliksiz surum duyarlilik
    olarak kosuluyor.
  * Dizi ozellikleri icin kolon secimi UC sekilde yapiliyor: (1) hazir
    `analysis_out/sdp_positions.csv` ilk N kolonu -- bu dosya kolonlari TUM
    tipleri gorerek sectigi icin bir SECIM SIZINTISI tasir ve boyle etiketlenir;
    (2) ayni olcut her katin YALNIZCA egitim tipleriyle yeniden hesaplanarak
    (ic-kat secim, sizintisiz); (3) hic secim yapmadan dolu kolonlarin tamami.
    Ucunun gruplanmis dogrulugu yan yana veriliyor, yani hazir listeyi
    kullanmanin bedeli sayiyla gosteriliyor.

NE OLCULDU
  1. Reaksiyon sinifi (7 sinif, chemistry.csv) dizinin kendisinden -- hizalama
     kolonlarindaki kalinti kimlikleri (genomic_context/cand_aln.sto match
     kolonlari, model uzunlugu 426).
  2. Hangi kolonlar sinyali tasiyor: agac bolmeleri, L1 katsayilari, orman
     onemleri ve TIP duzeyinde tek degiskenli Cramer's V. Her kolon
     `analysis_out/motif_stats.json`'daki sekiz tanimlayici kolona gore isaretli
     uzaklikla ve bolge adiyla (rieske / baglanti / katalitik) raporlanir,
     boylece "aktif bolgede mi" sorusu biyolog tarafindan sorulabilir.
  3. Substrat sinifi (3 sinif, cluster_ecology.csv) GENOMIK BAGLAMDAN: operon
     bilesimi, elektron transfer ortak tipleri, duzenleyici ailesi, plazmit
     durumu ve taksonomi. Taksonomi AYRI blok tutulur ve ablasyonla kosulur,
     cunku taksonomi ile baglam ic ice gecmistir: "komsuluk kimyayi ongoruyor"
     iddiasi ancak taksonomi cikarildiginda sinanabilir.
  4. Ek olarak: (a) POZITIF KONTROL -- `ro.ro_group` (HMM grubu) diziden
     ongorulur; bu etiket dizi benzerliginden turedigi icin yontemin sinyali
     YAKALAYABILDIGINI gosterir, yani dusuk dogruluklar yontemin korlugu degil;
     (b) modalite karsilastirmasi -- ayni etiket ayni katlarla bir kez diziden
     bir kez baglamdan ongorulur, boylece bilginin NEREDE oldugu sorulur;
     (c) tipin ozelligi OLMAYAN bir etiket: beta alt biriminin varligi (alpha3 -
     alpha3beta3 mimarisi), tip icinde degistigi icin girislerin kendi bilgisi
     tasidigi tek baslik.

NE OLCULMUYOR / VERININ DESTEKLEMEDIGI SONUCLAR
  * Etkin ornek buyuklugu 11.422 degil 61'dir (tip sayisi). Yedi sinifli bir
    problemde 61 bagimsiz birim cok azdir; guven araliklari genistir ve kucuk
    siniflar (C_N_cleavage 2 tip, O_demethylation 3 tip, angular_dioxygenation
    3 tip) icin hicbir genelleme iddiasi kurulamaz.
  * Ayirt edici bulunan kolonlar YAPISAL olarak dogrulanmis substrat baglama
    cebi kalintilari DEGILDIR. Veritabaninda yapi yok; olculen sey "siniflari
    ayiran hizalama kolonlari". Motif kolonlarina uzaklik bir HIPOTEZ uretir,
    kanit uretmez.
  * Reaksiyon sinifi da substrat sinifi da HMM grubuyla (derin filogenetik
    bolunme) guclu bicimde yapilanmis. Bu yuzden bir dizi modeli kladi tanimakla
    da ayni dogruluga ulasabilir; klad kahini taban cizgisi tam bunu gorunur
    kilmak icin var.
  * Substrat ve reaksiyon etiketleri tip duzeyinde kuratorludur; kanit kademesi
    `distant`/`novel` olan uyelerin kimyasi aslinda bilinmiyor. Bu modul etiketi
    dogru varsayar; dogrulugunu sinamaz.
  * Baglamdan substrat sinifi ongorusu bir NEDENSELLIK iddiasi degildir. Ayni
    genomda bulunmak ortak yol kaniti degil; ustelik dizilenmis genomlar orneklem
    yanlilidir (PAH yikan izolatlar fazla temsil edilir).

Cikti:
    analysis_out/learning.json
"""

import argparse
import csv
import json
import math
import os
import sqlite3
import sys
import time
from collections import Counter, defaultdict

import numpy as np

from ro_motif import DEFAULT_MODEL, MODEL_COLUMNS, read_stockholm_matchcols

# --- Sabitler -------------------------------------------------------------

# Grup = enzim tipi. Bu modulde her yerde boyle; degistirilecek bir ayar degil,
# metodolojinin kendisi.
GROUP_FIELD = "ro.ro_cluster"

N_GROUP_FOLDS = 10        # tip bazli GroupKFold kat sayisi
N_RANDOM_FOLDS = 10       # sizintili rastgele sema ayni kat sayisini kullanir
MIN_RESIDUE_COUNT = 25    # bir (kolon, kalinti) ikilisi ozellik olmak icin
TREE_MIN_LEAF = 5
DEPTH_SWEEP = (1, 2, 3, 4, 5, 6, 8, 12, 20)
NESTED_DEPTHS = (2, 3, 4, 5, 6, 8)
PRIMARY_DEPTH = 5         # tek-derinlikli sayilar icin; ic ice CV ayrica kosulur

# variant_signature.py ile AYNI konvansiyon; bolge sinirlari yaklasiktir.
RIESKE_END = 120
CATALYTIC_START = 191

# Rastgelelik tek yerden; her kosu ayni sayilari uretsin.
SEED = 20261006


# --- Kucuk yardimcilar ----------------------------------------------------

def entropy(counter):
    total = sum(counter.values())
    if total <= 1:
        return 0.0
    return -sum((n / total) * math.log2(n / total) for n in counter.values() if n)


def zone_of(column):
    """Kolonun bolgesi -- variant_signature.py'deki ayni konvansiyon."""
    if column <= RIESKE_END:
        return "rieske_domain"
    if column < CATALYTIC_START:
        return "linker"
    return "catalytic_domain"


def cramers_v(table):
    """Chi-kare tabanli Cramer's V; bos satir ve sutun atilir."""
    table = [row for row in table if sum(row) > 0]
    if len(table) < 2:
        return 0.0
    cols = len(table[0])
    col_sums = [sum(row[j] for row in table) for j in range(cols)]
    keep = [j for j in range(cols) if col_sums[j] > 0]
    if len(keep) < 2:
        return 0.0
    table = [[row[j] for j in keep] for row in table]
    rows_n, cols_n = len(table), len(table[0])
    n = sum(sum(r) for r in table)
    if n == 0:
        return 0.0
    row_sums = [sum(r) for r in table]
    col_sums = [sum(table[i][j] for i in range(rows_n)) for j in range(cols_n)]
    chi2 = 0.0
    for i in range(rows_n):
        for j in range(cols_n):
            expected = row_sums[i] * col_sums[j] / n
            if expected > 0:
                chi2 += (table[i][j] - expected) ** 2 / expected
    return float((chi2 / (n * min(rows_n - 1, cols_n - 1))) ** 0.5)


def codes(values):
    """Dizgi listesini (tamsayi kodlar, sirali etiketler) ikilisine cevir."""
    order = sorted(set(values))
    index = {value: i for i, value in enumerate(order)}
    return np.array([index[v] for v in values], dtype=np.int32), order


# --- Karar agaci (numpy) --------------------------------------------------
#
# NEDEN KENDI UYGULAMASI. Projenin python ortaminda scikit-learn KURULU DEGIL
# (`python3 -c "import sklearn"` basarisiz; run_all.sh yalnizca biopython,
# pandas ve numpy sart kosuyor). Bagimlilik eklemek yerine ihtiyac duyulan iki
# model numpy ile yazildi; boylece modul pipeline'in geri kalaniyla ayni ortamda
# kosuyor ve her makinede ayni sayilari uretiyor. Ozellikler ikili (0/1) oldugu
# icin CART cok basitlesiyor: her dugumde tek bir BLAS carpimi butun
# ozelliklerin sinif sayimlarini birden veriyor.

def fit_tree(x, y_onehot, weight, max_depth, min_leaf=TREE_MIN_LEAF,
             feature_mask=None):
    """Gini olcutlu, ikili ozellikli CART. Donen: ic ice sozluk."""
    n_classes = y_onehot.shape[1]

    def build(index, depth):
        y_sub = y_onehot[index] * weight[index, None]
        totals = y_sub.sum(axis=0)
        node = {"n": int(index.size), "class_weight": totals,
                "prediction": int(np.argmax(totals))}
        if depth >= max_depth or index.size < 2 * min_leaf:
            return node
        if (totals > 0).sum() <= 1:
            return node
        x_sub = x[index]
        ones = np.ones((index.size, 1), dtype=np.float32)
        counts = np.concatenate([y_sub, ones], axis=1).T @ x_sub
        class_one = counts[:n_classes]            # agirlikli, ozellik == 1
        plain_one = counts[n_classes]             # agirliksiz adet, ozellik == 1
        weight_one = class_one.sum(axis=0)
        weight_total = totals.sum()
        weight_zero = weight_total - weight_one
        class_zero = totals[:, None] - class_one
        safe_one = np.where(weight_one > 0, weight_one, 1.0)
        safe_zero = np.where(weight_zero > 0, weight_zero, 1.0)
        gini_one = 1.0 - ((class_one / safe_one) ** 2).sum(axis=0)
        gini_zero = 1.0 - ((class_zero / safe_zero) ** 2).sum(axis=0)
        impurity = (weight_one * gini_one + weight_zero * gini_zero) / weight_total
        # Bir yana min_leaf'ten az GIRIS dusuyorsa o bolme yasak. Olcut agirlik
        # degil adet, cunku agirliklar 1/tip_boyutu ile cok kucuk olabiliyor ve
        # bir agirlik esigi buyuk tipleri kayiriyordu.
        too_small = (plain_one < min_leaf) | ((index.size - plain_one) < min_leaf)
        impurity = np.where(too_small, np.inf, impurity)
        if feature_mask is not None:
            impurity = np.where(feature_mask, impurity, np.inf)
        best = int(np.argmin(impurity))
        if not np.isfinite(impurity[best]):
            return node
        parent_gini = 1.0 - ((totals / weight_total) ** 2).sum()
        gain = parent_gini - impurity[best]
        if gain <= 1e-12:
            return node
        goes_right = x_sub[:, best] > 0.5
        node["feature"] = best
        node["gain"] = float(gain * weight_total)
        node["right"] = build(index[goes_right], depth + 1)
        node["left"] = build(index[~goes_right], depth + 1)
        return node

    return build(np.arange(x.shape[0]), 0)


def predict_tree(node, x):
    out = np.zeros(x.shape[0], dtype=np.int32)

    def walk(node, index):
        if "feature" not in node or index.size == 0:
            out[index] = node["prediction"]
            return
        goes_right = x[index, node["feature"]] > 0.5
        walk(node["right"], index[goes_right])
        walk(node["left"], index[~goes_right])

    walk(node, np.arange(x.shape[0]))
    return out


def tree_gains(node, n_features):
    """Bolme kazanclarini ozellik basina topla -- agacin kendi onem olcutu."""
    gains = np.zeros(n_features)

    def walk(node):
        if "feature" in node:
            gains[node["feature"]] += node["gain"]
            walk(node["left"])
            walk(node["right"])

    walk(node)
    return gains


def tree_rules(node, names, class_names, max_depth=3, sequence=True):
    """Agaci okunabilir satirlara cevir -- JSON'da model boyle denetlenir."""
    lines = []

    def label_of(index):
        if sequence:
            column, residue = names[index]
            return f"col{column}={residue}"
        return str(names[index])

    def walk(node, prefix, level):
        total = float(node["class_weight"].sum())
        share = (float(node["class_weight"][node["prediction"]]) / total
                 if total > 0 else 0.0)
        if "feature" not in node or level >= max_depth:
            lines.append({"path": prefix or "root", "n_entries": node["n"],
                          "predicted_class": class_names[node["prediction"]],
                          "class_share_of_weight": round(share, 3)})
            return
        test = label_of(node["feature"])
        joiner = " AND " if prefix else ""
        walk(node["right"], prefix + joiner + test, level + 1)
        walk(node["left"], prefix + joiner + "NOT " + test, level + 1)

    walk(node, "", 0)
    return lines


# --- L1 cezali cok sinifli lojistik regresyon (numpy, FISTA) --------------
#
# Yorumlanabilirlik icin L1: katsayilarin cogu tam sifir kaliyor, yani model
# "hangi kolonlar" sorusuna dogrudan cevap veriyor. Cozucu hizlandirilmis
# proksimal gradyan (FISTA); ozellikler ikili oldugu icin yalnizca egitim
# ortalamasi cikariliyor, ayrica olcekleme gerekmiyor. Kesisim terimi
# cezalandirilmiyor ve sinif onsel olasiliklarinin logaritmasiyla baslatiliyor.

def fit_l1_logistic(x, y_index, weight, n_classes, penalty=0.002,
                    iterations=150):
    n, n_features = x.shape
    total_weight = max(float(weight.sum()), 1e-12)
    mean = (x * weight[:, None]).sum(axis=0) / total_weight
    centered = (x - mean).astype(np.float32)
    target = np.zeros((n, n_classes), dtype=np.float32)
    target[np.arange(n), y_index] = 1.0
    norm = (weight / total_weight).astype(np.float32)

    # Lipschitz sabiti: guc iterasyonuyla en buyuk ozdegerin yarisi.
    vector = np.random.RandomState(0).randn(n_features).astype(np.float32)
    vector /= max(float(np.linalg.norm(vector)), 1e-12)
    eigenvalue = 1.0
    for _ in range(12):
        vector = centered.T @ ((centered @ vector) * norm)
        eigenvalue = float(np.linalg.norm(vector))
        if eigenvalue < 1e-20:
            break
        vector /= eigenvalue
    lipschitz = max(eigenvalue * 0.5, 1e-6)
    step = 1.0 / lipschitz

    priors = np.array([float((norm * (y_index == k)).sum())
                       for k in range(n_classes)], dtype=np.float32)
    bias = np.log(np.maximum(priors, 1e-6))
    coef = np.zeros((n_features, n_classes), dtype=np.float32)
    momentum = coef.copy()
    previous = coef.copy()
    theta = 1.0
    for _ in range(iterations):
        scores = centered @ momentum + bias
        scores -= scores.max(axis=1, keepdims=True)
        np.exp(scores, out=scores)
        scores /= scores.sum(axis=1, keepdims=True)
        residual = (scores - target) * norm[:, None]
        gradient = centered.T @ residual
        # Kesisim cezasiz: Hessian siniri 0,5 oldugu icin adim 1/0,5 = 2.
        bias = bias - 2.0 * residual.sum(axis=0)
        candidate = momentum - step * gradient
        coef = np.sign(candidate) * np.maximum(np.abs(candidate) - step * penalty,
                                               0.0)
        theta_next = (1.0 + math.sqrt(1.0 + 4.0 * theta * theta)) / 2.0
        momentum = coef + ((theta - 1.0) / theta_next) * (coef - previous)
        previous = coef.copy()
        theta = theta_next
    return {"coef": coef, "bias": bias, "mean": mean}


def predict_l1_logistic(model, x):
    scores = (x - model["mean"]) @ model["coef"] + model["bias"]
    return np.argmax(scores, axis=1).astype(np.int32)


# --- Rastgele orman (yalnizca degisken onemi icin) ------------------------

def forest_importances(x, y_onehot, weight, n_trees, max_depth, rng,
                       feature_fraction=0.4):
    n, n_features = x.shape
    total = np.zeros(n_features)
    keep = max(1, int(round(n_features * feature_fraction)))
    for _ in range(n_trees):
        rows = rng.randint(0, n, size=n)
        mask = np.zeros(n_features, dtype=bool)
        mask[rng.choice(n_features, size=keep, replace=False)] = True
        node = fit_tree(x[rows], y_onehot[rows], weight[rows], max_depth,
                        feature_mask=mask)
        total += tree_gains(node, n_features)
    return total / max(n_trees, 1)


# --- Capraz dogrulama semalari -------------------------------------------

def group_kfold(group_codes, n_folds):
    """Tip bazli kat ayrimi; katlar GIRIS sayisina gore dengelenir.

    scikit-learn'un GroupKFold'u ile ayni acgozlu kural: tipler buyukten kucuge
    siralanir ve her biri o an en az yuklu kata atanir. Bir tip tek bir kata
    gider, yani ayni tip hem egitimde hem testte ASLA bulunmaz.
    """
    sizes = Counter(group_codes.tolist())
    order = sorted(sizes, key=lambda g: (-sizes[g], g))
    load = [0] * n_folds
    assignment = {}
    for group in order:
        target = int(np.argmin(load))
        assignment[group] = target
        load[target] += sizes[group]
    fold_of = np.array([assignment[g] for g in group_codes.tolist()])
    return [(np.where(fold_of != k)[0], np.where(fold_of == k)[0])
            for k in range(n_folds) if (fold_of == k).any()]


def leave_one_group_out(group_codes):
    return [(np.where(group_codes != g)[0], np.where(group_codes == g)[0])
            for g in range(int(group_codes.max()) + 1)]


def random_kfold(n, n_folds, rng):
    """SIZINTILI sema: giris duzeyinde rastgele bolme, tip yapisi yoksayilir."""
    order = rng.permutation(n)
    return [(np.setdiff1d(order, order[k::n_folds]), order[k::n_folds])
            for k in range(n_folds)]


# --- Olcumler -------------------------------------------------------------

def score(predicted, actual, group_codes, class_names):
    n_classes = len(class_names)
    accuracy = float((predicted == actual).mean())
    confusion = np.zeros((n_classes, n_classes), dtype=int)
    np.add.at(confusion, (actual, predicted), 1)
    per_class, recalls = {}, []
    for k, name in enumerate(class_names):
        support = int(confusion[k].sum())
        predicted_k = int(confusion[:, k].sum())
        hit = int(confusion[k, k])
        recall = hit / support if support else 0.0
        precision = hit / predicted_k if predicted_k else 0.0
        f1 = (2 * precision * recall / (precision + recall)
              if precision + recall > 0 else 0.0)
        per_class[name] = {"support": support, "recall": round(recall, 4),
                           "precision": round(precision, 4), "f1": round(f1, 4)}
        if support:
            recalls.append(recall)
    # Tip duzeyi dogruluk: her tip bir kez sayilir, tahmini uyelerinin
    # cogunlugu. Etkin ornek buyuklugu bu oldugu icin asil okunacak sayi bu.
    hits, total = 0, 0
    for group in range(int(group_codes.max()) + 1):
        members = group_codes == group
        if not members.any():
            continue
        vote = np.bincount(predicted[members], minlength=n_classes).argmax()
        total += 1
        hits += int(vote == actual[members][0])
    return {"accuracy": round(accuracy, 4),
            "balanced_accuracy": round(float(np.mean(recalls)) if recalls else 0.0,
                                       4),
            "type_level_accuracy": round(hits / max(total, 1), 4),
            "n_types": total, "n_entries": int(len(actual)),
            "per_class": per_class,
            "confusion_matrix": confusion.tolist(),
            "confusion_rows_are_true_classes": True}


def nearest_neighbour_accuracy(x, y_index, folds, block=1024):
    """En yakin komsu (1-NN) dogrulugu -- sizintinin en ciplak gosterimi.

    Ozellikler ikili oldugu icin Hamming uzakligi d(i,j) = s_i + s_j - 2*p_ij
    (s satir toplami, p ic carpim). Sabit bir test satiri icin s_i degismedigi
    icin en yakin komsu, 2*p_ij - s_j ifadesini EN BUYUK yapan egitim
    satiridir; boylece tek bir matris carpimi yetiyor.

    Bu siniflandiricinin ogrenecek hicbir seyi yok: "test dizisine en cok
    benzeyen egitim dizisinin etiketini soyle" diyor. Rastgele bolmede yuksek
    cikiyorsa, sebebi veri setinde her dizinin neredeyse ayni bir ikizinin
    bulunmasidir -- yani tam olarak bu modulun engellemeye calistigi sey.
    """
    sums = x.sum(axis=1)
    predicted = np.empty(len(y_index), dtype=np.int32)
    for train, test in folds:
        train_x = np.ascontiguousarray(x[train])
        penalty = sums[train]
        labels = y_index[train]
        for start in range(0, len(test), block):
            chunk = test[start:start + block]
            scores = 2.0 * (x[chunk] @ train_x.T) - penalty
            predicted[chunk] = labels[np.argmax(scores, axis=1)]
    return float((predicted == y_index).mean())


def clade_oracle(folds, group_codes, clade_codes, y_index, n_classes):
    """KLAD KAHINI: test girisinin HMM grubunun egitim cogunlugunu soyle.

    Her klad icinde oy TIP basinadir, giris basina degil: aksi halde tek bir
    buyuk tip kendi kladinin cevabini belirlerdi. `ro_group` dizi benzerliginden
    turedigi icin yeni bir protein icin de bilinebilir, yani bu mesru bir
    dizi-tabanli taban cizgisidir -- ve asil sorusu sudur: hizalama kolonlari
    derin klad uyeliginin UZERINE bir sey katiyor mu.
    """
    predicted = np.zeros(len(y_index), dtype=np.int32)
    for train, test in folds:
        votes = defaultdict(Counter)
        seen = set()
        for i in train:
            group = int(group_codes[i])
            if group in seen:
                continue
            seen.add(group)
            votes[int(clade_codes[i])][int(y_index[i])] += 1
        overall = Counter()
        for counter in votes.values():
            overall.update(counter)
        fallback = overall.most_common(1)[0][0] if overall else 0
        for i in test:
            counter = votes.get(int(clade_codes[i]))
            predicted[i] = counter.most_common(1)[0][0] if counter else fallback
    return predicted


def majority_baselines(folds, y_index, group_codes, class_names):
    n_classes = len(class_names)
    global_majority = int(np.bincount(y_index, minlength=n_classes).argmax())
    global_rate = float((y_index == global_majority).mean())
    predicted = np.empty_like(y_index)
    for train, test in folds:
        predicted[test] = int(np.bincount(y_index[train],
                                          minlength=n_classes).argmax())
    per_fold = score(predicted, y_index, group_codes, class_names)
    type_counts = Counter()
    seen = set()
    for group, label in zip(group_codes.tolist(), y_index.tolist()):
        if group not in seen:
            seen.add(group)
            type_counts[label] += 1
    type_majority = max(type_counts.values()) / sum(type_counts.values())
    return {
        "global_majority_class": class_names[global_majority],
        "global_majority_rate": round(global_rate, 4),
        "type_level_majority_rate": round(type_majority, 4),
        "grouped_train_majority_accuracy": per_fold["accuracy"],
        "note": ("the per-fold train-majority baseline can fall far below the "
                 "global majority rate: holding out a whole type removes a large "
                 "block of one class, so the training majority changes from fold "
                 "to fold. The global rate is the fair reading of 'always answer "
                 "the commonest class'; the type-level rate is its counterpart "
                 "when every type counts once"),
    }


def permutation_null(x, folds, group_codes, y_index, n_classes, weight, depth,
                     reps, rng, level="type", log_every=50):
    """Etiketleri karistir, ayni gruplanmis semayi tekrar kos.

    level="type": tip-etiket eslesmesi karistirilir. Etiket tipin ozelligi
    oldugu icin dogru null budur; giris duzeyinde karistirmak null'i yapay
    olarak dusururdu ve her sonucu "anlamli" gosterirdi.
    level="entry": etiket tipin ozelligi DEGILSE (beta alt birimi) girisler
    arasinda karistirilir.
    """
    n_groups = int(group_codes.max()) + 1
    group_label = np.array([int(y_index[group_codes == g][0])
                            for g in range(n_groups)])
    accuracies, type_accuracies = [], []
    for rep in range(reps):
        if level == "type":
            shuffled = rng.permutation(group_label)
            fake = shuffled[group_codes]
        else:
            fake = rng.permutation(y_index)
        onehot = np.zeros((len(fake), n_classes), dtype=np.float32)
        onehot[np.arange(len(fake)), fake] = 1.0
        predicted = np.empty(len(fake), dtype=np.int32)
        for train, test in folds:
            node = fit_tree(x[train], onehot[train], weight[train], depth)
            predicted[test] = predict_tree(node, x[test])
        accuracies.append(float((predicted == fake).mean()))
        hits = 0
        for g in range(n_groups):
            members = group_codes == g
            vote = np.bincount(predicted[members], minlength=n_classes).argmax()
            hits += int(vote == fake[members][0])
        type_accuracies.append(hits / n_groups)
        if log_every and (rep + 1) % log_every == 0:
            print(f"      null {rep + 1}/{reps} (su ana kadar ortalama "
                  f"{np.mean(accuracies):.4f})", flush=True)
    accuracies = np.array(accuracies)
    type_accuracies = np.array(type_accuracies)
    return {"reps": reps, "max_depth": depth,
            "shuffled_at": ("type level, because the label is a property of the "
                            "type" if level == "type"
                            else "entry level, because this label varies within "
                                 "a type"),
            "mean": round(float(accuracies.mean()), 4),
            "sd": round(float(accuracies.std(ddof=1)) if reps > 1 else 0.0, 4),
            "p95": round(float(np.quantile(accuracies, 0.95)), 4),
            "max": round(float(accuracies.max()), 4),
            "type_level_mean": round(float(type_accuracies.mean()), 4),
            "type_level_p95": round(float(np.quantile(type_accuracies, 0.95)), 4),
            "type_level_max": round(float(type_accuracies.max()), 4),
            "_entry_values": accuracies.tolist(),
            "_type_values": type_accuracies.tolist()}


def empirical_p(observed, values):
    """Tek yanli ampirik p: null'da gozlenenden iyi veya esit kac tekrar var."""
    if not values:
        return None
    better = sum(1 for value in values if value >= observed)
    return round((better + 1) / (len(values) + 1), 4)


def verdict(result):
    """Sayilardan karar -- sabit sablon, yorum degil.

    Iki METRIK AYRI AYRI karara baglanir, cunku bu veri setinde ikisi farkli
    cevap verebiliyor: giris duzeyi dogruluk az sayida cok buyuk tipin
    agirligini tasidigi icin null'u genis; tip duzeyi dogrulukta her tip bir
    oy. Ikisini tek cumleye sikistirmak, biri anlamsiz digeri anlamliyken
    yanlis bir "anlamli" izlenimi verirdi.

    Sira onemli: bir dogruluk once null'u, SONRA klad kahinini gecmek zorunda.
    Null'u gecmek "rastgele degil" demektir; klad kahinini gecmek "derin klad
    uyeliginin otesinde bir sey var" demektir ve ikincisi cok daha agir bir
    iddiadir.
    """
    grouped = result["grouped_tree"]["accuracy"]
    grouped_type = result["grouped_tree"]["type_level_accuracy"]
    base = result["baselines"]
    compare = result.get("grouped_vs_null")
    if compare is None:
        return ("no permutation null was run for this secondary question, so "
                "only the baselines can be compared")
    parts = []
    for name, probability, value, majority, oracle in (
            ("per entry", compare["empirical_p"], grouped,
             base["global_majority_rate"], base["clade_oracle"]["accuracy"]),
            ("per type", compare["type_level_empirical_p"], grouped_type,
             base["type_level_majority_rate"],
             base["clade_oracle"]["type_level_accuracy"])):
        if value < majority:
            parts.append(f"{name}: below the majority-class rate "
                         f"({value:.3f} against {majority:.3f}), so the model "
                         f"is worse than always answering with the commonest "
                         f"class")
        elif probability > 0.05:
            parts.append(f"{name}: inside the label-permutation null "
                         f"(one-sided p = {probability}), so no signal is "
                         f"demonstrated")
        elif value <= oracle:
            parts.append(f"{name}: above the null (p = {probability}) but not "
                         f"above the clade oracle ({value:.3f} against "
                         f"{oracle:.3f}), so what is recovered is deep clade "
                         f"membership")
        else:
            parts.append(f"{name}: above the null (p = {probability}) and "
                         f"above the clade oracle ({value:.3f} against "
                         f"{oracle:.3f}), so the features add something clade "
                         f"membership alone does not give")
    return "; ".join(parts)


def nested_depth_cv(x, folds, y_index, onehot, weight, group_codes, depths,
                    n_folds_inner=5):
    """Derinligi test verisine BAKMADAN sec: her dis katin kendi ic gruplanmis
    CV'si. Boylece raporlanan sayi hiperparametre secimiyle sismiyor."""
    predicted = np.empty(len(y_index), dtype=np.int32)
    chosen = []
    for train, test in folds:
        inner_codes, _ = codes(group_codes[train].tolist())
        inner = group_kfold(inner_codes,
                            min(n_folds_inner, len(set(inner_codes.tolist()))))
        x_train = x[train]
        onehot_train = onehot[train]
        weight_train = weight[train]
        y_train = y_index[train]
        best_depth, best_score = depths[0], -1.0
        for depth in depths:
            inner_predicted = np.empty(len(train), dtype=np.int32)
            for inner_train, inner_test in inner:
                node = fit_tree(x_train[inner_train], onehot_train[inner_train],
                                weight_train[inner_train], depth)
                inner_predicted[inner_test] = predict_tree(node,
                                                           x_train[inner_test])
            value = float((inner_predicted == y_train).mean())
            if value > best_score:
                best_depth, best_score = depth, value
        node = fit_tree(x_train, onehot_train, weight_train, best_depth)
        predicted[test] = predict_tree(node, x[test])
        chosen.append(best_depth)
    return predicted, chosen


# --- Veri yukleme ---------------------------------------------------------

def read_csv_dict(path, key="cluster"):
    with open(path, newline="", encoding="utf-8") as handle:
        return {row[key]: row for row in csv.DictReader(handle)}


def load_entries(connection):
    """Dogrulanmis her RO icin kimlik, tip, klad ve replikon bilgisi."""
    rows = []
    for row in connection.execute("""
            SELECT r.candidate_id, r.ro_cluster, r.ro_group, r.nucleotide_id,
                   p.is_plasmid, p.organism, p.taxonomy
            FROM ro r LEFT JOIN replicon p USING(nucleotide_id)
            WHERE r.is_confirmed = 1 AND r.ro_cluster IS NOT NULL
              AND r.ro_cluster <> 'N/A'"""):
        rows.append({"candidate_id": row[0], "cluster": row[1],
                     "group": row[2] or "?", "nucleotide_id": row[3],
                     "is_plasmid": row[4], "organism": row[5] or "",
                     "taxonomy": row[6] or ""})
    return rows


def load_alignment_matrix(path, candidate_ids):
    """Hizalamanin match kolonlarini (n x kolon) byte matrisine cevir."""
    aligned = read_stockholm_matchcols(path)
    missing = [cid for cid in candidate_ids if cid not in aligned]
    if missing:
        raise SystemExit(f"[hata] {len(missing)} dogrulanmis RO hizalamada yok "
                         f"(ornek: {missing[:3]}); {path} bayat olabilir")
    joined = "".join(aligned[cid] for cid in candidate_ids)
    matrix = np.frombuffer(joined.encode("ascii"), dtype="S1")
    return matrix.reshape(len(candidate_ids), -1)


def onehot_columns(matrix, columns, min_count=MIN_RESIDUE_COUNT):
    """(kolon, kalinti) ikili ozellik matrisi. Bosluk da bir durumdur."""
    features, names = [], []
    n = matrix.shape[0]
    for column in columns:
        values, counts = np.unique(matrix[:, column - 1], return_counts=True)
        for value, count in zip(values, counts):
            if min_count <= count <= n - min_count:
                features.append(matrix[:, column - 1] == value)
                names.append((int(column), value.decode("ascii")))
    if not features:
        raise SystemExit("[hata] hicbir kolon ozellik esigini gecmedi")
    return np.stack(features, axis=1).astype(np.float32), names


def column_occupancy(matrix):
    return (matrix != b"-").mean(axis=0)


def sdp_scores(matrix, group_codes, columns, min_members=5, min_filled=0.5):
    """analyze_variants.py ile AYNI olcut: tipler arasi eksi tip ici entropi.

    Burada yeniden hesaplaniyor, cunku her katin yalnizca EGITIM tipleriyle
    secim yapmasi gerekiyor; hazir CSV tum tipleri gormus durumda.
    """
    by_group = defaultdict(list)
    for i, group in enumerate(group_codes.tolist()):
        by_group[group].append(i)
    total_groups = len(by_group)
    scores = {}
    for column in columns:
        values = matrix[:, column - 1]
        withins, consensus = [], Counter()
        filled = 0
        for indices in by_group.values():
            counter = Counter(v for v in values[indices] if v != b"-")
            if sum(counter.values()) < min_members:
                continue
            filled += 1
            withins.append(entropy(counter))
            consensus[counter.most_common(1)[0][0]] += 1
        if filled < total_groups * min_filled or not withins:
            continue
        scores[column] = entropy(consensus) - sum(withins) / len(withins)
    return scores


def build_context_features(connection, entries):
    """Genomik baglam ozellikleri; taksonomi AYRI blok olarak isaretlenir.

    Hicbiri proteinin DIZISINE bakmaz -- sorunun tamami bu: komsuluk kimyayi
    ongoruyor mu.
    """
    ids = [e["candidate_id"] for e in entries]
    position = {cid: i for i, cid in enumerate(ids)}
    n = len(ids)
    raw = defaultdict(lambda: [None] * n)

    for row in connection.execute("""
            SELECT candidate_id, has_beta, has_ferredoxin, has_reductase,
                   nearby_beta, nearby_ferredoxin, nearby_reductase,
                   completeness, n_genes FROM operon"""):
        i = position.get(row[0])
        if i is None:
            continue
        raw["operon_has_beta"][i] = int(bool(row[1]))
        raw["operon_has_ferredoxin"][i] = int(bool(row[2]))
        raw["operon_has_reductase"][i] = int(bool(row[3]))
        raw["nearby_beta"][i] = int(bool(row[4]))
        raw["nearby_ferredoxin"][i] = int(bool(row[5]))
        raw["nearby_reductase"][i] = int(bool(row[6]))
        raw["operon_complete"][i] = int(bool(row[7]))
        genes = int(row[8] or 0)
        raw["operon_size_class"][i] = ("1" if genes <= 1 else "2to3" if genes <= 3
                                       else "4to6" if genes <= 6 else "7plus")

    for row in connection.execute("""
            SELECT candidate_id, reductase_type, ferredoxin_type, has_beta
            FROM ro_etc"""):
        i = position.get(row[0])
        if i is None:
            continue
        raw["etc_reductase"][i] = row[1] or "none"
        raw["etc_ferredoxin"][i] = row[2] or "none"
        raw["etc_oligomer"][i] = "alpha3beta3" if row[3] else "alpha3"

    for row in connection.execute("""
            SELECT candidate_id, architecture, upstream_family,
                   upstream_category, upstream_divergent, intergenic_bp
            FROM ro_regulation"""):
        i = position.get(row[0])
        if i is None:
            continue
        raw["regulation_architecture"][i] = row[1] or "unknown"
        raw["regulator_family"][i] = row[2] or "none"
        raw["upstream_category"][i] = row[3] or "unknown"
        raw["upstream_divergent"][i] = int(bool(row[4]))
        gap = row[5]
        raw["intergenic_class"][i] = ("unknown" if gap is None
                                      else "le100" if gap <= 100
                                      else "101to400" if gap <= 400
                                      else "401to1000" if gap <= 1000
                                      else "gt1000")

    # Operon icindeki gen kategorileri ve bilesenleri -- varlik / yokluk.
    categories, components = defaultdict(set), defaultdict(set)
    for candidate_id, component, category in connection.execute(
            "SELECT candidate_id, component, category FROM operon_gene"):
        if candidate_id not in position:
            continue
        if category and category != "ro_alpha":
            categories[candidate_id].add(category)
        if component and component not in ("none", "alpha"):
            components[candidate_id].add(component)

    # +-10 kb penceresindeki bilesenler (operon disi komsuluk). Birlestirme
    # `protein_key` uzerinden: build_operons.py o anahtari
    # "nucleotide_id:start-end:strand" olarak kuruyor, yani SQL'de ifade
    # birlestirmesi yerine burada sozlukle eslenmesi daha ucuz ve daha acik.
    component_of = {key: value for key, value in connection.execute(
        "SELECT protein_key, component FROM neighbor_component "
        "WHERE component IS NOT NULL AND component <> 'none'")}
    window = defaultdict(set)
    for candidate_id, nucleotide_id, start, end, strand in connection.execute(
            "SELECT candidate_id, nucleotide_id, start, end, strand "
            "FROM neighbor"):
        if candidate_id not in position:
            continue
        component = component_of.get(f"{nucleotide_id}:{start}-{end}:{strand}")
        if component:
            window[candidate_id].add(component)

    features, names, blocks = [], [], []

    def add(name, vector, block):
        features.append(np.asarray(vector, dtype=np.float32))
        names.append(name)
        blocks.append(block)

    regulation_keys = ("regulation_architecture", "regulator_family",
                       "upstream_category", "upstream_divergent",
                       "intergenic_class")
    for key in sorted(raw):
        values = raw[key]
        block = "regulation" if key in regulation_keys else "operon"
        if all(v is None or isinstance(v, int) for v in values):
            add(f"{key}=1", [1.0 if v else 0.0 for v in values], block)
        else:
            for level in sorted({v for v in values if v is not None}):
                add(f"{key}={level}",
                    [1.0 if v == level else 0.0 for v in values], block)
    for level in sorted({c for s in categories.values() for c in s}):
        add(f"operon_gene_category={level}",
            [1.0 if level in categories.get(cid, ()) else 0.0 for cid in ids],
            "operon")
    for level in sorted({c for s in components.values() for c in s}):
        add(f"operon_component={level}",
            [1.0 if level in components.get(cid, ()) else 0.0 for cid in ids],
            "operon")
    for level in sorted({c for s in window.values() for c in s}):
        add(f"window_component={level}",
            [1.0 if level in window.get(cid, ()) else 0.0 for cid in ids],
            "neighbourhood")
    add("replicon_is_plasmid=1",
        [1.0 if e["is_plasmid"] else 0.0 for e in entries], "mobility")

    # Taksonomi AYRI blok: ablasyonla cikarilabilsin.
    lineages = [[part.strip() for part in e["taxonomy"].split(";") if part.strip()]
                for e in entries]
    for rank in (1, 2, 3):
        counter = Counter(lineage[rank] if len(lineage) > rank else "unknown"
                          for lineage in lineages)
        for level, count in sorted(counter.items()):
            if count >= 50:
                add(f"lineage_rank{rank}={level}",
                    [1.0 if (len(l) > rank and l[rank] == level) else 0.0
                     for l in lineages], "taxonomy")
    genera = [e["organism"].split()[0] if e["organism"] else "unknown"
              for e in entries]
    for level, count in sorted(Counter(genera).items()):
        if count >= 25:
            add(f"genus={level}", [1.0 if g == level else 0.0 for g in genera],
                "taxonomy")

    return np.stack(features, axis=1), names, blocks


# --- Bir basligin tam degerlendirmesi ------------------------------------

def evaluate_task(x, names, y_labels, group_codes, clade_codes, class_names,
                  rng, args, label="task", depth=PRIMARY_DEPTH, null_reps=0,
                  null_level="type", with_logo=False, with_logistic=False,
                  with_depth_sweep=False, with_nested=False, with_nn=False,
                  sequence=True):
    index_of = {name: i for i, name in enumerate(class_names)}
    y_index = np.array([index_of[v] for v in y_labels], dtype=np.int32)
    n_classes = len(class_names)
    onehot = np.zeros((len(y_index), n_classes), dtype=np.float32)
    onehot[np.arange(len(y_index)), y_index] = 1.0

    sizes = Counter(group_codes.tolist())
    weight = np.array([1.0 / sizes[g] for g in group_codes.tolist()],
                      dtype=np.float32)
    weight_plain = np.ones(len(y_index), dtype=np.float32)

    grouped = group_kfold(group_codes, args.group_folds)
    random_folds = random_kfold(len(y_index), N_RANDOM_FOLDS,
                                np.random.RandomState(SEED))

    def run_tree(folds, sample_weight, tree_depth):
        predicted = np.empty(len(y_index), dtype=np.int32)
        for train, test in folds:
            node = fit_tree(x[train], onehot[train], sample_weight[train],
                            tree_depth)
            predicted[test] = predict_tree(node, x[test])
        return predicted

    type_labels = {}
    for group, value in zip(group_codes.tolist(), y_labels):
        type_labels.setdefault(group, value)

    result = {"n_entries": int(len(y_index)), "n_types": len(sizes),
              "n_features": int(x.shape[1]), "classes": list(class_names),
              "class_entry_counts": {c: int((y_index == i).sum())
                                     for i, c in enumerate(class_names)},
              "class_type_counts": dict(Counter(type_labels.values())),
              "tree_max_depth": depth}

    print(f"  [{label}] gruplanmis agac (derinlik {depth}) ...", flush=True)
    result["grouped_tree"] = score(run_tree(grouped, weight, depth), y_index,
                                   group_codes, class_names)
    result["grouped_tree"]["scheme"] = (
        f"GroupKFold, {args.group_folds} entry-balanced folds, "
        f"group = {GROUP_FIELD}")

    print(f"  [{label}] SIZINTILI rastgele sema ...", flush=True)
    result["leaked_random_tree"] = score(run_tree(random_folds, weight, depth),
                                         y_index, group_codes, class_names)
    result["leaked_random_tree"]["scheme"] = (
        f"random {N_RANDOM_FOLDS}-fold split over entries; INVALID here because "
        "near-identical sequences of the same type land on both sides and the "
        "label is a property of the type")
    result["leakage_gap"] = {
        "grouped_accuracy": result["grouped_tree"]["accuracy"],
        "leaked_accuracy": result["leaked_random_tree"]["accuracy"],
        "absolute_gap": round(result["leaked_random_tree"]["accuracy"]
                              - result["grouped_tree"]["accuracy"], 4),
        "leaked_over_grouped": round(
            result["leaked_random_tree"]["accuracy"]
            / max(result["grouped_tree"]["accuracy"], 1e-9), 3)}

    result["baselines"] = majority_baselines(grouped, y_index, group_codes,
                                             class_names)
    oracle = clade_oracle(grouped, group_codes, clade_codes, y_index, n_classes)
    result["baselines"]["clade_oracle"] = score(oracle, y_index, group_codes,
                                                class_names)
    result["baselines"]["clade_oracle"]["what_it_is"] = (
        "predict the majority label of the training TYPES that share the HMM "
        "group (ro.ro_group) of the test entry. ro_group comes from sequence "
        "similarity, so it is available for a new protein; if the model does "
        "not beat this, the alignment columns add nothing beyond deep clade "
        "membership")

    if with_nn:
        print(f"  [{label}] 1-NN sizinti gosterimi ...", flush=True)
        grouped_nn = nearest_neighbour_accuracy(x, y_index, grouped)
        leaked_nn = nearest_neighbour_accuracy(x, y_index, random_folds)
        result["nearest_neighbour_leakage_demo"] = {
            "grouped_accuracy": round(grouped_nn, 4),
            "leaked_random_accuracy": round(leaked_nn, 4),
            "leaked_over_grouped": round(leaked_nn / max(grouped_nn, 1e-9), 3),
            "what_it_shows": (
                "a one-nearest-neighbour classifier has nothing to learn: it "
                "copies the label of the most similar training row. On a random "
                "split it is close to perfect, because almost every test row has "
                "a near-identical twin of the same type in training. Hold the "
                "type out and the same classifier collapses. This is the "
                "leakage mechanism in its plainest form, and it is the reason "
                "every other number in this file is grouped by type")}

    plain = run_tree(grouped, weight_plain, depth)
    result["sensitivity_unweighted"] = {
        "accuracy": round(float((plain == y_index).mean()), 4),
        "note": ("entries weighted equally instead of 1/type-size, so the "
                 "largest types dominate training")}

    if with_depth_sweep:
        sweep = []
        for tree_depth in DEPTH_SWEEP:
            grouped_accuracy = float((run_tree(grouped, weight, tree_depth)
                                      == y_index).mean())
            leaked_accuracy = float((run_tree(random_folds, weight, tree_depth)
                                     == y_index).mean())
            sweep.append({"max_depth": tree_depth,
                          "grouped_accuracy": round(grouped_accuracy, 4),
                          "leaked_random_accuracy": round(leaked_accuracy, 4),
                          "ratio": round(leaked_accuracy
                                         / max(grouped_accuracy, 1e-9), 3)})
            print(f"      depth {tree_depth:>2}: grouped {grouped_accuracy:.4f} "
                  f"leaked {leaked_accuracy:.4f}", flush=True)
        result["depth_sweep"] = {
            "rows": sweep,
            "note": ("the gap between the two columns is what a random split "
                     "buys by memorising types: it widens with model capacity, "
                     "which is the signature of leakage rather than of "
                     "learning. The best grouped row of this sweep must NOT be "
                     "read as an achieved accuracy, because it was chosen by "
                     "looking at the held-out score; see nested_depth_selection")}

    if with_nested:
        print(f"  [{label}] ic ice derinlik secimi ...", flush=True)
        nested, chosen = nested_depth_cv(x, grouped, y_index, onehot, weight,
                                         group_codes, NESTED_DEPTHS)
        result["nested_depth_selection"] = score(nested, y_index, group_codes,
                                                 class_names)
        result["nested_depth_selection"]["depths_chosen_per_fold"] = chosen
        result["nested_depth_selection"]["what_it_is"] = (
            "the honest number: in every outer fold the tree depth was chosen "
            "by an inner grouped cross-validation over the training types only, "
            "so the held-out types took no part in the choice")

    if with_logo:
        print(f"  [{label}] her tipi sirayla disarida birak ...", flush=True)
        result["leave_one_type_out_tree"] = score(
            run_tree(leave_one_group_out(group_codes), weight, depth), y_index,
            group_codes, class_names)
        result["leave_one_type_out_tree"]["scheme"] = (
            "leave-one-type-out: each type is predicted by a model that never "
            "saw a single member of it")

    if with_logistic:
        print(f"  [{label}] L1 lojistik ...", flush=True)
        predicted = np.empty(len(y_index), dtype=np.int32)
        for train, test in grouped:
            model = fit_l1_logistic(x[train], y_index[train], weight[train],
                                    n_classes, args.l1_penalty,
                                    args.l1_iterations)
            predicted[test] = predict_l1_logistic(model, x[test])
        result["grouped_l1_logistic"] = score(predicted, y_index, group_codes,
                                              class_names)
        result["grouped_l1_logistic"]["penalty"] = args.l1_penalty
        leaked = np.empty(len(y_index), dtype=np.int32)
        for train, test in random_folds:
            model = fit_l1_logistic(x[train], y_index[train], weight[train],
                                    n_classes, args.l1_penalty,
                                    args.l1_iterations)
            leaked[test] = predict_l1_logistic(model, x[test])
        result["leaked_random_l1_logistic"] = {
            "accuracy": round(float((leaked == y_index).mean()), 4)}

    if null_reps > 0:
        print(f"  [{label}] permutasyon null'i ({null_reps}) ...", flush=True)
        started = time.time()
        null = permutation_null(x, grouped, group_codes, y_index, n_classes,
                                weight, depth, null_reps, rng, level=null_level)
        entry_values = null.pop("_entry_values")
        type_values = null.pop("_type_values")
        null["seconds"] = round(time.time() - started, 1)
        result["permutation_null"] = null
        observed = result["grouped_tree"]["accuracy"]
        observed_type = result["grouped_tree"]["type_level_accuracy"]
        result["grouped_vs_null"] = {
            "grouped_accuracy": observed,
            "null_mean": null["mean"], "null_sd": null["sd"],
            "null_max": null["max"],
            "z_against_null": (round((observed - null["mean"]) / null["sd"], 2)
                               if null["sd"] > 0 else None),
            "empirical_p": empirical_p(observed, entry_values),
            "type_level_accuracy": observed_type,
            "type_level_null_mean": null["type_level_mean"],
            "type_level_empirical_p": empirical_p(observed_type, type_values),
            "note": ("the empirical p is one-sided over the permutation "
                     "replicates. With 61 types of very unequal size the null "
                     "spread is wide, which is exactly why an accuracy without "
                     "a null is not interpretable here")}

    full_tree = fit_tree(x, onehot, weight, min(depth, 3))
    result["tree_rules_depth3_fitted_on_all_data"] = tree_rules(
        full_tree, names, class_names, max_depth=3, sequence=sequence)
    result["tree_rules_note"] = (
        "printed from a depth-3 tree fitted on every entry, so it is a "
        "description of the data and NOT a performance estimate")
    result["verdict"] = verdict(result)
    return result, onehot, weight, grouped


# --- Degisken onemi ve motif haritasi ------------------------------------

def motif_relation(column, motif_columns):
    """Kolonun sekiz tanimlayici kolona gore konumu -- biyolog bunu sorar."""
    nearest = min(motif_columns, key=lambda m: abs(m["column"] - column))
    offset = column - nearest["column"]
    return {"zone": zone_of(column),
            "is_motif_column": any(m["column"] == column for m in motif_columns),
            "nearest_catalytic_centre_column": nearest["column"],
            "nearest_catalytic_centre_label": nearest["label"],
            "offset_from_nearest": offset,
            "within_10_columns_of_centre": abs(offset) <= 10}


def importance_block(x, names, onehot, weight, grouped, rng, args,
                     motif_columns, matrix=None, group_codes=None,
                     y_labels=None, sequence=True, top=25):
    n_features = x.shape[1]
    n_classes = onehot.shape[1]
    y_index = np.argmax(onehot, axis=1).astype(np.int32)

    gains = np.zeros(n_features)
    for train, _ in grouped:
        node = fit_tree(x[train], onehot[train], weight[train], PRIMARY_DEPTH)
        gains += tree_gains(node, n_features)
    gains /= max(len(grouped), 1)

    model = fit_l1_logistic(x, y_index, weight, n_classes, args.l1_penalty,
                            args.l1_iterations)
    l1_weight = np.abs(model["coef"]).sum(axis=1)
    nonzero = int((l1_weight > 1e-8).sum())

    forest = forest_importances(x, onehot, weight, args.forest_trees,
                                PRIMARY_DEPTH + 1, rng)

    def rank(values):
        return {int(i): int(r) + 1 for r, i in enumerate(np.argsort(-values))}

    rank_tree, rank_l1, rank_forest = rank(gains), rank(l1_weight), rank(forest)
    combined = np.array([rank_tree[i] + rank_l1[i] + rank_forest[i]
                         for i in range(n_features)])
    rows = []
    for i in np.argsort(combined)[:top]:
        i = int(i)
        entry = {"feature": (f"column {names[i][0]} residue {names[i][1]}"
                             if sequence else names[i]),
                 "tree_gain": round(float(gains[i]), 5),
                 "l1_abs_weight": round(float(l1_weight[i]), 5),
                 "forest_gain": round(float(forest[i]), 5),
                 "rank_tree": rank_tree[i], "rank_l1": rank_l1[i],
                 "rank_forest": rank_forest[i]}
        if sequence:
            entry["column"] = names[i][0]
            entry["residue"] = names[i][1]
            # Bosluk durumu bir KALINTI degil: proteinin o kolonda hic harfi
            # olmadigini soyler, yani uzunluk ya da hizalama bilgisi tasir.
            # Ayri isaretleniyor, cunku biyolojik yorumu bambaskadir.
            entry["is_gap_state"] = names[i][1] == "-"
            entry.update(motif_relation(names[i][0], motif_columns))
        rows.append(entry)

    block = {"top_features": rows, "l1_nonzero_features": nonzero,
             "l1_total_features": n_features,
             "method": ("features ranked by the sum of three ranks: mean split "
                        "gain of the depth-limited tree over the grouped "
                        "training folds, absolute L1 logistic weight summed over "
                        "classes, and mean split gain of a random forest. The "
                        "forest is fitted ONLY to rank features; every accuracy "
                        "in this file comes from the tree or the L1 model, both "
                        "of which can be read")}

    if sequence and matrix is not None:
        per_column = []
        type_label = {}
        for group, label in zip(group_codes.tolist(), y_labels):
            type_label.setdefault(group, label)
        class_names = sorted(set(y_labels))
        group_masks = {g: (group_codes == g) for g in sorted(type_label)}
        for column in sorted({names[i][0] for i in range(n_features)}):
            values = matrix[:, column - 1]
            consensus = {}
            for group, mask in group_masks.items():
                counter = Counter(v for v in values[mask] if v != b"-")
                if counter:
                    consensus[group] = counter.most_common(1)[0][0].decode("ascii")
            residues = sorted(set(consensus.values()))
            table = [[sum(1 for g, r in consensus.items()
                          if r == residue and type_label[g] == name)
                      for name in class_names] for residue in residues]
            per_column.append({"column": column,
                               "n_types_scored": len(consensus),
                               "n_residue_states": len(residues),
                               "cramers_v_type_level": round(cramers_v(table), 3),
                               **motif_relation(column, motif_columns)})
        per_column.sort(key=lambda r: -r["cramers_v_type_level"])
        block["univariate_type_level"] = {
            "note": ("one row per alignment column. The residue of a TYPE is the "
                     "consensus of its members, so each type counts once and a "
                     "large sequenced clade cannot inflate the association. "
                     "Cramer's V on a 61-type table is noisy and is given for "
                     "ordering, not as a test"),
            "rows": per_column[:top]}
    return block


def load_motif_columns(path):
    """motif_stats.json'daki sekiz kolon; dosya yoksa ro_motif.py sabitleri."""
    if os.path.exists(path):
        try:
            with open(path, encoding="utf-8") as handle:
                payload = json.load(handle)
            sites = [{"column": int(s["column"]), "label": s["label"],
                      "expected": s["expected"],
                      "conserved_fraction": s.get("fraction")}
                     for s in payload.get("sites", [])]
            if sites:
                return sites, "analysis_out/motif_stats.json"
        except (ValueError, OSError, KeyError):
            pass
    model = MODEL_COLUMNS[DEFAULT_MODEL]
    sites = [{"column": c, "label": label, "expected": r,
              "conserved_fraction": None}
             for c, r, label in model["rieske"] + model["catalytic"]]
    sites.append({"column": model["bridging"][0],
                  "label": "inter-subunit bridging aspartate",
                  "expected": model["bridging"][1], "conserved_fraction": None})
    return sites, "ro_motif.py constants"


# --- Kolon secimi semalari -----------------------------------------------

def selection_variants(matrix, group_codes, y_labels, class_names, args,
                       precomputed_columns, occupancy):
    """Uc kolon secim semasinin gruplanmis dogrulugu -- secim sizintisinin bedeli.

    (1) hazir sdp_positions.csv: kolonlar TUM tipler gorulerek secilmis
    (2) ic-kat: ayni olcut her katin yalnizca egitim tipleriyle
    (3) secim yok: dolu kolonlarin tamami
    """
    index_of = {name: i for i, name in enumerate(class_names)}
    y_index = np.array([index_of[v] for v in y_labels], dtype=np.int32)
    n_classes = len(class_names)
    onehot = np.zeros((len(y_index), n_classes), dtype=np.float32)
    onehot[np.arange(len(y_index)), y_index] = 1.0
    sizes = Counter(group_codes.tolist())
    weight = np.array([1.0 / sizes[g] for g in group_codes.tolist()],
                      dtype=np.float32)
    grouped = group_kfold(group_codes, args.group_folds)
    occupied = [int(c) + 1 for c in np.where(occupancy >= args.min_occupancy)[0]]
    rows = []

    def accuracy_fixed(columns):
        x, _ = onehot_columns(matrix, columns)
        predicted = np.empty(len(y_index), dtype=np.int32)
        for train, test in grouped:
            node = fit_tree(x[train], onehot[train], weight[train], PRIMARY_DEPTH)
            predicted[test] = predict_tree(node, x[test])
        return float((predicted == y_index).mean())

    fixed = sorted(precomputed_columns[:args.top_sdp])
    rows.append({"selection": "precomputed_sdp_positions_csv",
                 "n_columns": len(fixed),
                 "grouped_accuracy": round(accuracy_fixed(fixed), 4),
                 "selection_leak": ("yes -- the ranking in sdp_positions.csv was "
                                    "computed with every type visible, including "
                                    "the types held out in each fold")})
    print(f"      precomputed: {rows[-1]['grouped_accuracy']}", flush=True)

    predicted = np.empty(len(y_index), dtype=np.int32)
    overlaps = []
    for train, test in grouped:
        scores = sdp_scores(matrix[train], group_codes[train], occupied)
        ranked = sorted(sorted(scores), key=lambda c: -scores[c])[:args.top_sdp]
        overlaps.append(len(set(ranked) & set(fixed)) / max(len(fixed), 1))
        x, _ = onehot_columns(matrix, sorted(ranked))
        node = fit_tree(x[train], onehot[train], weight[train], PRIMARY_DEPTH)
        predicted[test] = predict_tree(node, x[test])
    rows.append({"selection": "sdp_recomputed_inside_each_training_fold",
                 "n_columns": args.top_sdp,
                 "grouped_accuracy": round(float((predicted == y_index).mean()),
                                           4),
                 "mean_overlap_with_precomputed_list": round(
                     float(np.mean(overlaps)), 3),
                 "selection_leak": ("no -- the held-out types take no part in "
                                    "the scoring")})
    print(f"      in-fold: {rows[-1]['grouped_accuracy']}", flush=True)

    rows.append({"selection": "no_selection_all_occupied_columns",
                 "n_columns": len(occupied),
                 "grouped_accuracy": round(accuracy_fixed(occupied), 4),
                 "selection_leak": ("no -- no label was consulted; occupancy is "
                                    "the only filter")})
    print(f"      no selection: {rows[-1]['grouped_accuracy']}", flush=True)
    return {"rows": rows, "occupancy_threshold": args.min_occupancy,
            "note": ("if the precomputed row were far above the in-fold row, the "
                     "headline accuracy would be partly an artefact of choosing "
                     "the columns on the full data set")}


# --- Ana akis -------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--db", default="roar.sqlite")
    parser.add_argument("--out-dir", default="analysis_out")
    parser.add_argument("--alignment", default="genomic_context/cand_aln.sto")
    parser.add_argument("--chemistry", default="chemistry.csv")
    parser.add_argument("--ecology", default="cluster_ecology.csv")
    parser.add_argument("--sdp", default="analysis_out/sdp_positions.csv")
    parser.add_argument("--motif-stats", default="analysis_out/motif_stats.json")
    parser.add_argument("--top-sdp", type=int, default=60,
                        help="kac SDP kolonu ozellik havuzuna girsin")
    parser.add_argument("--min-occupancy", type=float, default=0.9,
                        help="secimsiz havuz icin kolon doluluk esigi")
    parser.add_argument("--group-folds", type=int, default=N_GROUP_FOLDS)
    parser.add_argument("--null-reps", type=int, default=200)
    parser.add_argument("--forest-trees", type=int, default=60)
    parser.add_argument("--l1-penalty", type=float, default=0.002)
    parser.add_argument("--l1-iterations", type=int, default=150)
    parser.add_argument("--skip-selection-variants", action="store_true")
    args = parser.parse_args()

    started = time.time()
    rng = np.random.RandomState(SEED)
    connection = sqlite3.connect(args.db)

    print("[okunuyor] etiketler ve girisler")
    chemistry = read_csv_dict(args.chemistry)
    ecology = read_csv_dict(args.ecology)
    entries = [e for e in load_entries(connection)
               if e["cluster"] in chemistry and e["cluster"] in ecology]
    candidate_ids = [e["candidate_id"] for e in entries]
    group_codes, group_names = codes([e["cluster"] for e in entries])
    clade_codes, clade_names = codes([e["group"] for e in entries])
    print(f"[bilgi] {len(entries)} giris, {len(group_names)} tip, "
          f"{len(clade_names)} HMM grubu")

    print("[okunuyor] hizalama")
    matrix = load_alignment_matrix(args.alignment, candidate_ids)
    occupancy = column_occupancy(matrix)
    print(f"[bilgi] {matrix.shape[1]} match kolonu, "
          f"{int((occupancy >= args.min_occupancy).sum())} tanesi "
          f">= {args.min_occupancy} dolu")

    motif_columns, motif_source = load_motif_columns(args.motif_stats)
    with open(args.sdp, newline="", encoding="utf-8") as handle:
        precomputed = [int(row["column"]) for row in csv.DictReader(handle)]
    sequence_columns = sorted(precomputed[:args.top_sdp])
    x_sequence, sequence_names = onehot_columns(matrix, sequence_columns)
    print(f"[bilgi] dizi ozellikleri: {x_sequence.shape[1]} ikili ozellik, "
          f"{len(sequence_columns)} kolon")

    print("[kuruluyor] genomik baglam ozellikleri")
    x_context, context_names, context_blocks = build_context_features(
        connection, entries)
    blocks = np.array(context_blocks)
    no_taxonomy = blocks != "taxonomy"
    only_taxonomy = ~no_taxonomy
    context_only_names = [n for n, keep in zip(context_names, no_taxonomy) if keep]
    print(f"[bilgi] baglam ozellikleri: {x_context.shape[1]} "
          f"({int(no_taxonomy.sum())} baglam + {int(only_taxonomy.sum())} "
          f"taksonomi)")

    clusters = [e["cluster"] for e in entries]
    reaction = [chemistry[c]["reaction_class"] for c in clusters]
    substrate = [ecology[c]["substrate_class"] for c in clusters]
    ro_group = [e["group"] for e in entries]
    reaction_classes = sorted(set(reaction))
    substrate_classes = sorted(set(substrate))

    results = {}

    print("\n=== Soru 1: reaksiyon sinifi <- dizi")
    q1, onehot1, weight1, grouped1 = evaluate_task(
        x_sequence, sequence_names, reaction, group_codes, clade_codes,
        reaction_classes, rng, args, label="reaction<-sequence",
        null_reps=args.null_reps, with_logo=True, with_logistic=True,
        with_depth_sweep=True, with_nested=True, with_nn=True)
    q1["question"] = ("can the reaction class of an enzyme be predicted from its "
                      "sequence alone, for a type the model has never seen")
    q1["features"] = {
        "kind": "residue identity at alignment match-state columns",
        "alignment": args.alignment,
        "column_pool": (f"top {args.top_sdp} rows of {args.sdp}, the "
                        "specificity-determining positions the pipeline already "
                        "computes"),
        "encoding": ("one binary feature per (column, residue) pair seen in at "
                     f"least {MIN_RESIDUE_COUNT} entries; the gap character is a "
                     "state of its own")}
    results["q1_reaction_class_from_sequence"] = q1

    print("\n=== Soru 2: hangi kolonlar")
    q2 = importance_block(x_sequence, sequence_names, onehot1, weight1, grouped1,
                          rng, args, motif_columns, matrix=matrix,
                          group_codes=group_codes, y_labels=reaction,
                          sequence=True)
    q2["catalytic_centre_columns"] = motif_columns
    q2["catalytic_centre_source"] = motif_source
    q2["zone_convention"] = {
        "rieske_domain": f"column <= {RIESKE_END}",
        "linker": f"{RIESKE_END} < column < {CATALYTIC_START}",
        "catalytic_domain": f"column >= {CATALYTIC_START}",
        "note": ("the boundaries are the convention of variant_signature.py and "
                 "are approximate; no structure was used, so a column close to a "
                 "metal ligand is a hypothesis about the active site and not "
                 "evidence of one")}
    q2["zone_counts_in_top_features"] = dict(
        Counter(row["zone"] for row in q2["top_features"]))
    q2["top_features_within_10_columns_of_a_centre_column"] = sum(
        1 for row in q2["top_features"] if row["within_10_columns_of_centre"])
    q2["zone_counts_in_whole_pool"] = dict(
        Counter(zone_of(c) for c in sequence_columns))
    q2["gap_state_features_in_top"] = sum(
        1 for row in q2["top_features"] if row["is_gap_state"])
    # Katalitik bolge zenginlesmesi, HAVUZA GORE. Havuz zaten katalitik
    # kolonlardan olustugu icin ham sayi tek basina bir sey soylemiyor.
    top_catalytic = q2["zone_counts_in_top_features"].get("catalytic_domain", 0)
    pool_catalytic = q2["zone_counts_in_whole_pool"].get("catalytic_domain", 0)
    top_share = top_catalytic / max(len(q2["top_features"]), 1)
    pool_share = pool_catalytic / max(len(sequence_columns), 1)
    q2["catalytic_domain_share_in_top_features"] = round(top_share, 3)
    q2["catalytic_domain_share_in_whole_pool"] = round(pool_share, 3)
    q2["catalytic_domain_enrichment_over_pool"] = round(
        top_share / pool_share, 3) if pool_share else None
    q2["how_to_read_the_offsets"] = (
        "the zone counts of the top features have to be read against the zone "
        "counts of the whole pool: the pool is already dominated by catalytic "
        "domain columns, so a catalytic majority among the top features is "
        "only informative if it exceeds the pool proportion. Two traps are "
        "flagged separately. First, a feature one or two columns away from a "
        "motif column may be the SAME residue displaced by an alignment shift "
        "rather than a different position: motif_stats.json records that the "
        "iron carboxylate sits at column 355 in 9,936 entries but is offset or "
        "absent in the rest, and the motif test itself uses a window of plus or "
        "minus two columns. Second, a gap state says the protein has no letter "
        "in that column at all, so it encodes length or alignment quality and "
        "not residue chemistry; is_gap_state marks those rows")
    results["q2_informative_columns"] = q2

    if not args.skip_selection_variants:
        print("\n=== Kolon secim semalari")
        results["q1_column_selection_variants"] = selection_variants(
            matrix, group_codes, reaction, reaction_classes, args, precomputed,
            occupancy)

    print("\n=== Soru 3: substrat sinifi <- genomik baglam")
    q3, onehot3, weight3, grouped3 = evaluate_task(
        x_context[:, no_taxonomy], context_only_names, substrate, group_codes,
        clade_codes, substrate_classes, rng, args, label="substrate<-context",
        null_reps=args.null_reps, with_logo=True, with_logistic=True,
        with_nested=True, with_nn=True, sequence=False)
    q3["question"] = ("does the genomic neighbourhood of the gene predict the "
                      "chemical class of its substrate, for a type never seen "
                      "in training")
    q3["features"] = {
        "kind": "genomic context only; no residue of the protein is used",
        "blocks": dict(Counter(blocks[no_taxonomy].tolist())),
        "sources": ["operon", "operon_gene", "neighbor_component", "ro_etc",
                    "ro_regulation", "replicon.is_plasmid"]}
    results["q3_substrate_class_from_context"] = q3

    print("\n--- Soru 3 ablasyonlari")
    ablations = {}
    for name, mask in (("context_without_taxonomy", no_taxonomy),
                       ("taxonomy_only", only_taxonomy),
                       ("context_plus_taxonomy", np.ones_like(no_taxonomy))):
        subset, _, _, _ = evaluate_task(
            x_context[:, mask],
            [n for n, keep in zip(context_names, mask) if keep], substrate,
            group_codes, clade_codes, substrate_classes, rng, args,
            label=f"ablation:{name}", sequence=False)
        ablations[name] = {
            "n_features": int(mask.sum()),
            "grouped_accuracy": subset["grouped_tree"]["accuracy"],
            "grouped_balanced_accuracy":
                subset["grouped_tree"]["balanced_accuracy"],
            "type_level_accuracy": subset["grouped_tree"]["type_level_accuracy"],
            "leaked_random_accuracy": subset["leaked_random_tree"]["accuracy"]}
    ablations["note"] = (
        "taxonomy and neighbourhood are entangled: a lineage carries both its "
        "gene repertoire and its habitat. If taxonomy alone reaches what the "
        "neighbourhood reaches, then 'the neighbourhood predicts the chemistry' "
        "is not established -- what is established is that the lineage does")
    results["q3_ablations"] = ablations

    results["q3_informative_context_features"] = importance_block(
        x_context[:, no_taxonomy], context_only_names, onehot3, weight3,
        grouped3, rng, args, motif_columns, sequence=False)

    print("\n=== Soru 4a: POZITIF KONTROL -- HMM grubu <- dizi")
    positive, _, _, _ = evaluate_task(
        x_sequence, sequence_names, ro_group, group_codes, clade_codes,
        sorted(set(ro_group)), rng, args, label="positive-control",
        null_reps=max(50, args.null_reps // 4))
    positive["why_this_control"] = (
        "ro_group is the coarse HMM group of the type, so it is itself derived "
        "from sequence similarity. A model that recovers it across held-out "
        "types proves the grouped scheme is not blind, so a low accuracy "
        "elsewhere says something about the biology and not about the method. "
        "It is NOT an independent biological result, and its clade-oracle "
        "baseline is 1.0 by construction")
    results["q4a_positive_control_ro_group_from_sequence"] = positive

    print("\n=== Soru 4b: modalite karsilastirmasi")
    modality = {}
    for name, data, feature_names, target, class_names, is_sequence in (
            ("reaction_class_from_sequence", x_sequence, sequence_names,
             reaction, reaction_classes, True),
            ("reaction_class_from_context", x_context[:, no_taxonomy],
             context_only_names, reaction, reaction_classes, False),
            ("substrate_class_from_sequence", x_sequence, sequence_names,
             substrate, substrate_classes, True),
            ("substrate_class_from_context", x_context[:, no_taxonomy],
             context_only_names, substrate, substrate_classes, False)):
        subset, _, _, _ = evaluate_task(
            data, feature_names, target, group_codes, clade_codes, class_names,
            rng, args, label=f"modality:{name}", sequence=is_sequence)
        modality[name] = {
            "grouped_accuracy": subset["grouped_tree"]["accuracy"],
            "grouped_balanced_accuracy":
                subset["grouped_tree"]["balanced_accuracy"],
            "type_level_accuracy": subset["grouped_tree"]["type_level_accuracy"],
            "global_majority_rate": subset["baselines"]["global_majority_rate"],
            "type_level_majority_rate":
                subset["baselines"]["type_level_majority_rate"],
            "clade_oracle_accuracy":
                subset["baselines"]["clade_oracle"]["accuracy"],
            "clade_oracle_type_level":
                subset["baselines"]["clade_oracle"]["type_level_accuracy"],
            "leaked_random_accuracy": subset["leaked_random_tree"]["accuracy"]}
    modality["note"] = ("the same folds, the same model and the same baselines "
                        "in every cell, so the four rows can be compared with "
                        "one another")
    results["q4b_modality_comparison"] = modality

    print("\n=== Soru 4c: tipin ozelligi OLMAYAN bir etiket (beta alt birimi)")
    beta_lookup = dict(connection.execute(
        "SELECT candidate_id, has_beta FROM ro_etc"))
    beta_label = ["alpha3beta3" if beta_lookup.get(e["candidate_id"]) else "alpha3"
                  for e in entries]
    beta, _, _, _ = evaluate_task(
        x_sequence, sequence_names, beta_label, group_codes, clade_codes,
        ["alpha3", "alpha3beta3"], rng, args, label="beta-subunit",
        null_reps=max(50, args.null_reps // 4), null_level="entry")
    beta["why_this_question"] = (
        "every other label in this file is a property of the type, so the "
        "effective sample size is the number of types. Whether a beta subunit "
        "sits next to the gene varies WITHIN a type (mean within-type entropy "
        "0.27 bit against 0.76 bit overall), so this is the one question where "
        "the individual entries carry information of their own. It also has a "
        "structural meaning: alpha3 against alpha3beta3 oligomer architecture")
    beta["limit"] = (
        "the label comes from the annotation of the neighbouring gene, so a "
        "missing beta may mean a missing annotation rather than a missing "
        "subunit. The permutation null for this question shuffles at ENTRY "
        "level, because the label is not a property of the type")
    results["q4c_beta_subunit_from_sequence"] = beta

    # Ozet blok: her baslik icin okunacak sayilar bir arada. Hicbir sayi elle
    # yazilmiyor, hepsi yukaridaki sonuclardan okunuyor.
    summary = {}
    for key in ("q1_reaction_class_from_sequence",
                "q3_substrate_class_from_context",
                "q4a_positive_control_ro_group_from_sequence",
                "q4c_beta_subunit_from_sequence"):
        task = results[key]
        null = task.get("permutation_null", {})
        summary[key] = {
            "grouped_accuracy": task["grouped_tree"]["accuracy"],
            "grouped_type_level_accuracy":
                task["grouped_tree"]["type_level_accuracy"],
            "grouped_balanced_accuracy":
                task["grouped_tree"]["balanced_accuracy"],
            "leaked_random_accuracy": task["leaked_random_tree"]["accuracy"],
            "leaked_over_grouped": task["leakage_gap"]["leaked_over_grouped"],
            "majority_rate": task["baselines"]["global_majority_rate"],
            "type_level_majority_rate":
                task["baselines"]["type_level_majority_rate"],
            "clade_oracle_accuracy":
                task["baselines"]["clade_oracle"]["accuracy"],
            "clade_oracle_type_level":
                task["baselines"]["clade_oracle"]["type_level_accuracy"],
            "null_mean": null.get("mean"),
            "null_sd": null.get("sd"),
            "null_type_level_mean": null.get("type_level_mean"),
            "nn_grouped": task.get("nearest_neighbour_leakage_demo", {}).get(
                "grouped_accuracy"),
            "nn_leaked_random": task.get(
                "nearest_neighbour_leakage_demo", {}).get(
                "leaked_random_accuracy"),
            "verdict": task["verdict"]}
    summary["note"] = (
        "the leaked_over_grouped column is the point of the whole file: it is "
        "above 1.3 wherever the label is a property of the type and essentially "
        "1.0 for the one label that varies within a type (the beta subunit). "
        "That contrast is the leakage mechanism caught in the act -- a random "
        "split helps only when the test entry has a near-identical twin of the "
        "same type in training")
    results["summary"] = summary

    results["method"] = {
        "why_this_module_exists": (
            "the data set is large enough to invite a sequence-to-reaction "
            "classifier, and a careless one would report an accuracy close to "
            "1.0 that means nothing"),
        "leakage": (
            "entries are not independent. 11,422 entries hold 10,019 distinct "
            "sequences, and the reaction class and the substrate class are "
            "curated per TYPE, not per entry. Under a random split a model "
            "recognises the type of a test sequence from its near-identical "
            "training twins and reads the label off it"),
        "how_it_was_avoided": (
            "every headline number comes from grouped cross-validation with the "
            f"group being {GROUP_FIELD}: GroupKFold with {args.group_folds} "
            "entry-balanced folds plus leave-one-type-out. The random-split "
            "number sits next to it and is labelled invalid, because the size "
            "of the gap is itself a result"),
        "baselines": (
            "four of them: the global majority class, the type-level majority "
            "class, the per-fold training majority, and a clade oracle that "
            "answers with the majority label of the training types sharing the "
            "HMM group of the test entry. The clade oracle is the demanding "
            "one, because it asks whether residue identities add anything to "
            "deep clade membership"),
        "weights": (
            "each entry carries weight 1/(size of its type) so that a type with "
            "1,345 members and a type with 2 members have the same say; the "
            "unweighted variant is reported as a sensitivity"),
        "models": (
            "a depth-limited CART and an L1-penalised multinomial logistic "
            "regression, both written here with numpy because scikit-learn is "
            "not installed in this project environment. A random forest is "
            "fitted ONLY to rank features and never used for an accuracy"),
        "hyperparameters": (
            "tree depth was chosen inside each training fold by a nested "
            "grouped cross-validation, so the reported accuracy is not inflated "
            "by tuning on held-out types. The full depth sweep is also given, "
            "but its best row is not an achieved accuracy"),
        "effective_sample_size": (
            "61 types, not 11,422 entries. For a 7-class problem this is small: "
            "no claim about the small classes can be supported"),
        "primary_metric": (
            "two numbers everywhere: accuracy over entries, and type-level "
            "accuracy where each type is scored once by the majority vote of "
            "its members. The second is the honest one"),
        "seed": SEED,
    }
    results["cannot_be_supported"] = [
        "no claim about a reaction class represented by 2 or 3 types "
        "(C_N_cleavage, O_demethylation, angular_dioxygenation): under the "
        "grouped scheme those classes have one or two test folds in total",
        "no claim that the columns reported here are substrate-binding pocket "
        "residues. There is no structure in this database; what was measured is "
        "which alignment columns separate the classes",
        "no claim that a residue identity explains the chemistry rather than "
        "marking the clade. Reaction class and substrate class are both strongly "
        "structured by HMM group, so a sequence model can reach its accuracy by "
        "recognising the clade alone; the clade-oracle baseline is there "
        "precisely to make that visible",
        "no causal claim in either direction for question 3. Operon composition "
        "and substrate class are both properties of a lineage with a history, "
        "and a shared genome is not evidence of a shared pathway",
        "no claim that the curated labels are correct. Most entries sit in the "
        "distant or novel evidence tiers, where the reaction is unknown and was "
        "inherited from the type. This module treats the label as given",
        "no usable sequence-to-substrate predictor for a new protein. The "
        "grouped accuracies here are far from a decision rule, and the "
        "reference-pair calibration already showed that no global identity "
        "threshold guarantees the same substrate",
    ]
    results["inputs"] = {
        "db": os.path.abspath(args.db),
        "alignment": os.path.abspath(args.alignment),
        "chemistry": os.path.abspath(args.chemistry),
        "ecology": os.path.abspath(args.ecology),
        "sdp_positions": os.path.abspath(args.sdp),
        "motif_stats": motif_source,
        "entries": len(entries), "types": len(group_names),
        "hmm_groups": len(clade_names),
        "alignment_match_columns": int(matrix.shape[1]),
        "sequence_features": int(x_sequence.shape[1]),
        "context_features": int(no_taxonomy.sum()),
        "taxonomy_features": int(only_taxonomy.sum()),
    }
    results["runtime_seconds"] = round(time.time() - started, 1)

    os.makedirs(args.out_dir, exist_ok=True)
    path = os.path.join(args.out_dir, "learning.json")
    with open(path, "w", encoding="utf-8") as handle:
        json.dump(results, handle, indent=1)

    report(results)
    print(f"\n[yazildi] {path}")
    connection.close()


def report(results):
    """stdout ozeti -- gruplanmis, sizintili, taban cizgisi, null, kolonlar."""
    line = "=" * 78
    print("\n" + line)
    print("VERIDEN OGRENME -- GRUPLANMIS (TIP BAZLI) CAPRAZ DOGRULAMA")
    print(line)

    def block(key, title):
        task = results.get(key)
        if not task:
            return
        grouped = task["grouped_tree"]
        base = task["baselines"]
        print(f"\n{title}")
        print(f"  sinif {len(task['classes'])} | tip {task['n_types']} | "
              f"giris {task['n_entries']} | ozellik {task['n_features']}")
        rows = [("gruplanmis dogruluk (tip bazli CV)", grouped["accuracy"]),
                ("  tip duzeyi dogruluk", grouped["type_level_accuracy"]),
                ("  dengeli dogruluk (makro duyarlilik)",
                 grouped["balanced_accuracy"])]
        if "nested_depth_selection" in task:
            rows.append(("ic ice derinlik secimi (durust sayi)",
                         task["nested_depth_selection"]["accuracy"]))
        if "leave_one_type_out_tree" in task:
            rows.append(("her tipi disarida birak",
                         task["leave_one_type_out_tree"]["accuracy"]))
        if "grouped_l1_logistic" in task:
            rows.append(("L1 lojistik, gruplanmis",
                         task["grouped_l1_logistic"]["accuracy"]))
        rows += [("SIZINTILI rastgele bolme (GECERSIZ)",
                  task["leaked_random_tree"]["accuracy"]),
                 ("taban: global cogunluk sinifi", base["global_majority_rate"]),
                 ("taban: tip duzeyi cogunluk", base["type_level_majority_rate"]),
                 ("taban: kat basina egitim cogunlugu",
                  base["grouped_train_majority_accuracy"]),
                 ("taban: KLAD KAHINI (HMM grubu)",
                  base["clade_oracle"]["accuracy"]),
                 ("  klad kahini, tip duzeyi",
                  base["clade_oracle"]["type_level_accuracy"])]
        for name, value in rows:
            print(f"  {name:48} {value:.4f}")
        print(f"  {'sizinti ucurumu (kat)':48} "
              f"{task['leakage_gap']['leaked_over_grouped']:.2f}x")
        demo = task.get("nearest_neighbour_leakage_demo")
        if demo:
            print(f"  {'1-NN gruplanmis':48} {demo['grouped_accuracy']:.4f}")
            print(f"  {'1-NN SIZINTILI rastgele (GECERSIZ)':48} "
                  f"{demo['leaked_random_accuracy']:.4f}  "
                  f"({demo['leaked_over_grouped']:.2f}x)")
        null = task.get("permutation_null")
        if null:
            print(f"  {'null (etiket karistirma)':48} "
                  f"{null['mean']:.4f} +- {null['sd']:.4f}")
            print(f"    p95 {null['p95']:.4f}  azami {null['max']:.4f}  "
                  f"{null['reps']} tekrar  [{null['shuffled_at'][:28]}]")
            compare = task.get("grouped_vs_null", {})
            print(f"    tip duzeyi null ortalamasi "
                  f"{compare.get('type_level_null_mean')}")
            print(f"    z = {compare.get('z_against_null')}   "
                  f"ampirik p = {compare.get('empirical_p')}   "
                  f"tip duzeyi p = {compare.get('type_level_empirical_p')}")
        for part in str(task.get("verdict")).split("; "):
            print(f"  KARAR: {part}")
        print("  sinif basina:")
        for name, stats in sorted(grouped["per_class"].items(),
                                  key=lambda kv: -kv[1]["support"]):
            print(f"    {name:30} n={stats['support']:6}  "
                  f"recall={stats['recall']:.3f}  "
                  f"precision={stats['precision']:.3f}")
        print("  karisiklik matrisi (satir = gercek sinif):")
        width = max(len(c) for c in task["classes"]) + 1
        print("    " + " " * width
              + " ".join(f"{c[:7]:>7}" for c in task["classes"]))
        for name, row in zip(task["classes"], grouped["confusion_matrix"]):
            print(f"    {name:{width}}" + " ".join(f"{v:>7}" for v in row))

    block("q1_reaction_class_from_sequence",
          "SORU 1 -- reaksiyon sinifi <- dizi (hizalama kolonlari)")
    block("q3_substrate_class_from_context",
          "SORU 3 -- substrat sinifi <- genomik baglam")
    block("q4a_positive_control_ro_group_from_sequence",
          "SORU 4a -- POZITIF KONTROL: HMM grubu <- dizi")
    block("q4c_beta_subunit_from_sequence",
          "SORU 4c -- beta alt birimi <- dizi (tip ici degisen etiket)")

    task = results.get("q1_reaction_class_from_sequence", {})
    sweep = task.get("depth_sweep")
    if sweep:
        print("\nDERINLIK TARAMASI -- sizinti ucurumu kapasiteyle genisliyor")
        print(f"  {'derinlik':>9} {'gruplanmis':>12} {'sizintili':>12} {'kat':>7}")
        for row in sweep["rows"]:
            print(f"  {row['max_depth']:>9} {row['grouped_accuracy']:>12.4f} "
                  f"{row['leaked_random_accuracy']:>12.4f} {row['ratio']:>7.2f}")
    nested = task.get("nested_depth_selection")
    if nested:
        print(f"  ic ice secilen derinlikler: {nested['depths_chosen_per_fold']}")

    variants = results.get("q1_column_selection_variants")
    if variants:
        print("\nKOLON SECIMI -- hazir SDP listesi ne kadar sisiriyor")
        for row in variants["rows"]:
            print(f"  {row['selection']:44} {row['n_columns']:>4} kolon  "
                  f"{row['grouped_accuracy']:.4f}  "
                  f"secim sizintisi: {row['selection_leak'].split(' --')[0]}")

    q2 = results.get("q2_informative_columns")
    if q2:
        print("\nSORU 2 -- en bilgilendirici hizalama kolonlari")
        print(f"  {'kolon':>6} {'kalinti':>8} {'bolge':>17} "
              f"{'merkeze uzaklik':>18} {'agac kazanci':>13}")
        for row in q2["top_features"][:15]:
            print(f"  {row['column']:>6} {row['residue']:>8} {row['zone']:>17} "
                  f"{row['offset_from_nearest']:>+11} "
                  f"({row['nearest_catalytic_centre_column']:>3}) "
                  f"{row['tree_gain']:>12.5f}")
        print(f"  ilk {len(q2['top_features'])} ozelligin bolge dagilimi: "
              f"{q2['zone_counts_in_top_features']}")
        print(f"  havuzun bolge dagilimi: {q2['zone_counts_in_whole_pool']}")
        print(f"  katalitik bolge payi: ilk ozelliklerde "
              f"{q2['catalytic_domain_share_in_top_features']}, havuzda "
              f"{q2['catalytic_domain_share_in_whole_pool']} -> zenginlesme "
              f"{q2['catalytic_domain_enrichment_over_pool']}x")
        print(f"  bosluk durumu olan ozellik: "
              f"{q2['gap_state_features_in_top']}/{len(q2['top_features'])}")
        print(f"  katalitik merkez kolonuna +-10 icinde: "
              f"{q2['top_features_within_10_columns_of_a_centre_column']}"
              f"/{len(q2['top_features'])}")
        print(f"  L1 sifir olmayan katsayi: {q2['l1_nonzero_features']}"
              f"/{q2['l1_total_features']}")
        univariate = q2.get("univariate_type_level")
        if univariate:
            print("  tip duzeyi tek degiskenli en guclu kolonlar:")
            for row in univariate["rows"][:8]:
                print(f"    kolon {row['column']:>4}  "
                      f"V={row['cramers_v_type_level']:.3f}  {row['zone']:17}  "
                      f"merkeze {row['offset_from_nearest']:>+5} "
                      f"({row['nearest_catalytic_centre_column']})")

    ablations = results.get("q3_ablations")
    if ablations:
        print("\nSORU 3 ABLASYONLARI -- baglam mi, taksonomi mi")
        for name in ("context_without_taxonomy", "taxonomy_only",
                     "context_plus_taxonomy"):
            row = ablations.get(name)
            if row:
                print(f"  {name:28} {row['n_features']:>4} ozellik  "
                      f"gruplanmis {row['grouped_accuracy']:.4f}  "
                      f"dengeli {row['grouped_balanced_accuracy']:.4f}  "
                      f"tip {row['type_level_accuracy']:.4f}")

    modality = results.get("q4b_modality_comparison")
    if modality:
        print("\nSORU 4b -- bilgi nerede")
        print("  (giris duzeyi | tip duzeyi) ciftleri")
        print(f"  {'baslik':32} {'gruplanmis':>17} {'cogunluk':>17} "
              f"{'klad kahini':>17} {'sizintili':>10}")
        for name, row in modality.items():
            if name == "note":
                continue
            print(f"  {name:32} "
                  f"{row['grouped_accuracy']:>8.3f}|{row['type_level_accuracy']:<8.3f} "
                  f"{row['global_majority_rate']:>8.3f}|"
                  f"{row['type_level_majority_rate']:<8.3f} "
                  f"{row['clade_oracle_accuracy']:>8.3f}|"
                  f"{row['clade_oracle_type_level']:<8.3f} "
                  f"{row['leaked_random_accuracy']:>10.3f}")

    context_importance = results.get("q3_informative_context_features")
    if context_importance:
        print("\nSORU 3 -- en bilgilendirici baglam ozellikleri")
        for row in context_importance["top_features"][:12]:
            print(f"  {row['feature']:46} agac {row['tree_gain']:.5f}  "
                  f"L1 {row['l1_abs_weight']:.4f}")

    print("\nVERININ DESTEKLEMEDIGI SONUCLAR")
    for item in results.get("cannot_be_supported", []):
        print(f"  - {item}")


if __name__ == "__main__":
    sys.exit(main() or 0)
