"""
Veriden ogrenme -- "iyi bir veri setimiz var, daha iyi korelasyonlar bulabiliriz"
sorusunu DENETLENEBILIR bir denetimli ogrenme problemine cevirir.

NEDEN BU MODUL VAR
  Veritabaninda 11.422 dogrulanmis alpha alt birimi var ve her biri bir tipe,
  her tip de kuratorlu bir reaksiyon sinifina ve substrat sinifina bagli. Bu
  haliyle akla gelen ilk sey "dizi -> reaksiyon sinifi" siniflandiricisi
  kurmaktir. Ama bu veri setinde o is, dogru yapilmazsa, sonucu BASTAN
  BELLI bir sayi uretir.

SIZINTI (LEAKAGE) PROBLEMI -- bu modulun asil konusu
  Girisler BAGIMSIZ DEGIL. Iki ayri sebepten:
    1. Fazlalik: 11.422 giris yalnizca 10.019 tekil dizi (bkz. redundancy.json).
       Ayni protein birbirine cok yakin suslardan tekrar tekrar geliyor.
    2. Etiket tipin OZELLIGI: reaksiyon sinifi ve substrat sinifi tek tek
       girislere degil `ro.ro_cluster` tipine atanmis. Yani bir tipin bir
       uyesi egitimde, baska bir uyesi testte oldugunda model dizi-reaksiyon
       iliskisini DEGIL, "bu dizi hangi tipe benziyor" sorusunu cozuyor ve
       cevabi egitim kumesinden okuyor.
  Sonuc: rastgele bolunmus bir train/test ayrimi 1,0'a yakin ve ANLAMSIZ bir
  dogruluk verir. Bu sayi olculdu ve raporlaniyor -- ama SIZINTILI diye
  etiketlenerek, cunku onu gizlemek yerine gostermek ogreticidir.

NASIL ONLENDI
  * Butun basliklar GRUPLANMIS capraz dogrulama ile olculuyor; grup = enzim
    TIPI (`ro.ro_cluster`). Ayni tip hem egitimde hem testte asla bulunmuyor.
    Iki sema kosuluyor: GroupKFold (tip bazli, boyut dengeli 10 kat) ve
    her tipi sirayla disarida birakma (leave-one-type-out, 61 kat).
  * Ayni model ve ayni hiperparametreler rastgele bolunmus semayla da
    kosuluyor; iki sayi arasindaki UCURUM bir sonuc olarak veriliyor.
  * Her baslikta uc taban cizgisi var: (a) global cogunluk sinifi orani,
    (b) kat basina egitim cogunlugunu tahmin eden taban cizgisi, (c) etiket
    permutasyon null'i -- etiketler TIP duzeyinde karistirilip ayni gruplanmis
    sema tekrar kosulur, yuzlerce tekrar, ve null dagiliminin ortalamasi,
    standart sapmasi, %95'lik dilimi ve azamisi yazilir. Permutasyon tip
    duzeyinde yapilir, cunku etiket tipin ozelligi; giris duzeyinde
    karistirmak null'i yapay olarak dusururdu.
  * Model ACIK: derinligi sinirli karar agaci (bolme kurallari JSON'a
    yaziliyor) ve L1 cezali cok sinifli lojistik regresyon. Kapali bir
    topluluk (rastgele orman) YALNIZCA degisken onemi icin kuruluyor ve
    ciktida boyle isaretli.
  * Her girisin agirligi 1/(tipinin buyuklugu): 1.345 uyeli bir tip ile 2
    uyeli bir tip egitimde esit soz sahibi olsun diye. Agirliksiz surum de
    duyarlilik olarak kosuluyor.
  * Dizi ozellikleri icin kolon secimi UC sekilde yapiliyor: (1) hazir
    `analysis_out/sdp_positions.csv` ilk N kolonu -- bu dosya kolonlari TUM
    tipleri gorerek sectigi icin bir SECIM SIZINTISI tasir ve boyle
    etiketlenir; (2) ayni olcut her katin YALNIZCA egitim tipleriyle yeniden
    hesaplanarak (ic-kat secim, sizintisiz); (3) hic secim yapmadan dolu
    kolonlarin tamami. Ucunun gruplanmis dogrulugu yan yana veriliyor, yani
    hazir listeyi kullanmanin bedeli sayiyla gosteriliyor.

NE OLCULDU
  1. Reaksiyon sinifi (7 sinif, chemistry.csv) dizinin kendisinden -- hizalama
     kolonlarindaki kalinti kimlikleri (genomic_context/cand_aln.sto match
     kolonlari, model uzunlugu 426).
  2. Hangi kolonlar sinyali tasiyor: agac bolmeleri, L1 katsayilari, orman
     onemleri ve TIP duzeyinde tek degiskenli Cramer's V. Her kolon
     `analysis_out/motif_stats.json`'daki sekiz tanimlayici kolona gore
     isaretli uzaklikla ve bolge adiyla (rieske / baglanti / katalitik)
     raporlanir, boylece "aktif bolgede mi" sorusu biyolog tarafindan
     sorulabilir.
  3. Substrat sinifi (3 sinif, cluster_ecology.csv) GENOMIK BAGLAMDAN: operon
     bilesimi, elektron transfer ortak tipleri, duzenleyici ailesi, plazmit
     durumu ve taksonomi. Taksonomi AYRI bir blok tutulur ve ablasyonla
     kosulur, cunku taksonomi ile baglam ic ice gecmistir: "komsuluk kimyayi
     ongoruyor" iddiasi ancak taksonomi cikarildiginda sinanabilir.
  4. Ek olarak: (a) POZITIF KONTROL -- `ro.ro_group` (HMM grubu) diziden
     ongorulur; bu etiket dizi benzerliginden turedigi icin yontemin sinyali
     YAKALAYABILDIGINI gosterir, yani dusuk dogruluklar yontemin korlugu
     degil; (b) modalite karsilastirmasi -- ayni etiket ayni katlarla bir
     kez diziden bir kez baglamdan ongorulur, boylece bilginin NEREDE
     oldugu sorulur; (c) tipin ozelligi OLMAYAN bir etiket: beta alt
     biriminin varligi (alpha3 - alpha3beta3 mimarisi) -- tip icinde degistigi
     icin gruplanmis semada anlamli olcude ogrenilebilir tek baslik.

NE OLCULMUYOR / VERININ DESTEKLEMEDIGI SONUCLAR
  * Etkin ornek buyuklugu 11.422 degil 61'dir (tip sayisi). Yedi sinifli bir
    problemde 61 bagimsiz birim cok azdir; guven araliklari genistir ve
    kucuk siniflar (C_N_cleavage 2 tip, O_demethylation 3 tip,
    angular_dioxygenation 3 tip) icin hicbir genelleme iddiasi kurulamaz.
  * Ayirt edici bulunan kolonlar YAPISAL olarak dogrulanmis substrat baglama
    cebi kalintilari DEGILDIR. Veritabaninda yapi yok; olculen sey
    "siniflari ayiran hizalama kolonlari". Motif kolonlarina uzaklik bir
    HIPOTEZ uretir, kanit uretmez.
  * Substrat ve reaksiyon etiketleri tip duzeyinde kuratorludur; kanit
    kademesi `distant`/`novel` olan uyelerin (%58,2) gercek kimyasi
    bilinmiyor. Bu modul etiketi dogru varsayar; dogrulugunu sinamaz.
  * Baglamdan substrat sinifi ongorusu bir NEDENSELLIK iddiasi degildir.
    Ayni genomda bulunmak ortak yol kaniti degil; ustelik dizilenmis
    genomlar orneklem yanlilidir (PAH yikan izolatlar fazla temsil edilir).

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
PRIMARY_DEPTH = 5         # derinlik taramasindan sonra sabitlenen birincil derinlik

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
    """Chi-kare tabanli Cramer's V; bos satir/sutun atilir."""
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


# --- Karar agaci (numpy) --------------------------------------------------
#
# NEDEN KENDI UYGULAMASI. Projenin python ortaminda scikit-learn KURULU DEGIL
# (`python3 -c "import sklearn"` basarisiz). Bagimlilik eklemek yerine ihtiyac
# duyulan iki model numpy ile yazildi; boylece modul pipeline'in geri kalaniyla
# ayni ortamda kosuyor ve sayilar her makinede ayni cikiyor. scikit-learn
# kurulu BIR ortam varsa `--cross-check` ile ayni gruplanmis sema orada da
# kosulur ve iki uygulamanin dogrulugu JSON'a yan yana yazilir.
#
# Ozellikler ikili (0/1) oldugu icin CART cok basitlesiyor: her dugumde tek bir
# BLAS carpimi butun ozelliklerin sinif sayimlarini birden veriyor.

def _node_counts(x_sub, weighted_y, ones):
    """(sinif x ozellik) agirlikli sayim + agirlik toplami + agirliksiz sayim."""
    stacked = np.concatenate([weighted_y, ones], axis=1)
    return stacked.T @ x_sub


def fit_tree(x, y_onehot, weight, max_depth, min_leaf=TREE_MIN_LEAF,
             feature_mask=None):
    """Gini olcutlu, ikili ozellikli CART. Donen: ic ice sozluk.

    feature_mask verilirse yalnizca o ozellikler bolme adayi olur (rastgele
    orman icin gerekli).
    """
    n_features = x.shape[1]
    n_classes = y_onehot.shape[1]

    def build(index, depth):
        y_sub = y_onehot[index] * weight[index, None]
        totals = y_sub.sum(axis=0)
        node = {"n": int(index.size),
                "class_weight": totals,
                "prediction": int(np.argmax(totals))}
        if depth >= max_depth or index.size < 2 * min_leaf:
            return node
        if (totals > 0).sum() <= 1:
            return node
        x_sub = x[index]
        ones = np.ones((index.size, 1), dtype=np.float32)
        counts = _node_counts(x_sub, y_sub, ones)
        class_one = counts[:n_classes]            # agirlikli, ozellik==1
        plain_one = counts[n_classes]             # agirliksiz adet, ozellik==1
        weight_one = class_one.sum(axis=0)
        weight_total = totals.sum()
        weight_zero = weight_total - weight_one
        class_zero = totals[:, None] - class_one
        safe_one = np.where(weight_one > 0, weight_one, 1.0)
        safe_zero = np.where(weight_zero > 0, weight_zero, 1.0)
        gini_one = 1.0 - ((class_one / safe_one) ** 2).sum(axis=0)
        gini_zero = 1.0 - ((class_zero / safe_zero) ** 2).sum(axis=0)
        impurity = (weight_one * gini_one + weight_zero * gini_zero) / weight_total
        # Bir yana min_leaf'ten az GIRIS dusuyorsa o bolme yasak. Olcut
        # agirlik degil adet, cunku agirliklar 1/tip_boyutu ile cok kucuk
        # olabiliyor ve bir agirlik esigi buyuk tipleri kayiriyordu.
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

    assert n_features > 0
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


def tree_rules(node, names, class_names, depth=0, max_depth=3):
    """Agaci okunabilir satirlara cevir -- JSON'da model boyle denetlenir."""
    lines = []

    def walk(node, prefix, level):
        total = float(node["class_weight"].sum())
        share = (float(node["class_weight"][node["prediction"]]) / total
                 if total > 0 else 0.0)
        if "feature" not in node or level >= max_depth:
            lines.append({"path": prefix or "root", "n": node["n"],
                          "predicted_class": class_names[node["prediction"]],
                          "class_share": round(share, 3)})
            return
        column, residue = names[node["feature"]]
        test = f"{column}={residue}"
        walk(node["right"], (prefix + " AND " if prefix else "") + test, level + 1)
        walk(node["left"],
             (prefix + " AND " if prefix else "") + "NOT " + test, level + 1)

    walk(node, "", 0)
    return lines


# --- L1 cezali cok sinifli lojistik regresyon (numpy, FISTA) --------------
#
# Yorumlanabilirlik icin L1: katsayilarin cogu tam sifir kaliyor, yani model
# "hangi kolonlar" sorusuna dogrudan cevap veriyor. Cozucu yakinsama hizlandirici
# ile proksimal gradyan (FISTA); ozellikler ikili oldugu icin yalnizca egitim
# ortalamasi cikariliyor, ayrica olcekleme gerekmiyor.

def fit_l1_logistic(x, y_index, weight, n_classes, penalty=0.002,
                    iterations=250):
    n, n_features = x.shape
    mean = (x * weight[:, None]).sum(axis=0) / max(weight.sum(), 1e-12)
    centered = x - mean
    target = np.zeros((n, n_classes), dtype=np.float32)
    target[np.arange(n), y_index] = 1.0
    norm = weight / max(weight.sum(), 1e-12)
    # Lipschitz sabiti: guc iterasyonuyla en buyuk tekil degerin karesi.
    vector = np.random.RandomState(0).randn(n_features).astype(np.float32)
    for _ in range(12):
        projected = centered @ vector
        vector = centered.T @ (projected * norm)
        length = np.linalg.norm(vector)
        if length < 1e-20:
            break
        vector /= length
    lipschitz = max(float(length) * 0.5, 1e-6)
    step = 1.0 / lipschitz

    coef = np.zeros((n_features, n_classes), dtype=np.float32)
    bias = np.zeros(n_classes, dtype=np.float32)
    momentum = coef.copy()
    previous = coef.copy()
    theta = 1.0
    for iteration in range(iterations):
        scores = centered @ momentum + bias
        scores -= scores.max(axis=1, keepdims=True)
        np.exp(scores, out=scores)
        scores /= scores.sum(axis=1, keepdims=True)
        residual = (scores - target) * norm[:, None]
        gradient = centered.T @ residual
        bias -= step * residual.sum(axis=0) * max(lipschitz, 1.0)
        candidate = momentum - step * gradient
        # proksimal adim: yumusak esikleme
        coef = np.sign(candidate) * np.maximum(np.abs(candidate) - step * penalty, 0.0)
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

def group_kfold(groups, n_folds):
    """Tip bazli kat ayrimi; katlar GIRIS sayisina gore dengelenir.

    scikit-learn'un GroupKFold'u ile ayni acgozlu kural: tipler buyukten
    kucuge siralanir ve her biri o an en az yuklu kata atanir. Bir tip tek bir
    kata gider, yani ayni tip hem egitimde hem testte ASLA bulunmaz.
    """
    sizes = Counter(groups)
    order = sorted(sizes, key=lambda g: (-sizes[g], g))
    load = [0] * n_folds
    assignment = {}
    for group in order:
        target = int(np.argmin(load))
        assignment[group] = target
        load[target] += sizes[group]
    fold_of = np.array([assignment[g] for g in groups])
    return [(np.where(fold_of != k)[0], np.where(fold_of == k)[0])
            for k in range(n_folds)]


def leave_one_group_out(groups):
    unique = sorted(set(groups))
    arr = np.asarray(groups)
    return [(np.where(arr != g)[0], np.where(arr == g)[0]) for g in unique]


def random_kfold(n, n_folds, rng):
    """SIZINTILI sema: giris duzeyinde rastgele bolme, tip yapisi yoksayilir."""
    order = rng.permutation(n)
    return [(np.setdiff1d(order, order[k::n_folds]), order[k::n_folds])
            for k in range(n_folds)]


# --- Olcumler -------------------------------------------------------------

def score(predicted, actual, groups, class_names):
    n_classes = len(class_names)
    accuracy = float((predicted == actual).mean())
    confusion = np.zeros((n_classes, n_classes), dtype=int)
    for a, p in zip(actual, predicted):
        confusion[a, p] += 1
    per_class = {}
    recalls = []
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
    type_hits, type_total = 0, 0
    for group in sorted(set(groups)):
        members = np.asarray(groups) == group
        vote = np.bincount(predicted[members], minlength=n_classes).argmax()
        type_total += 1
        type_hits += int(vote == actual[members][0])
    return {"accuracy": round(accuracy, 4),
            "balanced_accuracy": round(float(np.mean(recalls)) if recalls else 0.0, 4),
            "type_level_accuracy": round(type_hits / type_total, 4),
            "n_types": type_total,
            "n_entries": int(len(actual)),
            "per_class": per_class,
            "confusion_matrix": confusion.tolist(),
            "confusion_rows_are_true_classes": True}


def majority_baselines(folds, y_index, groups, class_names):
    """Iki taban cizgisi: global cogunluk ve kat basina egitim cogunlugu."""
    n_classes = len(class_names)
    global_majority = int(np.bincount(y_index, minlength=n_classes).argmax())
    global_rate = float((y_index == global_majority).mean())
    predicted = np.empty_like(y_index)
    for train, test in folds:
        majority = int(np.bincount(y_index[train], minlength=n_classes).argmax())
        predicted[test] = majority
    per_fold = score(predicted, y_index, groups, class_names)
    type_counts = Counter()
    seen = set()
    for g, k in zip(groups, y_index):
        if g not in seen:
            seen.add(g)
            type_counts[k] += 1
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


def permutation_null(x, folds, groups, y_index, n_classes, weight, depth,
                     reps, rng, log_every=25):
    """Etiketleri TIP duzeyinde karistir, ayni gruplanmis semayi tekrar kos.

    Giris duzeyinde karistirmak null'i yapay olarak dusururdu: etiket tipin
    ozelligi oldugu icin asil soru "tip-etiket eslesmesi rastgele olsaydi bu
    sema ne kadar dogruluk uretirdi" sorusudur.
    """
    unique = sorted(set(groups))
    index_of = {g: i for i, g in enumerate(unique)}
    group_index = np.array([index_of[g] for g in groups])
    group_label = np.array([y_index[group_index == i][0]
                            for i in range(len(unique))])
    accuracies, type_accuracies = [], []
    for rep in range(reps):
        shuffled = rng.permutation(group_label)
        fake = shuffled[group_index]
        onehot = np.zeros((len(fake), n_classes), dtype=np.float32)
        onehot[np.arange(len(fake)), fake] = 1.0
        predicted = np.empty(len(fake), dtype=np.int32)
        for train, test in folds:
            node = fit_tree(x[train], onehot[train], weight[train], depth)
            predicted[test] = predict_tree(node, x[test])
        accuracies.append(float((predicted == fake).mean()))
        hits = 0
        for i in range(len(unique)):
            members = group_index == i
            vote = np.bincount(predicted[members], minlength=n_classes).argmax()
            hits += int(vote == shuffled[i])
        type_accuracies.append(hits / len(unique))
        if log_every and (rep + 1) % log_every == 0:
            print(f"      null {rep + 1}/{reps} ...", flush=True)
    accuracies = np.array(accuracies)
    type_accuracies = np.array(type_accuracies)
    return {"reps": reps,
            "shuffled_at": "type level (the label is a property of the type)",
            "mean": round(float(accuracies.mean()), 4),
            "sd": round(float(accuracies.std(ddof=1)) if reps > 1 else 0.0, 4),
            "p95": round(float(np.quantile(accuracies, 0.95)), 4),
            "max": round(float(accuracies.max()), 4),
            "type_level_mean": round(float(type_accuracies.mean()), 4),
            "type_level_max": round(float(type_accuracies.max()), 4)}


def empirical_p(observed, null):
    """Tek yanli ampirik p: null'da gozlenenden iyi veya esit kac tekrar var."""
    if not null:
        return None
    better = sum(1 for value in null if value >= observed)
    return round((better + 1) / (len(null) + 1), 4)


# --- Veri yukleme ---------------------------------------------------------

def read_csv_dict(path, key="cluster"):
    with open(path, newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    return {row[key]: row for row in rows}


def load_entries(connection):
    """Dogrulanmis her RO icin kimlik, tip, grup ve replikon bilgisi."""
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
    length = len(aligned[candidate_ids[0]])
    joined = "".join(aligned[cid] for cid in candidate_ids)
    matrix = np.frombuffer(joined.encode("ascii"), dtype="S1")
    return matrix.reshape(len(candidate_ids), length)


def onehot_columns(matrix, columns, min_count=MIN_RESIDUE_COUNT):
    """(kolon, kalinti) ikili ozellik matrisi. Bosluk da bir durumdur."""
    features, names = [], []
    n = matrix.shape[0]
    for column in columns:
        values, counts = np.unique(matrix[:, column - 1], return_counts=True)
        for value, count in zip(values, counts):
            if count >= min_count and count <= n - min_count:
                features.append((matrix[:, column - 1] == value))
                names.append((int(column), value.decode("ascii")))
    if not features:
        raise SystemExit("[hata] hicbir kolon ozellik esigini gecmedi")
    return np.stack(features, axis=1).astype(np.float32), names


def column_occupancy(matrix):
    return (matrix != b"-").mean(axis=0)


def sdp_scores(matrix, groups, columns, min_members=5, min_filled=0.5):
    """analyze_variants.py ile AYNI olcut: kumeler arasi eksi kume ici entropi.

    Burada yeniden hesaplaniyor, cunku her katin yalnizca EGITIM tipleriyle
    secim yapmasi gerekiyor; hazir CSV tum tipleri gormus durumda.
    """
    by_group = defaultdict(list)
    for i, group in enumerate(groups):
        by_group[group].append(i)
    scores = {}
    total_groups = len(by_group)
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
        raw["operon_completeness"][i] = int(row[7] or 0)
        genes = int(row[8] or 0)
        raw["operon_size_class"][i] = ("1" if genes <= 1 else "2-3" if genes <= 3
                                       else "4-6" if genes <= 6 else "7+")

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
                                      else "<=100" if gap <= 100
                                      else "101-400" if gap <= 400
                                      else "401-1000" if gap <= 1000 else ">1000")

    # Operon icindeki gen kategorileri ve bilesenleri -- varlik/yokluk.
    categories = defaultdict(set)
    components = defaultdict(set)
    for candidate_id, component, category in connection.execute(
            "SELECT candidate_id, component, category FROM operon_gene"):
        if candidate_id not in position:
            continue
        if category and category != "ro_alpha":
            categories[candidate_id].add(category)
        if component and component not in ("none", "alpha"):
            components[candidate_id].add(component)
    category_vocabulary = sorted({c for s in categories.values() for c in s})
    component_vocabulary = sorted({c for s in components.values() for c in s})

    # +-10 kb penceresindeki bilesenler (operon disi komsuluk).
    window = defaultdict(set)
    for candidate_id, component in connection.execute("""
            SELECT n.candidate_id, c.component
            FROM neighbor n JOIN neighbor_protein p
              ON p.neighbor_id = n.neighbor_id
            JOIN neighbor_component c ON c.protein_key = p.protein_key
            WHERE c.component IS NOT NULL AND c.component <> 'none'"""):
        if candidate_id in position:
            window[candidate_id].add(component)
    window_vocabulary = sorted({c for s in window.values() for c in s})

    features, names, blocks = [], [], []

    def add(name, vector, block):
        features.append(np.asarray(vector, dtype=np.float32))
        names.append(name)
        blocks.append(block)

    for key in sorted(raw):
        values = raw[key]
        if all(isinstance(v, int) or v is None for v in values):
            add(f"{key}=1", [1.0 if v else 0.0 for v in values], "operon")
        else:
            for level in sorted({v for v in values if v is not None}):
                add(f"{key}={level}", [1.0 if v == level else 0.0 for v in values],
                    "regulation" if "regul" in key or "upstream" in key
                    or "intergenic" in key else "operon")
    for level in category_vocabulary:
        add(f"operon_gene_category={level}",
            [1.0 if level in categories.get(cid, ()) else 0.0 for cid in ids],
            "operon")
    for level in component_vocabulary:
        add(f"operon_component={level}",
            [1.0 if level in components.get(cid, ()) else 0.0 for cid in ids],
            "operon")
    for level in window_vocabulary:
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
        for level, count in counter.items():
            if count >= 50:
                add(f"lineage_rank{rank}={level}",
                    [1.0 if (len(l) > rank and l[rank] == level) else 0.0
                     for l in lineages], "taxonomy")
    genera = [e["organism"].split()[0] if e["organism"] else "unknown"
              for e in entries]
    for level, count in Counter(genera).items():
        if count >= 25:
            add(f"genus={level}", [1.0 if g == level else 0.0 for g in genera],
                "taxonomy")

    matrix = np.stack(features, axis=1)
    return matrix, names, blocks


# --- Bir basligin tam degerlendirmesi ------------------------------------

def evaluate_task(x, names, y_labels, groups, class_names, rng, args,
                  with_logo=False, with_null=True, depth=PRIMARY_DEPTH,
                  depth_sweep=False, label="task"):
    index_of = {name: i for i, name in enumerate(class_names)}
    y_index = np.array([index_of[v] for v in y_labels], dtype=np.int32)
    n_classes = len(class_names)
    onehot = np.zeros((len(y_index), n_classes), dtype=np.float32)
    onehot[np.arange(len(y_index)), y_index] = 1.0

    sizes = Counter(groups)
    weight_balanced = np.array([1.0 / sizes[g] for g in groups], dtype=np.float32)
    weight_plain = np.ones(len(y_index), dtype=np.float32)

    grouped = group_kfold(groups, args.group_folds)
    random_folds = random_kfold(len(y_index), N_RANDOM_FOLDS,
                                np.random.RandomState(SEED))

    def run_tree(folds, weight, tree_depth):
        predicted = np.empty(len(y_index), dtype=np.int32)
        for train, test in folds:
            node = fit_tree(x[train], onehot[train], weight[train], tree_depth)
            predicted[test] = predict_tree(node, x[test])
        return predicted

    result = {"n_entries": int(len(y_index)), "n_types": len(set(groups)),
              "n_features": int(x.shape[1]), "classes": list(class_names),
              "class_entry_counts": {c: int((y_index == i).sum())
                                     for i, c in enumerate(class_names)},
              "class_type_counts": dict(Counter(
                  y_labels[list(groups).index(g)] for g in sorted(set(groups))))}

    print(f"  [{label}] gruplanmis agac (derinlik {depth}) ...", flush=True)
    grouped_predicted = run_tree(grouped, weight_balanced, depth)
    result["grouped_tree"] = score(grouped_predicted, y_index, groups, class_names)
    result["grouped_tree"]["scheme"] = (
        f"GroupKFold, {args.group_folds} folds, group = {GROUP_FIELD}")

    print(f"  [{label}] SIZINTILI rastgele sema ...", flush=True)
    leaked_predicted = run_tree(random_folds, weight_balanced, depth)
    result["leaked_random_tree"] = score(leaked_predicted, y_index, groups,
                                         class_names)
    result["leaked_random_tree"]["scheme"] = (
        f"random {N_RANDOM_FOLDS}-fold split over entries; INVALID here because "
        f"near-identical sequences of the same type land on both sides and the "
        f"label is a property of the type")
    gap = (result["leaked_random_tree"]["accuracy"]
           - result["grouped_tree"]["accuracy"])
    result["leakage_gap"] = {
        "grouped_accuracy": result["grouped_tree"]["accuracy"],
        "leaked_accuracy": result["leaked_random_tree"]["accuracy"],
        "absolute_gap": round(gap, 4),
        "leaked_over_grouped": round(
            result["leaked_random_tree"]["accuracy"]
            / max(result["grouped_tree"]["accuracy"], 1e-9), 3)}

    result["baselines"] = majority_baselines(grouped, y_index, groups, class_names)

    # Agirliksiz duyarlilik: buyuk tipler egitimde baskin olunca ne degisir.
    plain_predicted = run_tree(grouped, weight_plain, depth)
    result["sensitivity_unweighted"] = {
        "accuracy": round(float((plain_predicted == y_index).mean()), 4),
        "note": ("entries weighted equally instead of 1/type-size, so the "
                 "largest types dominate training")}

    if depth_sweep:
        sweep = []
        for tree_depth in DEPTH_SWEEP:
            grouped_accuracy = float(
                (run_tree(grouped, weight_balanced, tree_depth) == y_index).mean())
            leaked_accuracy = float(
                (run_tree(random_folds, weight_balanced, tree_depth)
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
                     "buys you by memorising types: it widens with model "
                     "capacity, which is the signature of leakage rather than "
                     "of learning")}

    if with_logo:
        print(f"  [{label}] her tipi sirayla disarida birak ({len(set(groups))}) ...",
              flush=True)
        logo_predicted = run_tree(leave_one_group_out(groups), weight_balanced,
                                  depth)
        result["leave_one_type_out_tree"] = score(logo_predicted, y_index, groups,
                                                  class_names)
        result["leave_one_type_out_tree"]["scheme"] = (
            "leave-one-type-out: each type is predicted by a model that never "
            "saw a single member of it")

    print(f"  [{label}] L1 lojistik ...", flush=True)
    logistic_predicted = np.empty(len(y_index), dtype=np.int32)
    for train, test in grouped:
        model = fit_l1_logistic(x[train], y_index[train], weight_balanced[train],
                                n_classes, penalty=args.l1_penalty,
                                iterations=args.l1_iterations)
        logistic_predicted[test] = predict_l1_logistic(model, x[test])
    result["grouped_l1_logistic"] = score(logistic_predicted, y_index, groups,
                                          class_names)
    logistic_leaked = np.empty(len(y_index), dtype=np.int32)
    for train, test in random_folds:
        model = fit_l1_logistic(x[train], y_index[train], weight_balanced[train],
                                n_classes, penalty=args.l1_penalty,
                                iterations=args.l1_iterations)
        logistic_leaked[test] = predict_l1_logistic(model, x[test])
    result["leaked_random_l1_logistic"] = {
        "accuracy": round(float((logistic_leaked == y_index).mean()), 4)}

    if with_null and args.null_reps > 0:
        print(f"  [{label}] permutasyon null'i ({args.null_reps}) ...", flush=True)
        started = time.time()
        null = permutation_null(x, grouped, groups, y_index, n_classes,
                                weight_balanced, depth, args.null_reps, rng)
        null["seconds"] = round(time.time() - started, 1)
        result["permutation_null"] = null
        result["grouped_vs_null"] = {
            "grouped_accuracy": result["grouped_tree"]["accuracy"],
            "null_mean": null["mean"],
            "null_sd": null["sd"],
            "null_max": null["max"],
            "z_against_null": (round((result["grouped_tree"]["accuracy"]
                                      - null["mean"]) / null["sd"], 2)
                               if null["sd"] > 0 else None),
            "above_null_max": bool(result["grouped_tree"]["accuracy"] > null["max"]),
        }

    # Model icerigini goster: agac tamamen veriyle kurulur ve denetlenebilir.
    full_tree = fit_tree(x, onehot, weight_balanced, min(depth, 3))
    result["tree_rules_depth3_fitted_on_all_data"] = tree_rules(
        full_tree, names, class_names, max_depth=3)
    result["tree_rules_note"] = (
        "printed from a depth-3 tree fitted on every entry, so it is a "
        "description of the data and NOT a performance estimate")
    return result, y_index, onehot, weight_balanced, grouped


# --- Degisken onemi ve motif haritasi ------------------------------------

def importance_block(x, names, y_onehot, weight, grouped, rng, args,
                     motif_columns, matrix=None, groups=None, y_labels=None,
                     sequence=True, top=25):
    n_features = x.shape[1]

    # 1. Gruplanmis katlarda agac bolme kazanclari (kat disi, sema ile tutarli).
    gains = np.zeros(n_features)
    for train, _ in grouped:
        node = fit_tree(x[train], y_onehot[train], weight[train], PRIMARY_DEPTH)
        gains += tree_gains(node, n_features)
    gains /= max(len(grouped), 1)

    # 2. L1 katsayilari (tum veriye uydurulmus, betimleyici).
    n_classes = y_onehot.shape[1]
    y_index = np.argmax(y_onehot, axis=1)
    model = fit_l1_logistic(x, y_index, weight, n_classes,
                            penalty=args.l1_penalty,
                            iterations=args.l1_iterations)
    l1_weight = np.abs(model["coef"]).sum(axis=1)
    nonzero = int((l1_weight > 1e-8).sum())

    # 3. Rastgele orman -- SADECE onem icin, karar modeli olarak kullanilmiyor.
    forest = forest_importances(x, y_onehot, weight, args.forest_trees,
                                PRIMARY_DEPTH + 1, rng)

    def rank(values):
        order = np.argsort(-values)
        return {int(i): int(r) + 1 for r, i in enumerate(order)}

    rank_tree, rank_l1, rank_forest = rank(gains), rank(l1_weight), rank(forest)
    combined = np.array([rank_tree[i] + rank_l1[i] + rank_forest[i]
                         for i in range(n_features)])
    order = np.argsort(combined)

    rows = []
    for i in order[:top]:
        entry = {"feature": (f"column {names[i][0]} residue {names[i][1]}"
                             if sequence else names[i]),
                 "tree_gain": round(float(gains[i]), 5),
                 "l1_abs_weight": round(float(l1_weight[i]), 5),
                 "forest_gain": round(float(forest[i]), 5),
                 "rank_tree": rank_tree[i], "rank_l1": rank_l1[i],
                 "rank_forest": rank_forest[i]}
        if sequence:
            column = names[i][0]
            entry.update(motif_relation(column, motif_columns))
        rows.append(entry)

    block = {"top_features": rows,
             "l1_nonzero_features": nonzero,
             "l1_total_features": n_features,
             "method": ("features ranked by the sum of three ranks: mean split "
                        "gain of the depth-limited tree across the grouped "
                        "training folds, absolute L1 logistic weight summed over "
                        "classes, and mean split gain of a random forest. The "
                        "forest is fitted ONLY to rank features; every accuracy "
                        "reported in this file comes from the tree or the L1 "
                        "model, which can be read")}

    # 4. Kolon duzeyinde TIP bazli tek degiskenli olcut: her tip bir kez sayilir,
    #    yani buyuk tipler ne sinyali ne de gurultuyu sisirmiyor.
    if sequence and matrix is not None:
        per_column = []
        type_label = {}
        type_residue = defaultdict(Counter)
        for group, label, row in zip(groups, y_labels, range(matrix.shape[0])):
            type_label[group] = label
        unique_columns = sorted({names[i][0] for i in range(n_features)})
        class_names = sorted(set(y_labels))
        for column in unique_columns:
            type_residue.clear()
            values = matrix[:, column - 1]
            consensus = {}
            for group in sorted(set(groups)):
                members = np.asarray(groups) == group
                counter = Counter(v for v in values[members] if v != b"-")
                if counter:
                    consensus[group] = counter.most_common(1)[0][0].decode("ascii")
            residues = sorted(set(consensus.values()))
            table = [[sum(1 for g, r in consensus.items()
                          if r == residue and type_label[g] == label)
                      for label in class_names] for residue in residues]
            per_column.append({"column": column,
                               "n_types_scored": len(consensus),
                               "n_residue_states": len(residues),
                               "cramers_v_type_level": round(cramers_v(table), 3),
                               **motif_relation(column, motif_columns)})
        per_column.sort(key=lambda r: -r["cramers_v_type_level"])
        block["univariate_type_level"] = {
            "note": ("one row per alignment column. The residue of a TYPE is the "
                     "consensus of its members, so each of the types counts once "
                     "and a large sequenced clade cannot inflate the association. "
                     "Cramer's V on a 61-type table is noisy and is given for "
                     "ordering, not as a test"),
            "rows": per_column[:top]}
    return block


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
    sites.append({"column": model["bridging"][0], "label": "bridging carboxylate",
                  "expected": model["bridging"][1], "conserved_fraction": None})
    return sites, "ro_motif.py constants"


# --- Kolon secimi semalari -----------------------------------------------

def selection_variants(matrix, groups, y_labels, class_names, args, rng,
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
    sizes = Counter(groups)
    weight = np.array([1.0 / sizes[g] for g in groups], dtype=np.float32)
    grouped = group_kfold(groups, args.group_folds)
    occupied = [int(c) + 1 for c in np.where(occupancy >= args.min_occupancy)[0]]
    rows = []

    def grouped_accuracy(columns_per_fold, fixed_columns=None):
        predicted = np.empty(len(y_index), dtype=np.int32)
        for k, (train, test) in enumerate(grouped):
            columns = fixed_columns if fixed_columns else columns_per_fold[k]
            x, _ = onehot_columns(matrix, columns)
            node = fit_tree(x[train], onehot[train], weight[train], PRIMARY_DEPTH)
            predicted[test] = predict_tree(node, x[test])
        return float((predicted == y_index).mean())

    fixed = sorted(precomputed_columns[:args.top_sdp])
    rows.append({"selection": "precomputed_sdp_positions_csv",
                 "n_columns": len(fixed),
                 "grouped_accuracy": round(grouped_accuracy(None, fixed), 4),
                 "leak": ("yes -- the column list in sdp_positions.csv was "
                          "scored with every type visible, including the types "
                          "that are held out in each fold")})
    print(f"      selection precomputed: {rows[-1]['grouped_accuracy']}", flush=True)

    per_fold = []
    for train, _ in grouped:
        train_groups = [groups[i] for i in train]
        scores = sdp_scores(matrix[train], train_groups, occupied)
        ranked = sorted(scores, key=lambda c: -scores[c])[:args.top_sdp]
        per_fold.append(sorted(ranked))
    overlap = [len(set(per_fold[k]) & set(fixed)) / max(len(fixed), 1)
               for k in range(len(per_fold))]
    rows.append({"selection": "sdp_recomputed_inside_each_training_fold",
                 "n_columns": args.top_sdp,
                 "grouped_accuracy": round(grouped_accuracy(per_fold), 4),
                 "mean_overlap_with_precomputed_list": round(
                     float(np.mean(overlap)), 3),
                 "leak": "no -- the held-out types never take part in the scoring"})
    print(f"      selection in-fold: {rows[-1]['grouped_accuracy']}", flush=True)

    rows.append({"selection": "no_selection_all_occupied_columns",
                 "n_columns": len(occupied),
                 "grouped_accuracy": round(grouped_accuracy(None, occupied), 4),
                 "leak": ("no -- no label was consulted; occupancy is the only "
                          "filter")})
    print(f"      selection none: {rows[-1]['grouped_accuracy']}", flush=True)
    return {"rows": rows,
            "note": ("if the precomputed row were far above the in-fold row, "
                     "the headline accuracy would be partly an artefact of "
                     "choosing the columns on the full data set")}


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
    parser.add_argument("--l1-iterations", type=int, default=250)
    parser.add_argument("--skip-selection-variants", action="store_true")
    args = parser.parse_args()

    started = time.time()
    rng = np.random.RandomState(SEED)
    connection = sqlite3.connect(args.db)

    print("[okunuyor] etiketler ve girisler")
    chemistry = read_csv_dict(args.chemistry)
    ecology = read_csv_dict(args.ecology)
    entries = load_entries(connection)
    entries = [e for e in entries
               if e["cluster"] in chemistry and e["cluster"] in ecology]
    candidate_ids = [e["candidate_id"] for e in entries]
    groups = [e["cluster"] for e in entries]
    print(f"[bilgi] {len(entries)} giris, {len(set(groups))} tip")

    print("[okunuyor] hizalama")
    matrix = load_alignment_matrix(args.alignment, candidate_ids)
    occupancy = column_occupancy(matrix)
    print(f"[bilgi] {matrix.shape[1]} match kolonu, "
          f"{int((occupancy >= args.min_occupancy).sum())} tanesi "
          f">= {args.min_occupancy} dolu")

    motif_columns, motif_source = load_motif_columns(args.motif_stats)

    with open(args.sdp, newline="", encoding="utf-8") as handle:
        sdp_rows = list(csv.DictReader(handle))
    precomputed = [int(row["column"]) for row in sdp_rows]
    sequence_columns = sorted(precomputed[:args.top_sdp])
    x_sequence, sequence_names = onehot_columns(matrix, sequence_columns)
    print(f"[bilgi] dizi ozellikleri: {x_sequence.shape[1]} ikili ozellik, "
          f"{len(sequence_columns)} kolon")

    print("[kuruluyor] genomik baglam ozellikleri")
    x_context_all, context_names, context_blocks = build_context_features(
        connection, entries)
    blocks = np.array(context_blocks)
    no_taxonomy = blocks != "taxonomy"
    only_taxonomy = blocks == "taxonomy"
    print(f"[bilgi] baglam ozellikleri: {x_context_all.shape[1]} "
          f"({int(no_taxonomy.sum())} baglam + {int(only_taxonomy.sum())} taksonomi)")

    reaction = [chemistry[g]["reaction_class"] for g in groups]
    substrate = [ecology[g]["substrate_class"] for g in groups]
    ro_group = [e["group"] for e in entries]
    reaction_classes = sorted(set(reaction))
    substrate_classes = sorted(set(substrate))
    group_classes = sorted(set(ro_group))

    results = {}

    # --- Soru 1 + 2: reaksiyon sinifi diziden
    print("\n=== Soru 1: reaksiyon sinifi <- dizi")
    q1, y1, onehot1, weight1, grouped1 = evaluate_task(
        x_sequence, sequence_names, reaction, groups, reaction_classes, rng, args,
        with_logo=True, depth_sweep=True, label="reaction<-sequence")
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
                          rng, args, motif_columns, matrix=matrix, groups=groups,
                          y_labels=reaction, sequence=True)
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
    zones = Counter(row["zone"] for row in q2["top_features"])
    q2["zone_counts_in_top_features"] = dict(zones)
    near = sum(1 for row in q2["top_features"] if row["within_10_columns_of_centre"])
    q2["top_features_within_10_columns_of_a_centre_column"] = near
    results["q2_informative_columns"] = q2

    if not args.skip_selection_variants:
        print("\n=== Secim semalari (secim sizintisinin bedeli)")
        results["q1_column_selection_variants"] = selection_variants(
            matrix, groups, reaction, reaction_classes, args, rng,
            precomputed, occupancy)

    # --- Soru 3: substrat sinifi genomik baglamdan
    print("\n=== Soru 3: substrat sinifi <- genomik baglam")
    q3, y3, onehot3, weight3, grouped3 = evaluate_task(
        x_context_all[:, no_taxonomy], list(np.array(context_names)[no_taxonomy]),
        substrate, groups, substrate_classes, rng, args, with_logo=True,
        label="substrate<-context")
    q3["question"] = ("does the genomic neighbourhood of a gene predict the "
                      "chemical class of its substrate, for a type never seen "
                      "in training")
    q3["features"] = {
        "kind": "genomic context only, no residue of the protein is used",
        "blocks": dict(Counter(np.array(context_blocks)[no_taxonomy])),
        "sources": ["operon", "operon_gene", "neighbor_component", "ro_etc",
                    "ro_regulation", "replicon.is_plasmid"]}
    results["q3_substrate_class_from_context"] = q3

    print("\n--- Soru 3 ablasyonlari")
    ablations = {}
    for name, mask in (("context_without_taxonomy", no_taxonomy),
                       ("taxonomy_only", only_taxonomy),
                       ("context_plus_taxonomy", np.ones_like(no_taxonomy))):
        subset, _, _, _, _ = evaluate_task(
            x_context_all[:, mask], list(np.array(context_names)[mask]),
            substrate, groups, substrate_classes, rng, args, with_null=False,
            label=f"ablation:{name}")
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
        "neighbourhood reaches, 'the neighbourhood predicts the chemistry' is "
        "not established -- what is established is that the lineage does")
    results["q3_ablations"] = ablations

    q3_importance = importance_block(
        x_context_all[:, no_taxonomy],
        list(np.array(context_names)[no_taxonomy]), onehot3, weight3, grouped3,
        rng, args, motif_columns, sequence=False)
    results["q3_informative_context_features"] = q3_importance

    # --- Soru 4: kontroller ve modalite karsilastirmasi
    print("\n=== Soru 4a: POZITIF KONTROL -- HMM grubu <- dizi")
    positive, _, _, _, _ = evaluate_task(
        x_sequence, sequence_names, ro_group, groups, group_classes, rng, args,
        label="positive-control")
    positive["why_this_control"] = (
        "ro_group is the coarse HMM group of the type, so it is itself derived "
        "from sequence similarity. A model that can recover it across held-out "
        "types proves the grouped scheme is not blind: a low accuracy elsewhere "
        "then says something about the biology, not about the method. It is NOT "
        "an independent biological result")
    results["q4a_positive_control_ro_group_from_sequence"] = positive

    print("\n=== Soru 4b: modalite karsilastirmasi")
    modality = {}
    for label, matrix_x, names_x, target, classes in (
            ("reaction_class_from_sequence", x_sequence, sequence_names,
             reaction, reaction_classes),
            ("reaction_class_from_context", x_context_all[:, no_taxonomy],
             list(np.array(context_names)[no_taxonomy]), reaction,
             reaction_classes),
            ("substrate_class_from_sequence", x_sequence, sequence_names,
             substrate, substrate_classes),
            ("substrate_class_from_context", x_context_all[:, no_taxonomy],
             list(np.array(context_names)[no_taxonomy]), substrate,
             substrate_classes)):
        subset, _, _, _, _ = evaluate_task(
            matrix_x, names_x, target, groups, classes, rng, args,
            with_null=False, label=f"modality:{label}")
        modality[label] = {
            "grouped_accuracy": subset["grouped_tree"]["accuracy"],
            "grouped_balanced_accuracy":
                subset["grouped_tree"]["balanced_accuracy"],
            "type_level_accuracy": subset["grouped_tree"]["type_level_accuracy"],
            "global_majority_rate":
                subset["baselines"]["global_majority_rate"],
            "leaked_random_accuracy": subset["leaked_random_tree"]["accuracy"]}
    modality["note"] = ("the same folds, the same model and the same baselines "
                        "for every cell, so the four numbers can be compared "
                        "with each other")
    results["q4b_modality_comparison"] = modality

    print("\n=== Soru 4c: tipin ozelligi OLMAYAN bir etiket (beta alt birimi)")
    beta_label = []
    beta_lookup = dict(connection.execute(
        "SELECT candidate_id, has_beta FROM ro_etc"))
    for entry in entries:
        beta_label.append("alpha3beta3" if beta_lookup.get(entry["candidate_id"])
                          else "alpha3")
    beta, _, _, _, _ = evaluate_task(
        x_sequence, sequence_names, beta_label, groups, ["alpha3", "alpha3beta3"],
        rng, args, label="beta-subunit")
    beta["why_this_question"] = (
        "every other label in this file is a property of the type, so the "
        "effective sample size is 61. Whether a beta subunit sits next to the "
        "gene varies WITHIN a type (mean within-type entropy 0.27 bit against "
        "0.76 bit overall), so this is the one question where the entries carry "
        "information of their own. It also has a structural meaning: alpha3 "
        "against alpha3beta3 oligomer architecture")
    beta["limit"] = (
        "the label comes from annotation of the neighbouring gene, so a missing "
        "beta may mean a missing annotation rather than a missing subunit; the "
        "accuracy is therefore an underestimate of the sequence signal and an "
        "overestimate of how clean the label is")
    results["q4c_beta_subunit_from_sequence"] = beta

    # --- Bilgi bloklari
    results["method"] = {
        "why_this_module_exists": (
            "the data set is large enough to invite a sequence to reaction "
            "classifier, and a careless one would report an accuracy near 1.0 "
            "that means nothing"),
        "leakage": (
            "entries are not independent. 11,422 entries hold 10,019 distinct "
            "sequences, and the reaction class and the substrate class are "
            "curated per TYPE, not per entry. Under a random split a model can "
            "recognise the type of a test sequence from its near-identical "
            "training twins and read the label off it"),
        "how_it_was_avoided": (
            "every headline number comes from grouped cross-validation with the "
            f"group being {GROUP_FIELD}; GroupKFold with "
            f"{args.group_folds} entry-balanced folds plus leave-one-type-out. "
            "The random-split number is reported next to it and labelled "
            "invalid, because the size of the gap is itself a result"),
        "weights": (
            "each entry carries weight 1/(size of its type) so that a type with "
            "1,345 members and a type with 2 members have the same say; the "
            "unweighted variant is reported as a sensitivity"),
        "models": (
            "a depth-limited CART and an L1-penalised multinomial logistic "
            "regression, both implemented here with numpy because scikit-learn "
            "is not installed in this project environment. A random forest is "
            "fitted ONLY to rank features and never used for an accuracy"),
        "effective_sample_size": (
            "61 types, not 11,422 entries. For a 7-class problem this is small: "
            "no claim about the small classes can be supported"),
        "primary_metric": (
            "two numbers are given everywhere: accuracy over entries, and "
            "type-level accuracy where each type is scored once by the majority "
            "vote of its members. The second is the honest one"),
        "seed": SEED,
    }
    results["cannot_be_supported"] = [
        "no claim about a reaction class represented by 2 or 3 types "
        "(C_N_cleavage, O_demethylation, angular_dioxygenation): the grouped "
        "scheme gives those classes one or two test folds in total",
        "no claim that the columns found here are substrate-binding pocket "
        "residues. There is no structure in this database; what was measured is "
        "which alignment columns separate the classes",
        "no causal claim in either direction for question 3. Operon composition "
        "and substrate class are both properties of a lineage with a history; a "
        "shared genome is not evidence of a shared pathway",
        "no claim that the curated labels are correct. 58 % of entries sit in "
        "the distant or novel evidence tiers, where the reaction is unknown and "
        "was inherited from the type. This module treats the label as given",
        "no sequence-to-substrate prediction for a new protein. Even the best "
        "grouped accuracy here is far from a usable decision rule, and the "
        "reference pair calibration already showed that no global identity "
        "threshold guarantees the same substrate",
    ]
    results["inputs"] = {
        "db": os.path.abspath(args.db),
        "alignment": os.path.abspath(args.alignment),
        "chemistry": os.path.abspath(args.chemistry),
        "ecology": os.path.abspath(args.ecology),
        "sdp_positions": os.path.abspath(args.sdp),
        "motif_stats": motif_source,
        "entries": len(entries), "types": len(set(groups)),
        "alignment_match_columns": int(matrix.shape[1]),
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
    line = "=" * 76
    print("\n" + line)
    print("VERIDEN OGRENME -- GRUPLANMIS (TIP BAZLI) CAPRAZ DOGRULAMA")
    print(line)

    def block(key, title):
        task = results.get(key)
        if not task:
            return
        grouped = task["grouped_tree"]
        leaked = task["leaked_random_tree"]
        base = task["baselines"]
        print(f"\n{title}")
        print(f"  siniflar: {len(task['classes'])}   tip: {task['n_types']}   "
              f"giris: {task['n_entries']}   ozellik: {task['n_features']}")
        print(f"  {'gruplanmis dogruluk (tip bazli CV)':48} "
              f"{grouped['accuracy']:.4f}")
        print(f"  {'  tip duzeyi dogruluk':48} "
              f"{grouped['type_level_accuracy']:.4f}")
        print(f"  {'  dengeli dogruluk (makro duyarlilik)':48} "
              f"{grouped['balanced_accuracy']:.4f}")
        if "leave_one_type_out_tree" in task:
            print(f"  {'her tipi disarida birak':48} "
                  f"{task['leave_one_type_out_tree']['accuracy']:.4f}")
        print(f"  {'L1 lojistik, gruplanmis':48} "
              f"{task['grouped_l1_logistic']['accuracy']:.4f}")
        print(f"  {'SIZINTILI rastgele bolme (GECERSIZ)':48} "
              f"{leaked['accuracy']:.4f}")
        print(f"  {'  sizinti ucurumu (kat)':48} "
              f"{task['leakage_gap']['leaked_over_grouped']:.2f}x")
        print(f"  {'taban: global cogunluk sinifi':48} "
              f"{base['global_majority_rate']:.4f}  "
              f"({base['global_majority_class']})")
        print(f"  {'taban: tip duzeyi cogunluk':48} "
              f"{base['type_level_majority_rate']:.4f}")
        print(f"  {'taban: kat basina egitim cogunlugu':48} "
              f"{base['grouped_train_majority_accuracy']:.4f}")
        null = task.get("permutation_null")
        if null:
            print(f"  {'null (tip duzeyi etiket karistirma)':48} "
                  f"{null['mean']:.4f} +- {null['sd']:.4f}  "
                  f"(p95 {null['p95']:.4f}, azami {null['max']:.4f}, "
                  f"{null['reps']} tekrar)")
            z = task.get("grouped_vs_null", {}).get("z_against_null")
            if z is not None:
                print(f"  {'  null karsisinda z':48} {z}")
        print("  sinif basina duyarlilik:")
        for name, stats in sorted(grouped["per_class"].items(),
                                  key=lambda kv: -kv[1]["support"]):
            print(f"    {name:30} n={stats['support']:6}  "
                  f"recall={stats['recall']:.3f}  "
                  f"precision={stats['precision']:.3f}")
        print("  karisiklik matrisi (satir = gercek sinif):")
        width = max(len(c) for c in task["classes"]) + 1
        print("    " + " " * width
              + " ".join(f"{c[:6]:>6}" for c in task["classes"]))
        for name, row in zip(task["classes"], grouped["confusion_matrix"]):
            print(f"    {name:{width}}" + " ".join(f"{v:>6}" for v in row))

    block("q1_reaction_class_from_sequence",
          "SORU 1 -- reaksiyon sinifi <- dizi (hizalama kolonlari)")
    block("q3_substrate_class_from_context",
          "SORU 3 -- substrat sinifi <- genomik baglam")
    block("q4a_positive_control_ro_group_from_sequence",
          "SORU 4a -- POZITIF KONTROL: HMM grubu <- dizi")
    block("q4c_beta_subunit_from_sequence",
          "SORU 4c -- beta alt birimi <- dizi (tip ici degisen etiket)")

    sweep = results.get("q1_reaction_class_from_sequence", {}).get("depth_sweep")
    if sweep:
        print("\nDERINLIK TARAMASI -- sizinti ucurumu kapasiteyle genisliyor")
        print(f"  {'derinlik':>9} {'gruplanmis':>12} {'sizintili':>12} {'kat':>7}")
        for row in sweep["rows"]:
            print(f"  {row['max_depth']:>9} {row['grouped_accuracy']:>12.4f} "
                  f"{row['leaked_random_accuracy']:>12.4f} {row['ratio']:>7.2f}")

    variants = results.get("q1_column_selection_variants")
    if variants:
        print("\nKOLON SECIMI -- hazir SDP listesi ne kadar sisiriyor")
        for row in variants["rows"]:
            print(f"  {row['selection']:42} {row['n_columns']:>5} kolon  "
                  f"{row['grouped_accuracy']:.4f}  sizinti: {row['leak'][:3]}")

    q2 = results.get("q2_informative_columns")
    if q2:
        print("\nSORU 2 -- en bilgilendirici hizalama kolonlari")
        print(f"  {'kolon':>6} {'kalinti':>8} {'bolge':>17} "
              f"{'merkeze uzaklik':>17} {'agac kazanci':>13}")
        for row in q2["top_features"][:15]:
            text = row["feature"].replace("column ", "").replace(" residue ", " ")
            column, residue = text.split()
            print(f"  {column:>6} {residue:>8} {row['zone']:>17} "
                  f"{row['offset_from_nearest']:>+9} "
                  f"({row['nearest_catalytic_centre_column']:>3})"
                  f"{row['tree_gain']:>13.5f}")
        print(f"  bolge dagilimi: {q2['zone_counts_in_top_features']}")
        print(f"  katalitik merkez kolonuna +-10 icinde: "
              f"{q2['top_features_within_10_columns_of_a_centre_column']}"
              f"/{len(q2['top_features'])}")
        print(f"  L1 sifir olmayan katsayi: {q2['l1_nonzero_features']}"
              f"/{q2['l1_total_features']}")

    ablations = results.get("q3_ablations")
    if ablations:
        print("\nSORU 3 ABLASYONLARI -- baglam mi, taksonomi mi")
        for name in ("context_without_taxonomy", "taxonomy_only",
                     "context_plus_taxonomy"):
            row = ablations.get(name)
            if row:
                print(f"  {name:30} {row['n_features']:>4} ozellik  "
                      f"gruplanmis {row['grouped_accuracy']:.4f}  "
                      f"dengeli {row['grouped_balanced_accuracy']:.4f}  "
                      f"tip {row['type_level_accuracy']:.4f}")

    modality = results.get("q4b_modality_comparison")
    if modality:
        print("\nSORU 4b -- bilgi nerede")
        print(f"  {'baslik':38} {'gruplanmis':>11} {'taban':>8} {'sizintili':>11}")
        for name, row in modality.items():
            if name == "note":
                continue
            print(f"  {name:38} {row['grouped_accuracy']:>11.4f} "
                  f"{row['global_majority_rate']:>8.4f} "
                  f"{row['leaked_random_accuracy']:>11.4f}")

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
