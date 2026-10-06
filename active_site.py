"""
ROAR-DB -- kristal yapisi olan tiplerde AKTIF BOLGE cebinin cikarimi.

NEDEN BU MODUL VAR
  Veritabaninin merkezi bulgusu su: bu ailede substrat secimi global dizi
  benzerligiyle belirlenmiyor. EdoA1 ile cumA1 %99,8 kimlikte ama biri
  etilbenzen biri kumen okside ediyor; ayni substrati paylasan en uzak cift
  %32,4 kimlikte. Demek ki karar az sayida pozisyonda veriliyor.
  variant_signature.py (adim 12b) o pozisyonlari ISTATISTIKSEL olarak ariyor:
  varyantlari ayiran hizalama kolonlari. Ama o modulun kendi uyarisi acik --
  "bunlar yapisal olarak dogrulanmis baglanma cebi kalintilari DEGILDIR;
  veritabaninda yapi yok". Bu modul tam o boslugu kapatir: 17 tipin
  kuratorlu PDB kimligi var, yani o 17 tip icin ceb TAHMIN EDILMEK zorunda
  degil, OLCULEBILIR.

NE OLCULDU
  1. Yapi RCSB'den cekilir ve structures/ altinda onbelleklenir.
  2. Katalitik mononukleer demir, Rieske [2Fe-2S] demirlerinden AYIRT EDILIR.
     Ayrim geometrik ve kimyasal, varsayima dayanmiyor:
       - Rieske demiri: 2 inorganik kopru kukurdu (<=3,0 A) ve <=3,5 A'da
         ortak bir ikinci demir. Cift, bir yanda 2 Cys SG, diger yanda 2 His
         N ile baglanir.
       - Katalitik demir: inorganik kukurt YOK, yakininda baska demir YOK,
         2 His azotu + 1 karboksilat oksijeni (2-His-1-karboksilat yuz triadi).
     Her yapi icin bu ayrimin hangi sayilarla yapildigi JSON'a yazilir.
     Guvenle ayirt edilemeyen yapi TAHMIN EDILMEZ, "undetermined" raporlanir.
  3. Cebi doseyen kalintilar: katalitik demire belirtilen yaricap icinde
     herhangi bir atomu olan amino asitler. Iki yaricapta birden raporlanir
     (bkz. POCKET_RADII_A), cunku tek bir keyfi yaricapa guvenmek yerine
     okurun duyarliligi gormesi gerekir.
  4. Kalintilar PDB yazar numaralamasindan HIZALAMA KOLONUNA tasinir. Yazar
     numaralamasi dizi indeksiyle ayni DEGILDIR, bu yuzden yapinin kendi
     dizisi referans diziye hizalanir (BLOSUM62, serbest uc bosluklari) ve
     kolon haritasi oradan kurulur. Haritanin dogrulugu katalitik triad ve
     Rieske ligand kolonlari uzerinde SINANIR: kolon 212/217 gercekten iki
     His'e, kolon 355 karboksilata denk gelmeli VE o kalintilar yapida
     gercekten demiri bagliyor olmali. Sinanmayan harita raporlanmaz.
  5. Her ceb kolonu icin o tipin UYELERINDEKI kalinti dagilimi: hangi ceb
     pozisyonu degismez, hangisi varyantlar arasinda oynuyor.
  6. "Basit elektrostatik" fikri: Poisson-Boltzmann YAPILMADI ve yapilamaz
     (uzman yazilim, protonasyon atamasi, cozucu modeli gerekir). Yapilan sey
     notr pH'ta FORMEL YUK KOMPOZISYONU ozetidir: asidik/bazik/polar/hidrofobik
     sayimlari, net yuk, aromatik sayisi. JSON bunu acik acik bir kompozisyon
     ozeti olarak etiketler, elektrostatik hesap olarak degil.

NE OLCULMEDI / NE SOYLENEMEZ
  - Baglanma enerjisi, elektrostatik potansiyel, pKa, protonasyon durumu yok.
  - Ceb kalintilari KRISTALDEKI apo/holo konformasyondan okunur; substrat
    baglanirken yan zincir hareketi hesaba katilmaz.
  - 17 tip 13 ayri PDB kaydindan gelir (1NDO uc tipe, 1Z03 ve 1WW9 ikiser
    tipe bagli). Bu yuzden tip bazinda sayim BAGIMSIZ gozlem sayisi degildir;
    her toplam icin ayri bir "independent" sayimi verilir.
  - Substrat sinifi x kompozisyon karsilastirmasinda grup basina n = 1-4.
    Betimleyici, istatistik degil.

Girdi : chemistry.csv (pdb kolonu), cluster_ecology.csv (substrate_class),
        ROs_71_Clean/refs71.fasta, ROs_71_Clean/refs71_hmmaln.sto,
        genomic_context/cand_aln.sto, roar.sqlite (ro, ro_leaf, cluster_sdp)
Cikti : analysis_out/active_site.json
"""

import argparse
import csv
import json
import math
import os
import hashlib
import sqlite3
import sys
import time
import urllib.error
import urllib.request
import warnings
from collections import Counter, defaultdict
from datetime import datetime, timezone

from ro_motif import (COLUMN_WINDOW, DEFAULT_MODEL, MODEL_COLUMNS,
                      read_stockholm_matchcols)

# --- Indirme: iyi vatandaslik ---
# Tek seferde tek istek, her indirmeden sonra kisa bir bekleme, ve kim oldugunu
# soyleyen bir User-Agent. Onbellek sayesinde tekrar kosuslarda hic istek gitmez.
RCSB_CIF_URL = "https://files.rcsb.org/download/{0}.cif"
USER_AGENT = ("ROAR-DB/1.0 active_site.py (Rieske oxygenase database; "
              "academic use; contact via repository)")
DOWNLOAD_PAUSE_S = 1.5
HTTP_TIMEOUT_S = 60
STRUCTURE_CACHE_DIR = "structures"

# --- Demir atama esikleri ---
# Hepsi koordinasyon kimyasindan, yuvarlak sayi secimi degil:
#   Fe-N(His)        2,0-2,2 A
#   Fe-O(karboksilat) 2,0-2,3 A
#   Fe-S(Cys)        2,2-2,3 A
# 2,8 A tavani bu uc baga rahat yer acar ama ikinci kabuga (>3 A) tasmaz.
FE_LIGAND_MAX_A = 2.8
# [2Fe-2S] merkezinde Fe-Fe mesafesi 2,65-2,80 A olcilur (1NDO'da 2,69-2,77).
# 3,5 A tavani bu cifti yakalar; mononukleer demirin yaninda baska demir yoktur.
FE_FE_MAX_A = 3.5
# Kopru kukurdu Fe-S 2,15-2,35 A. 3,0 A tavani kristalografik hatayi tolere eder.
INORGANIC_S_MAX_A = 3.0
# Katalitik merkez icin asgari sart: 2 His azotu + en az 1 karboksilat oksijeni.
MIN_HIS_NITROGEN_FOR_TRIAD = 2
MIN_CARBOXYLATE_OXYGEN_FOR_TRIAD = 1
# Rieske cifti icin asgari sart: iki demirin toplaminda 2 Cys SG + 2 His N.
MIN_CYS_FOR_RIESKE_PAIR = 2
MIN_HIS_FOR_RIESKE_PAIR = 2
# Kopru karboksilati (ro_motif kolon 209) metali BAGLAMAZ -- 1NDO'da olculdu:
# katalitik demire 6,03 A. Yaptigi sey katalitik His ile KOMSU alt birimin
# Rieske His'i arasinda hidrojen bagi kurmaktir. 3,5 A donor-akseptor mesafesi
# hidrojen bagi icin standart ust sinirdir. Bu yuzden o kolon icin "metali
# bagliyor mu" diye sorulmaz; geometrisi OLCULUP raporlanir, gecme/kalma
# testi yapilmaz.
HBOND_MAX_A = 3.5

# --- Ceb yaricaplari ---
# 5,0 A: BIRINCI KABUK. Demirin koordinasyon kuresi (~2,2 A) artik bir van der
#        Waals temasi (~3 A) kadar. Bu yaricap metal ligandlarini ve onlara
#        dogrudan degen kalintilari verir; substratin tamamini vermez.
# 8,0 A: SUBSTRAT CEBI. Bu enzimlerin substratlari naftalen (uzun eksen ~7 A),
#        bifenil (~9 A), kolesterol halkalari gibi uzun molekuller ve demire
#        yakin halka karbonundan uzak uca kadar olan mesafe 6-9 A'ya ulasir.
#        Yani baglanan substratla temas edebilecek kalintilarin kapsanmasi icin
#        yaricap bu araliga cikmak zorundadir.
# Ikisi birlikte raporlanir: 5,0 A "kesin olan", 8,0 A "ilgili olan". Aradaki
# fark okura yaricap secimine ne kadar bagli oldugunu gosterir.
POCKET_RADII_A = (5.0, 8.0)
PRIMARY_RADIUS_A = 8.0

# --- Yapi dizisi <-> referans dizi hizalamasi ---
# Global hizalama, uc bosluklari serbest (yapi dizisi referansin bir parcasi
# olabilir: kristalde gorulmeyen uclar ve ilmekler eksiktir).
ALIGN_GAP_OPEN = -11.0
ALIGN_GAP_EXTEND = -1.0
# Bir zincirin "alpha alt birimi" sayilmasi icin asgari sartlar. Beta alt
# birimi ve ferredoksin zincirleri bunu gecmez, yani kolon haritasi yanlis
# zincire kurulmaz.
MIN_CHAIN_IDENTITY = 0.30
MIN_CHAIN_LENGTH = 100
# Demiri tasiyan zincirin referansa kimligi bunun altindaysa harita
# RAPORLANMAZ: hizalama guvenilmezse kolon numaralari da guvenilmezdir.
MIN_OWNER_IDENTITY = 0.30

# --- Varyant analizi ---
# variant_signature.py ile ayni esik: bir varyantin konsensusundan soz
# edebilmek icin en az 3 uye gerekir.
MIN_LEAF_MEMBERS = 3
# Betimleyici korunma kovalari. Bunlar KARAR esigi degil; "kolonun en sik
# kalintisi uyelerin yuzde kacinda" sorusunun okunabilir ozetidir.
CONSERVATION_BINS = (1.0, 0.99, 0.90)

# --- Kompozisyon siniflari (notr pH, FORMEL yuk) ---
# His AYRI tutulur: imidazol pKa'si 6,0-7,0 arasindadir, yani notr pH'ta
# protonasyon durumu belirsizdir; ustelik katalitik His'ler demire bagli
# oldugu icin "serbest baz" gibi davranmaz. Onu bazik saymak net yuku
# gerekcesiz sisirirdi.
ACIDIC = "DE"
BASIC = "KR"
HISTIDINE = "H"
POLAR = "STNQCY"
HYDROPHOBIC = "AVLIMFWPG"
AROMATIC = "FWY"

AMINO_ACID_3TO1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    # modifiye kalintilar: dizi acisindan ana kalinti gibi davranirlar
    "MSE": "M", "CSO": "C", "CME": "C", "KCX": "K", "MLY": "K",
    "SEP": "S", "TPO": "T", "PTR": "Y", "HIC": "H", "NEP": "H",
}
CARBOXYLATE_OXYGENS = {("ASP", "OD1"), ("ASP", "OD2"),
                       ("GLU", "OE1"), ("GLU", "OE2")}
HISTIDINE_NITROGENS = {("HIS", "ND1"), ("HIS", "NE2")}
CYSTEINE_SULFUR = {("CYS", "SG")}


# ---------------------------------------------------------------- yardimcilar

def shannon_entropy(counter):
    total = sum(counter.values())
    if total <= 1:
        return 0.0
    return -sum((n / total) * math.log2(n / total)
                for n in counter.values() if n > 0)


def read_fasta(path):
    """Basit FASTA okuyucu. Donen: {baslik: dizi (buyuk harf)}"""
    out, name, chunks = {}, None, []
    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    out[name] = "".join(chunks).upper()
                name, chunks = line[1:].split()[0], []
            else:
                chunks.append(line)
    if name is not None:
        out[name] = "".join(chunks).upper()
    return out


def read_stockholm_full(path):
    """Stockholm'u TAM genislikte oku (insertion kolonlari dahil).

    read_stockholm_matchcols() insertion kolonlarini atar; bu dogru davranistir
    ama kalinti indeksinden kolona harita kurmak icin yetmez, cunku insertion
    state'teki kalintilar o dizgede hic gorunmez. Burada tam hizalama tutulur ve
    match kolonlari ro_motif ile AYNI kuralla belirlenir (ilk dizide '-' veya
    buyuk harf olan kolonlar). Ayni kural kullanildigi icin iki okuyucunun
    uyustugu main() icinde ayrica sinanir.
    """
    blocks, order = {}, []
    with open(path) as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line or line.startswith("#") or line.startswith("//"):
                continue
            parts = line.split()
            if len(parts) != 2:
                continue
            name, chunk = parts
            if name not in blocks:
                blocks[name] = []
                order.append(name)
            blocks[name].append(chunk)
    return {n: "".join(v) for n, v in blocks.items()}, order


def build_column_maps(sto_path, fasta_sequences):
    """Referans dizi indeksi -> hizalama match kolonu haritasi.

    DIKKAT -- bu adimda sessiz bir hata tuzagi var: refs71_hmmaln.sto
    'hmmalign --trim' ile uretilmis, yani hizalanan dizi FASTA dizisinin uc
    kisimlari KIRPILMIS halidir (ornek: NahAc'ta N-ucundan 4 kalinti, 449 -> 433).
    Hizalamadaki kalinti siralamasini dogrudan FASTA indeksi saymak butun
    kolon haritasini 4 kalinti kaydirirdi. Bu yuzden kirpilmis dizi FASTA
    dizisinin ICINDE aranir ve bulunan ofset haritaya eklenir. Ofset
    bulunamazsa o referans icin harita URETILMEZ, raporlanir.

    Donen: (maps, matchcol_sequences, problems)
      maps[isim] = {fasta_indeksi_0_tabanli: kolon_1_tabanli}
    """
    full, order = read_stockholm_full(sto_path)
    if not full:
        return {}, {}, ["reference alignment is empty"]
    first = full[order[0]]
    match_indices = [i for i, ch in enumerate(first)
                     if ch == "-" or ch.isupper()]
    match_set = set(match_indices)
    matchcols = {n: "".join(s[i] for i in match_indices)
                 for n, s in full.items()}

    maps, problems = {}, []
    for name, aligned in full.items():
        ungapped = "".join(c for c in aligned if c.isalpha()).upper()
        reference = fasta_sequences.get(name)
        if reference is None:
            problems.append("%s: aligned but absent from the reference FASTA"
                            % name)
            continue
        offset = reference.find(ungapped)
        if offset < 0:
            problems.append("%s: trimmed alignment is not a contiguous "
                            "substring of the FASTA sequence; column map "
                            "not built" % name)
            continue
        mapping, seen, column = {}, 0, 0
        for i, ch in enumerate(aligned):
            if i in match_set:
                column += 1
            if ch.isalpha():
                if i in match_set:
                    mapping[offset + seen] = column
                seen += 1
        maps[name] = mapping
    return maps, matchcols, problems


def reference_name_for_cluster(cluster, reference_names):
    """chemistry.csv 'cluster' degerinden FASTA basligini bul.

    Baslik tip kimligine bir sonek ekler ('3_309_NahAc_dioxygenase_pro').
    Baslik KESILMEZ, cunku bazi gen adlari harf-rakam karisimi; onun yerine
    kuratorlu tip kimligi ONEK olarak aranir.
    """
    hits = [n for n in reference_names
            if n == cluster or n.startswith(cluster + "_")]
    if len(hits) == 1:
        return hits[0], None
    if not hits:
        return None, "no reference sequence matches the type id as a prefix"
    return None, ("type id is an ambiguous prefix of %d reference headers: %s"
                  % (len(hits), ", ".join(sorted(hits))))


def motif_columns():
    """ro_motif.py'nin tanimladigi sekiz kolon + etiketleri."""
    model = MODEL_COLUMNS[DEFAULT_MODEL]
    out = []
    for column, accepted, label in model["rieske"]:
        out.append((column, accepted, label, "rieske_cluster_ligand"))
    for column, accepted, label in model["catalytic"]:
        out.append((column, accepted, label, "catalytic_iron_ligand"))
    column, accepted, label = model["bridging"]
    out.append((column, accepted, label, "inter_subunit_bridge"))
    return out


# ---------------------------------------------------------------- indirme

def fetch_structure(pdb_id, cache_dir, offline=False):
    """Yapi dosyasini onbellekten ver, yoksa RCSB'den indir.

    Donen: (yol veya None, durum dizgisi)
    """
    pdb_id = pdb_id.strip().upper()
    path = os.path.join(cache_dir, pdb_id + ".cif")
    if os.path.exists(path) and os.path.getsize(path) > 0:
        return path, "cached"
    if offline:
        return None, "absent from cache and --offline was given"
    os.makedirs(cache_dir, exist_ok=True)
    request = urllib.request.Request(RCSB_CIF_URL.format(pdb_id),
                                     headers={"User-Agent": USER_AGENT})
    try:
        with urllib.request.urlopen(request, timeout=HTTP_TIMEOUT_S) as resp:
            payload = resp.read()
    except (urllib.error.URLError, urllib.error.HTTPError, OSError) as exc:
        return None, "download failed: %s" % exc
    if not payload:
        return None, "download returned an empty file"
    temporary = path + ".part"
    with open(temporary, "wb") as handle:
        handle.write(payload)
    os.replace(temporary, path)
    # Sunucuya nazik davran: her indirmeden sonra bekle (tek seferde tek istek).
    time.sleep(DOWNLOAD_PAUSE_S)
    return path, "downloaded"


# ---------------------------------------------------------------- yapi analizi

def residue_label(residue):
    """'HIS208' gibi okunabilir etiket (insertion kodu varsa eklenir)."""
    hetero, number, icode = residue.id
    suffix = icode.strip()
    return "%s%d%s" % (residue.get_resname(), number, suffix)


def residue_auth_id(residue):
    """Yazar numarasi, insertion kodu ile birlikte dizgi olarak."""
    _, number, icode = residue.id
    return "%d%s" % (number, icode.strip())


def is_amino_acid(residue):
    return residue.get_resname().upper() in AMINO_ACID_3TO1


def iron_atoms(model):
    out = []
    for chain in model:
        for residue in chain:
            for atom in residue:
                if (atom.element or "").upper() == "FE":
                    out.append((chain, residue, atom))
    return out


def describe_iron(atom, neighbor_search):
    """Bir demir atomunun koordinasyon cevresini say.

    Ayrimin tamami bu sayilarda: inorganik kopru kukurdu ve ortak demir
    [2Fe-2S] merkezini isaretler; onlarin YOKLUGU artik 2 His + karboksilat
    mononukleer katalitik merkezi isaretler.
    """
    radius = max(FE_FE_MAX_A, INORGANIC_S_MAX_A, FE_LIGAND_MAX_A)
    counts = {"inorganic_sulfur": 0, "partner_iron": 0,
              "cysteine_sulfur": 0, "histidine_nitrogen": 0,
              "carboxylate_oxygen": 0}
    detail = {"inorganic_sulfur": [], "partner_iron": [],
              "cysteine_sulfur": [], "histidine_nitrogen": [],
              "carboxylate_oxygen": []}
    for other in neighbor_search.search(atom.coord, radius, level="A"):
        if other is atom:
            continue
        distance = float(atom - other)
        element = (other.element or "").upper()
        parent = other.get_parent()
        resname = parent.get_resname().upper()
        key = (resname, other.get_name().strip().upper())
        chain_id = parent.get_parent().id
        tag = "%s/%s:%s" % (chain_id, residue_label(parent),
                            other.get_name().strip())
        if element == "FE" and distance <= FE_FE_MAX_A:
            counts["partner_iron"] += 1
            detail["partner_iron"].append((tag, round(distance, 2)))
        elif element == "S" and resname not in AMINO_ACID_3TO1 \
                and distance <= INORGANIC_S_MAX_A:
            # Amino asit olmayan bir kalintinin kukurdu = inorganik kopru
            # kukurdu (FES, SF4, F3S ... ). Cys SG ve Met SD bu dala girmez.
            counts["inorganic_sulfur"] += 1
            detail["inorganic_sulfur"].append((tag, round(distance, 2)))
        elif distance <= FE_LIGAND_MAX_A:
            if key in CYSTEINE_SULFUR:
                counts["cysteine_sulfur"] += 1
                detail["cysteine_sulfur"].append((tag, round(distance, 2)))
            elif key in HISTIDINE_NITROGENS:
                counts["histidine_nitrogen"] += 1
                detail["histidine_nitrogen"].append((tag, round(distance, 2)))
            elif key in CARBOXYLATE_OXYGENS:
                counts["carboxylate_oxygen"] += 1
                detail["carboxylate_oxygen"].append((tag, round(distance, 2)))
    for value in detail.values():
        value.sort()
    return counts, detail


def classify_irons(model, neighbor_search):
    """Yapidaki her demiri siniflandir.

    Donen: liste; her oge dict (chain, residue, atom, role, counts, detail).
    role in {"rieske_cluster_iron", "other_iron_sulfur_cluster_iron",
             "catalytic_mononuclear_iron", "unassigned"}
    """
    out = []
    for chain, residue, atom in iron_atoms(model):
        counts, detail = describe_iron(atom, neighbor_search)
        if counts["inorganic_sulfur"] >= 2 and counts["partner_iron"] >= 1:
            role = "rieske_cluster_iron"
        elif counts["inorganic_sulfur"] >= 2:
            role = "other_iron_sulfur_cluster_iron"
        elif counts["inorganic_sulfur"] == 0 and counts["partner_iron"] == 0 \
                and counts["histidine_nitrogen"] >= MIN_HIS_NITROGEN_FOR_TRIAD \
                and counts["carboxylate_oxygen"] >= MIN_CARBOXYLATE_OXYGEN_FOR_TRIAD:
            role = "catalytic_mononuclear_iron"
        else:
            role = "unassigned"
        out.append({"chain": chain.id, "residue": residue, "atom": atom,
                    "label": residue_label(residue),
                    "atom_name": atom.get_name().strip(),
                    "role": role, "counts": counts, "detail": detail})
    return out


def rieske_pair_check(irons):
    """Rieske demirlerini ciftlere ayir ve 2 Cys + 2 His sartini sina.

    Beklenen yapi: cift, bir demirinde iki sisteine, diger demirinde iki
    histidine baglanir -- Rieske merkezini klasik ferredoksinden ayiran sey
    tam olarak bu. Sart saglanmazsa cift raporlanir ama ELENMEZ, cunku
    eksik modellenmis yan zincirler de ayni sonucu verebilir.
    """
    rieske = [i for i in irons if i["role"] == "rieske_cluster_iron"]
    used, pairs = set(), []
    for index, iron in enumerate(rieske):
        if index in used:
            continue
        partner_index = None
        best = None
        for other_index, other in enumerate(rieske):
            if other_index == index or other_index in used:
                continue
            distance = float(iron["atom"] - other["atom"])
            if distance <= FE_FE_MAX_A and (best is None or distance < best):
                best, partner_index = distance, other_index
        if partner_index is None:
            pairs.append({"irons": [iron["label"] + ":" + iron["atom_name"]],
                          "chain": iron["chain"],
                          "fe_fe_distance_A": None,
                          "cysteine_sulfur": iron["counts"]["cysteine_sulfur"],
                          "histidine_nitrogen": iron["counts"]["histidine_nitrogen"],
                          "pair_complete": False})
            used.add(index)
            continue
        partner = rieske[partner_index]
        used.update({index, partner_index})
        cys = iron["counts"]["cysteine_sulfur"] + partner["counts"]["cysteine_sulfur"]
        his = iron["counts"]["histidine_nitrogen"] + partner["counts"]["histidine_nitrogen"]
        pairs.append({
            "irons": sorted([iron["label"] + ":" + iron["atom_name"],
                             partner["label"] + ":" + partner["atom_name"]]),
            "chain": iron["chain"],
            "fe_fe_distance_A": round(float(best), 2),
            "cysteine_sulfur": cys,
            "histidine_nitrogen": his,
            "pair_complete": (cys >= MIN_CYS_FOR_RIESKE_PAIR
                              and his >= MIN_HIS_FOR_RIESKE_PAIR),
        })
    return pairs


def ligating_residue_roles(irons):
    """Hangi kalinti HANGI metali bagliyor. Donen: {zincir/etiket: {rol}}

    Kolon dogrulamasinda "metali bagliyor mu" sorusu yetmez: Rieske ligand
    kolonu Rieske demirini, katalitik ligand kolonu KATALITIK demiri baglamak
    zorundadir. Ikisini ayri tutmak kontrolu belirgin olarak sertlestirir.
    """
    out = defaultdict(set)
    for iron in irons:
        for key in ("histidine_nitrogen", "carboxylate_oxygen",
                    "cysteine_sulfur"):
            for tag, _ in iron["detail"][key]:
                out[tag.split(":")[0]].add(iron["role"])
    return out


def chain_sequence(chain):
    """Zincirin GOZLENEN dizisi + kalinti listesi (yazar numaralariyla)."""
    residues = [r for r in chain if is_amino_acid(r)]
    sequence = "".join(AMINO_ACID_3TO1[r.get_resname().upper()]
                       for r in residues)
    return sequence, residues


def make_aligner():
    from Bio import Align
    from Bio.Align import substitution_matrices
    aligner = Align.PairwiseAligner()
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = ALIGN_GAP_OPEN
    aligner.extend_gap_score = ALIGN_GAP_EXTEND
    # Uc bosluklari bedava: yapida gorulmeyen uclar referans diziyi
    # cezalandirmamali.
    aligner.target_end_gap_score = 0.0
    aligner.query_end_gap_score = 0.0
    return aligner


def sanitize(sequence):
    """BLOSUM62'de olmayan harfleri X'e cevir."""
    allowed = set("ACDEFGHIKLMNPQRSTVWY")
    return "".join(c if c in allowed else "X" for c in sequence.upper())


def align_chain_to_reference(aligner, chain_seq, reference_seq):
    """Yapi zinciri -> referans dizi indeks haritasi.

    Donen: (harita {zincir_indeksi: referans_indeksi}, kimlik, hizalanan_sayi)
    """
    if not chain_seq or not reference_seq:
        return {}, 0.0, 0
    alignments = aligner.align(sanitize(reference_seq), sanitize(chain_seq))
    best = alignments[0]
    mapping, identical, aligned = {}, 0, 0
    for (ref_start, ref_end), (chain_start, chain_end) in zip(*best.aligned):
        for offset in range(ref_end - ref_start):
            ref_index = ref_start + offset
            chain_index = chain_start + offset
            mapping[chain_index] = ref_index
            aligned += 1
            if reference_seq[ref_index].upper() == chain_seq[chain_index].upper():
                identical += 1
    denominator = min(len(chain_seq), len(reference_seq))
    identity = identical / denominator if denominator else 0.0
    return mapping, identity, aligned


def pocket_residues(model, neighbor_search, iron_atom, radius):
    """Demire 'radius' icinde HERHANGI bir atomu olan amino asitler.

    Su/iyon/ligand kalintilari disarida: olculen sey proteinin cebi.
    Donen: [(residue, min_mesafe)] mesafeye gore sirali.
    """
    best = {}
    for atom in neighbor_search.search(iron_atom.coord, radius, level="A"):
        residue = atom.get_parent()
        if not is_amino_acid(residue):
            continue
        distance = float(iron_atom - atom)
        key = id(residue)
        if key not in best or distance < best[key][1]:
            best[key] = (residue, distance)
    return sorted(best.values(), key=lambda item: (item[1],
                                                   item[0].get_parent().id,
                                                   item[0].id[1]))


# ---------------------------------------------------------------- kompozisyon

def composition_summary(one_letter_residues, ligand_one_letter):
    """Notr pH'ta FORMEL YUK KOMPOZISYONU ozeti -- elektrostatik hesap DEGIL.

    ligand_one_letter: demiri bagladigi yapisal olarak dogrulanmis kalintilar.
    Onlarin yuku ayri raporlanir, cunku demire koordine olan karboksilatin
    'serbest eksi yuk' gibi sayilmasi yanlis olurdu.
    """
    counts = Counter(one_letter_residues)
    ligands = Counter(ligand_one_letter)

    def tally(source, alphabet):
        return sum(source[c] for c in alphabet)

    acidic = tally(counts, ACIDIC)
    basic = tally(counts, BASIC)
    ligand_acidic = tally(ligands, ACIDIC)
    ligand_basic = tally(ligands, BASIC)
    total = sum(counts.values())
    return {
        "n_residues": total,
        "acidic_D_E": acidic,
        "basic_K_R": basic,
        "histidine": counts["H"],
        "polar_S_T_N_Q_C_Y": tally(counts, POLAR),
        "hydrophobic_A_V_L_I_M_F_W_P_G": tally(counts, HYDROPHOBIC),
        "aromatic_F_W_Y": tally(counts, AROMATIC),
        "aromatic_including_histidine": tally(counts, AROMATIC) + counts["H"],
        "net_formal_charge": basic - acidic,
        "net_formal_charge_excluding_iron_ligands":
            (basic - ligand_basic) - (acidic - ligand_acidic),
        "fraction_hydrophobic": (round(tally(counts, HYDROPHOBIC) / total, 3)
                                 if total else None),
        "fraction_charged": (round((acidic + basic) / total, 3)
                             if total else None),
        "residue_counts": dict(sorted(counts.items())),
    }


# ---------------------------------------------------------------- varyantlar

def column_variation(column, member_sequences):
    """Bir hizalama kolonunda uyelerdeki kalinti dagilimi."""
    raw = Counter()
    for sequence in member_sequences:
        raw[sequence[column - 1] if len(sequence) >= column else "-"] += 1
    gaps = raw.pop("-", 0)
    gaps += raw.pop(".", 0)
    observed = sum(raw.values())
    if observed == 0:
        return {"n_members_with_residue": 0, "n_members_with_gap": gaps,
                "distribution": {}, "top_residue": None,
                "top_fraction": None, "entropy_bits": None,
                "n_distinct_residues": 0, "invariant": None}
    top_residue, top_count = raw.most_common(1)[0]
    return {
        "n_members_with_residue": observed,
        "n_members_with_gap": gaps,
        "distribution": dict(sorted(raw.items(),
                                    key=lambda kv: (-kv[1], kv[0]))),
        "top_residue": top_residue,
        "top_fraction": round(top_count / observed, 4),
        "entropy_bits": round(shannon_entropy(raw), 4),
        "n_distinct_residues": len(raw),
        "invariant": top_count == observed,
    }


def variant_breakdown(columns, leaves, aligned):
    """Ceb kolonlarinda VARYANT (leaf) duzeyi dagilim.

    Her yaprak icin kolon konsensusu alinir; sonra kolon basina kac farkli
    varyant durumu oldugu ve her durumu kac varyantin tasidigi raporlanir.
    """
    usable = {leaf_id: members for leaf_id, members in leaves.items()
              if len(members) >= MIN_LEAF_MEMBERS}
    per_column = {}
    for column in columns:
        states = Counter()
        for leaf_id, members in usable.items():
            counter = Counter(aligned[m][column - 1] for m in members
                              if len(aligned[m]) >= column
                              and aligned[m][column - 1] not in "-.")
            if sum(counter.values()) < MIN_LEAF_MEMBERS:
                continue
            states[counter.most_common(1)[0][0]] += 1
        per_column[str(column)] = {
            "n_variants_scored": sum(states.values()),
            "variant_consensus_states": dict(
                sorted(states.items(), key=lambda kv: (-kv[1], kv[0]))),
            "n_distinct_variant_states": len(states),
        }
    return {
        "n_variants_total": len(leaves),
        "n_variants_with_enough_members": len(usable),
        "min_members_per_variant": MIN_LEAF_MEMBERS,
        "per_column": per_column,
    }


# ---------------------------------------------------------------- ana akis

def analyse_structure(cluster, pdb_id, structure_path, reference_name,
                      reference_sequence, column_map, aligner, parser):
    """Tek bir (tip, yapi) cifti icin tum olcumler. Donen: dict."""
    record = {
        "type": cluster,
        "pdb_id": pdb_id,
        "reference_name": reference_name,
        "reference_length": len(reference_sequence),
        # Kuratorlu set birebir ayni diziyi iki adla iceriyor (bkz. README:
        # 3_309_NahAc = 3_315_NDO). Bagimsiz gozlem sayarken ADA degil DIZIYE
        # bakmak gerekir, yoksa ayni yapi ayni diziyle iki kez sayilir.
        "reference_sequence_sha1_12": hashlib.sha1(
            reference_sequence.encode()).hexdigest()[:12],
        "status": "ok",
        "problems": [],
    }
    from Bio.PDB import NeighborSearch

    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            structure = parser.get_structure(pdb_id, structure_path)
    except Exception as exc:                      # noqa: BLE001
        record["status"] = "undetermined"
        record["problems"].append("mmCIF parsing failed: %s" % exc)
        return record

    model = next(structure.get_models())
    atoms = [a for a in model.get_atoms()]
    neighbor_search = NeighborSearch(atoms)

    irons = classify_irons(model, neighbor_search)
    roles = Counter(i["role"] for i in irons)
    pairs = rieske_pair_check(irons)
    catalytic = [i for i in irons if i["role"] == "catalytic_mononuclear_iron"]

    record["iron_inventory"] = {
        "n_iron_atoms": len(irons),
        "n_rieske_cluster_irons": roles.get("rieske_cluster_iron", 0),
        "n_other_iron_sulfur_cluster_irons":
            roles.get("other_iron_sulfur_cluster_iron", 0),
        "n_catalytic_mononuclear_irons": roles.get(
            "catalytic_mononuclear_iron", 0),
        "n_unassigned_irons": roles.get("unassigned", 0),
        "rieske_pairs": pairs,
        "irons": [{
            "chain": i["chain"], "residue": i["label"],
            "atom": i["atom_name"], "assigned_role": i["role"],
            "inorganic_sulfur_within_%.1fA" % INORGANIC_S_MAX_A:
                i["counts"]["inorganic_sulfur"],
            "partner_iron_within_%.1fA" % FE_FE_MAX_A:
                i["counts"]["partner_iron"],
            "cysteine_sulfur_within_%.1fA" % FE_LIGAND_MAX_A:
                i["counts"]["cysteine_sulfur"],
            "histidine_nitrogen_within_%.1fA" % FE_LIGAND_MAX_A:
                i["counts"]["histidine_nitrogen"],
            "carboxylate_oxygen_within_%.1fA" % FE_LIGAND_MAX_A:
                i["counts"]["carboxylate_oxygen"],
            "contacts": {k: v for k, v in i["detail"].items() if v},
        } for i in sorted(irons, key=lambda x: (x["chain"], x["label"],
                                                x["atom_name"]))],
    }

    # Ayrimin NASIL yapildigi, gercek sayilarla, yapi basina.
    pair_distances = [p["fe_fe_distance_A"] for p in pairs
                      if p["fe_fe_distance_A"] is not None]
    triad_shape = sorted({
        "%dHis+%dcarboxylate" % (i["counts"]["histidine_nitrogen"],
                                 i["counts"]["carboxylate_oxygen"])
        for i in catalytic})
    record["how_irons_were_distinguished"] = (
        "Rieske cluster irons: %d atoms in %d pair(s); each pair is bridged by "
        "two non-protein (inorganic) sulfur atoms within %.1f A and the two "
        "irons sit %s A apart, one ligated by cysteine sulfur and the other by "
        "histidine nitrogen (%d/%d pairs show the full 2-Cys + 2-His pattern). "
        "Catalytic mononuclear iron: %d atom(s) with zero inorganic sulfur and "
        "zero partner iron within %.1f A, ligated instead by %s, i.e. the "
        "2-His-1-carboxylate facial triad. The discriminating observation is "
        "the presence or absence of bridging inorganic sulfur plus a second "
        "iron, not residue numbering or sequence."
        % (roles.get("rieske_cluster_iron", 0), len(pairs),
           INORGANIC_S_MAX_A,
           ("%.2f-%.2f" % (min(pair_distances), max(pair_distances))
            if pair_distances else "n/a"),
           sum(1 for p in pairs if p["pair_complete"]), len(pairs),
           len(catalytic), FE_FE_MAX_A,
           ", ".join(triad_shape) if triad_shape else "n/a"))

    if not catalytic:
        record["status"] = "undetermined"
        record["problems"].append(
            "no iron satisfies the catalytic criteria (no inorganic sulfur and "
            "no partner iron nearby, at least %d histidine nitrogen and %d "
            "carboxylate oxygen within %.1f A); the catalytic site is not "
            "reported rather than guessed"
            % (MIN_HIS_NITROGEN_FOR_TRIAD, MIN_CARBOXYLATE_OXYGEN_FOR_TRIAD,
               FE_LIGAND_MAX_A))
        return record

    # --- zincirleri referansa hizala ---
    chain_info = {}
    for chain in model:
        sequence, residues = chain_sequence(chain)
        if len(sequence) < MIN_CHAIN_LENGTH:
            continue
        mapping, identity, aligned_n = align_chain_to_reference(
            aligner, sequence, reference_sequence)
        chain_info[chain.id] = {
            "sequence": sequence, "residues": residues,
            "chain_to_reference": mapping, "identity": identity,
            "n_aligned": aligned_n,
            "alpha_like": identity >= MIN_CHAIN_IDENTITY,
        }
    record["chains"] = [{
        "chain": cid,
        "observed_residues": len(info["sequence"]),
        "identity_to_reference": round(info["identity"], 4),
        "aligned_positions": info["n_aligned"],
        "treated_as_oxygenase_alpha_subunit": info["alpha_like"],
    } for cid, info in sorted(chain_info.items())]

    # --- temsilci katalitik demir secimi ---
    # Oncelik: koordinasyonu en tam olan; esitlikte referansa kimligi en yuksek
    # zincirdeki; son olarak zincir kimligi alfabetik (belirlenimli olsun diye).
    def catalytic_rank(iron):
        info = chain_info.get(iron["chain"])
        identity = info["identity"] if info else 0.0
        return (-(iron["counts"]["histidine_nitrogen"]
                  + iron["counts"]["carboxylate_oxygen"]),
                -identity, iron["chain"], iron["label"])

    ordered = sorted(catalytic, key=catalytic_rank)
    chosen = ordered[0]
    owner = chain_info.get(chosen["chain"])
    record["catalytic_iron"] = {
        "chain": chosen["chain"],
        "residue": chosen["label"],
        "n_equivalent_copies_in_asymmetric_unit": len(catalytic),
        "coordinating_contacts": {k: v for k, v in chosen["detail"].items()
                                  if v},
    }
    if owner is None or not owner["alpha_like"] \
            or owner["identity"] < MIN_OWNER_IDENTITY:
        record["status"] = "undetermined"
        record["problems"].append(
            "the chain carrying the catalytic iron (%s) aligns to the "
            "reference sequence at only %.1f%% identity, below the %.0f%% "
            "floor; PDB-to-alignment-column mapping would not be trustworthy "
            "so it is not reported"
            % (chosen["chain"],
               100.0 * (owner["identity"] if owner else 0.0),
               100.0 * MIN_OWNER_IDENTITY))
        return record
    record["catalytic_iron"]["chain_identity_to_reference"] = round(
        owner["identity"], 4)
    # Secilen demiri YAPISAL olarak bagladigi dogrulanmis kalintilar.
    # Kompozisyon ozetinde yuklerinin ayrica cikarilmasi icin gerekli.
    chosen_ligand_labels = {tag.split(":")[0].split("/")[1]
                            for key in ("histidine_nitrogen",
                                        "carboxylate_oxygen")
                            for tag, _ in chosen["detail"][key]}
    record["catalytic_iron"]["ligating_residues"] = sorted(chosen_ligand_labels)

    # --- ceb kalintilari, iki yaricapta ---
    index_of_residue = {id(r): i for i, r in enumerate(owner["residues"])}

    def map_residue(residue):
        """Kalintiyi hizalama kolonuna tasi. Donen: (kolon veya None, neden)"""
        chain_id = residue.get_parent().id
        if chain_id != chosen["chain"]:
            return None, "residue belongs to another chain (%s)" % chain_id
        index = index_of_residue.get(id(residue))
        if index is None:
            return None, "residue not part of the aligned chain sequence"
        reference_index = owner["chain_to_reference"].get(index)
        if reference_index is None:
            return None, "position is a gap in the structure-to-reference alignment"
        column = column_map.get(reference_index)
        if column is None:
            return None, ("reference residue %d falls in an insert state of "
                          "the profile, so it has no match-state column"
                          % (reference_index + 1))
        return column, None

    record["pocket"] = {}
    for radius in POCKET_RADII_A:
        rows, foreign, unmapped = [], [], []
        for residue, distance in pocket_residues(model, neighbor_search,
                                                 chosen["atom"], radius):
            chain_id = residue.get_parent().id
            one_letter = AMINO_ACID_3TO1[residue.get_resname().upper()]
            column, reason = map_residue(residue)
            entry = {
                "chain": chain_id,
                "pdb_residue": residue_label(residue),
                "author_number": residue_auth_id(residue),
                "residue_3letter": residue.get_resname().upper(),
                "residue_1letter": one_letter,
                "min_distance_to_iron_A": round(distance, 2),
                "alignment_column": column,
            }
            if chain_id != chosen["chain"]:
                entry["note"] = reason
                foreign.append(entry)
                continue
            if column is None:
                entry["note"] = reason
                unmapped.append(entry)
            rows.append(entry)
        mapped = [r for r in rows if r["alignment_column"] is not None]
        record["pocket"]["%.1f" % radius] = {
            "radius_A": radius,
            "n_residues_same_chain": len(rows),
            "n_residues_mapped_to_alignment_columns": len(mapped),
            "n_residues_unmapped": len(unmapped),
            "n_contacts_from_other_chains": len(foreign),
            "residues": rows,
            "contacts_from_other_chains": foreign,
            "alignment_columns": sorted({r["alignment_column"]
                                         for r in mapped}),
            "composition": composition_summary(
                [r["residue_1letter"] for r in rows],
                [r["residue_1letter"] for r in rows
                 if r["pdb_residue"] in chosen_ligand_labels]),
        }

    # --- kopyalar arasi ic tutarlilik ---
    # Asimetrik birimde birden fazla katalitik demir varsa hepsinin ayni kolon
    # kumesini vermesi gerekir. Vermiyorsa ya hizalama ya da secim kirilgandir;
    # bu yuzden olculup raporlanir.
    consistency = {}
    for radius in POCKET_RADII_A:
        column_sets = []
        for iron in ordered:
            info = chain_info.get(iron["chain"])
            if info is None or not info["alpha_like"]:
                continue
            local_index = {id(r): i for i, r in enumerate(info["residues"])}
            columns = set()
            for residue, _ in pocket_residues(model, neighbor_search,
                                              iron["atom"], radius):
                if residue.get_parent().id != iron["chain"]:
                    continue
                index = local_index.get(id(residue))
                if index is None:
                    continue
                reference_index = info["chain_to_reference"].get(index)
                if reference_index is None:
                    continue
                column = column_map.get(reference_index)
                if column is not None:
                    columns.add(column)
            column_sets.append(columns)
        if column_sets:
            union = set().union(*column_sets)
            intersection = set.intersection(*column_sets)
            consistency["%.1f" % radius] = {
                "n_copies_compared": len(column_sets),
                "columns_in_every_copy": len(intersection),
                "columns_in_any_copy": len(union),
                "jaccard": (round(len(intersection) / len(union), 4)
                            if union else None),
                "identical_across_copies": len(intersection) == len(union),
            }
    record["copy_consistency"] = consistency

    # --- HARITANIN DOGRULANMASI ---
    # Bu adim modulun en kolay hata yapilan yeri, bu yuzden ayrica sinanir:
    # ro_motif.py'nin sekiz kolonu yapida hangi kalintiya dusuyor, o kalinti
    # beklenen tipte mi, VE gercekten dogru metali bagliyor mu.
    #
    # Pencere neden var: ro_motif.py'nin kendisi motif testini +-COLUMN_WINDOW
    # (=2) penceresiyle yapar ve gerekcesi modulunde yazili -- bazi alt
    # ailelerde profil hizalamasi bir iki kolon kayiyor. Burada hem TAM kolon
    # hem pencereli sonuc ayri ayri raporlanir ve kayma miktari yazilir, yani
    # tolerans sessiz degil: olculup gosteriliyor.
    ligand_roles = ligating_residue_roles(irons)
    required_role = {
        "catalytic_iron_ligand": "catalytic_mononuclear_iron",
        "rieske_cluster_ligand": "rieske_cluster_iron",
    }
    column_to_residue = {}
    for index, residue in enumerate(owner["residues"]):
        reference_index = owner["chain_to_reference"].get(index)
        if reference_index is None:
            continue
        column = column_map.get(reference_index)
        if column is not None:
            column_to_residue[column] = residue

    def inspect(column, accepted, role):
        """Bir kolondaki kalintiyi degerlendir. Donen: dict veya None."""
        residue = column_to_residue.get(column)
        if residue is None:
            return None
        one_letter = AMINO_ACID_3TO1[residue.get_resname().upper()]
        tag = "%s/%s" % (residue.get_parent().id, residue_label(residue))
        roles = ligand_roles.get(tag, set())
        wanted = required_role.get(role)
        ligates = wanted in roles if wanted else bool(roles)
        return {
            "structure_residue": residue_label(residue),
            "structure_residue_1letter": one_letter,
            "residue_matches_expectation": one_letter in accepted,
            "ligates_the_expected_metal": ligates,
            "metal_roles_observed": sorted(roles),
            "passed": (one_letter in accepted) and ligates,
        }

    checks = []
    exact_pass = window_pass = 0
    shifts = []
    for column, accepted, label, role in sorted(motif_columns()):
        entry = {"column": column, "expected_residue": accepted, "role": role,
                 "label_in_ro_motif": label}
        if role == "inter_subunit_bridge":
            # Bu kolon metal ligandi DEGIL; asagida ayrica olculur.
            residue = column_to_residue.get(column)
            entry.update({
                "structure_residue": (residue_label(residue) if residue
                                      else None),
                "residue_matches_expectation": (
                    AMINO_ACID_3TO1[residue.get_resname().upper()] in accepted
                    if residue else None),
                "evaluated_as": "not a metal ligand; geometry is reported "
                                "separately under bridging_carboxylate",
            })
            checks.append(entry)
            continue
        strict = inspect(column, accepted, role)
        entry["exact_column"] = strict
        if strict and strict["passed"]:
            exact_pass += 1
            window_pass += 1
            entry["observed_at_column"] = column
            entry["offset_from_model_column"] = 0
            entry["passed"] = True
            checks.append(entry)
            continue
        # Pencere taramasi: en yakin kolondan baslayarak
        order = [o for pair in zip(range(-1, -COLUMN_WINDOW - 1, -1),
                                   range(1, COLUMN_WINDOW + 1))
                 for o in pair]
        found = None
        for offset in order:
            candidate = inspect(column + offset, accepted, role)
            if candidate and candidate["passed"]:
                found = (offset, candidate)
                break
        if found:
            offset, candidate = found
            window_pass += 1
            entry["observed_at_column"] = column + offset
            entry["offset_from_model_column"] = offset
            entry["within_window"] = candidate
            entry["passed"] = True
            shifts.append({"model_column": column, "role": role,
                           "observed_at_column": column + offset,
                           "offset": offset,
                           "structure_residue": candidate["structure_residue"]})
        else:
            entry["observed_at_column"] = None
            entry["offset_from_model_column"] = None
            entry["passed"] = False
        checks.append(entry)

    metal_checks = [c for c in checks if c["role"] != "inter_subunit_bridge"]
    catalytic_checks = [c for c in checks
                        if c["role"] == "catalytic_iron_ligand"]
    triad_exact = all(c.get("exact_column") and c["exact_column"]["passed"]
                      for c in catalytic_checks)
    triad_window = all(c.get("passed") for c in catalytic_checks)
    record["column_mapping_verification"] = {
        "what_is_checked": (
            "every metal-ligand column of the motif model must land on a "
            "residue of the expected type AND that residue must actually "
            "ligate the expected metal in this structure (a Rieske ligand "
            "column must ligate a Rieske iron, a catalytic ligand column must "
            "ligate the mononuclear iron). The catalytic triad columns "
            "212/217/355 are the decisive ones. The inter-subunit bridging "
            "column is excluded from this test because it is not a metal "
            "ligand; its geometry is measured separately."),
        "window_columns": COLUMN_WINDOW,
        "why_a_window": (
            "ro_motif.py applies the same plus-or-minus %d column window to "
            "its own motif test, because the profile alignment shifts by a "
            "column or two in some subfamilies. Both the exact-column and the "
            "windowed result are reported, and every shift is listed, so the "
            "tolerance is measured rather than assumed." % COLUMN_WINDOW),
        "n_metal_ligand_columns": len(metal_checks),
        "n_passed_at_exact_column": exact_pass,
        "n_passed_within_window": window_pass,
        "catalytic_triad_verified_at_exact_columns": triad_exact,
        "catalytic_triad_verified_within_window": triad_window,
        "column_shifts": shifts,
        "checks": checks,
    }

    # --- kopru karboksilatinin geometrisi (olcum, test degil) ---
    bridge_column = MODEL_COLUMNS[DEFAULT_MODEL]["bridging"][0]
    bridge_residue = column_to_residue.get(bridge_column)
    catalytic_ligands = {t for t, r in ligand_roles.items()
                         if "catalytic_mononuclear_iron" in r}
    rieske_ligands = {t for t, r in ligand_roles.items()
                      if "rieske_cluster_iron" in r}
    bridging = {
        "model_column": bridge_column,
        "structure_residue": (residue_label(bridge_residue)
                              if bridge_residue else None),
        "hydrogen_bond_cutoff_A": HBOND_MAX_A,
        "distance_from_catalytic_iron_A": None,
        "side_chain_contacts_to_histidine": [],
        "note": (
            "This residue does not coordinate a metal; in 1NDO its carboxylate "
            "sits 6.0 A from the catalytic iron. Its documented role is to "
            "hydrogen-bond between the catalytic histidine and the Rieske "
            "histidine of the neighbouring subunit, so what is reported is "
            "measured geometry, not a pass or fail."),
    }
    if bridge_residue is not None:
        bridging["distance_from_catalytic_iron_A"] = round(
            min(float(chosen["atom"] - a) for a in bridge_residue), 2)
        resname = bridge_residue.get_resname().upper()
        accepted = MODEL_COLUMNS[DEFAULT_MODEL]["bridging"][1]
        one_letter = AMINO_ACID_3TO1[resname]
        bridging["residue_matches_expectation"] = one_letter in accepted
        if one_letter not in accepted:
            # Kolon beklenen Asp/Glu DEGIL. ro_motif testi +-COLUMN_WINDOW
            # penceresi kullandigi icin karboksilat komsu kolonda olabilir;
            # bu yuzden pencere TARANIR ve bulunan sey raporlanir.
            nearby = []
            for offset in range(-COLUMN_WINDOW, COLUMN_WINDOW + 1):
                if offset == 0:
                    continue
                other = column_to_residue.get(bridge_column + offset)
                if other is None:
                    continue
                if AMINO_ACID_3TO1[other.get_resname().upper()] in accepted:
                    nearby.append({
                        "column": bridge_column + offset,
                        "offset": offset,
                        "structure_residue": residue_label(other),
                        "distance_from_catalytic_iron_A": round(
                            min(float(chosen["atom"] - a) for a in other), 2),
                    })
            bridging["carboxylate_found_within_window"] = nearby
            bridging["expectation_note"] = (
                "the exact bridging column does not carry an aspartate or "
                "glutamate in this reference; %d carboxylate residue(s) were "
                "found within the plus-or-minus %d column window"
                % (len(nearby), COLUMN_WINDOW))
        for atom in bridge_residue:
            if (resname, atom.get_name().strip().upper()) \
                    not in CARBOXYLATE_OXYGENS:
                continue
            for other in neighbor_search.search(atom.coord, HBOND_MAX_A,
                                                level="A"):
                parent = other.get_parent()
                if (parent.get_resname().upper(),
                        other.get_name().strip().upper()) \
                        not in HISTIDINE_NITROGENS:
                    continue
                tag = "%s/%s" % (parent.get_parent().id,
                                 residue_label(parent))
                if tag in catalytic_ligands:
                    partner = "ligates the catalytic mononuclear iron"
                elif tag in rieske_ligands:
                    partner = "ligates a Rieske cluster iron"
                else:
                    partner = "histidine that ligates no iron here"
                bridging["side_chain_contacts_to_histidine"].append({
                    "from_atom": atom.get_name().strip(),
                    "to_residue": tag,
                    "to_atom": other.get_name().strip(),
                    "distance_A": round(float(atom - other), 2),
                    "partner_role": partner,
                    "same_chain_as_catalytic_iron":
                        parent.get_parent().id == chosen["chain"],
                })
        contacts = bridging["side_chain_contacts_to_histidine"]
        contacts.sort(key=lambda c: c["distance_A"])
        bridging["bridges_both_centres"] = (
            any("catalytic" in c["partner_role"] for c in contacts)
            and any("Rieske" in c["partner_role"] for c in contacts))
    record["bridging_carboxylate"] = bridging

    # --- merkezler arasi mesafe: demir atamasinin BAGIMSIZ capraz kontrolu ---
    # Katalitik demir ile en yakin Rieske demiri arasinda ~12 A beklenir ve bu
    # mesafe tipik olarak KOMSU alt birime gider; ayni zincir icinde 40 A'nin
    # uzerindedir. Bu sayi iki merkezin karistirilmadiginin bagimsiz kaniti.
    rieske_irons = [i for i in irons if i["role"] == "rieske_cluster_iron"]
    if rieske_irons:
        distances = sorted((round(float(chosen["atom"] - i["atom"]), 2),
                            i["chain"], i["label"] + ":" + i["atom_name"])
                           for i in rieske_irons)
        same_chain = [d for d in distances if d[1] == chosen["chain"]]
        record["centre_separation"] = {
            "nearest_rieske_iron_A": distances[0][0],
            "nearest_rieske_iron_chain": distances[0][1],
            "nearest_rieske_iron_is_in_the_same_chain":
                distances[0][1] == chosen["chain"],
            "nearest_rieske_iron_within_same_chain_A":
                same_chain[0][0] if same_chain else None,
            "note": (
                "Independent cross-check on the iron assignment: the catalytic "
                "iron and the Rieske cluster that feeds it sit about 12 A "
                "apart and that partner cluster usually belongs to the "
                "neighbouring subunit, while the Rieske cluster of the same "
                "chain is much further away. A short same-chain distance would "
                "mean the two centres had been confused."),
        }

    if not triad_window:
        record["problems"].append(
            "the catalytic triad columns did not all map onto residues of the "
            "expected type that ligate the mononuclear iron, even allowing the "
            "plus-or-minus %d column window; pocket column numbers for this "
            "structure should be treated as unverified" % COLUMN_WINDOW)
        record["status"] = "mapping_unverified"
    elif not triad_exact:
        record["problems"].append(
            "the catalytic triad is present and structurally verified, but "
            "not all of it at the exact model columns: %s. The pocket columns "
            "are usable; the shift is a property of the profile alignment in "
            "this subfamily, not a missing residue."
            % "; ".join("model column %d observed at %d (%s)"
                        % (s["model_column"], s["observed_at_column"],
                           s["structure_residue"]) for s in shifts))
    return record


def group_composition(records, key_function, label):
    """Tipleri bir anahtara gore grupla ve kompozisyon ortalamalarini ver."""
    buckets = defaultdict(list)
    for record in records:
        value = key_function(record)
        if value:
            buckets[value].append(record)
    out = {}
    for value, members in sorted(buckets.items()):
        pockets = [m["pocket"]["%.1f" % PRIMARY_RADIUS_A]["composition"]
                   for m in members]
        independent = {(m["pdb_id"], m["reference_sequence_sha1_12"])
                       for m in members}
        fields = ("n_residues", "acidic_D_E", "basic_K_R", "histidine",
                  "polar_S_T_N_Q_C_Y", "hydrophobic_A_V_L_I_M_F_W_P_G",
                  "aromatic_F_W_Y", "net_formal_charge")
        out[value] = {
            "n_types": len(members),
            "n_independent_structure_sequence_pairs": len(independent),
            "types": sorted(m["type"] for m in members),
            "mean": {f: round(sum(p[f] for p in pockets) / len(pockets), 2)
                     for f in fields},
            "range": {f: [min(p[f] for p in pockets),
                          max(p[f] for p in pockets)] for f in fields},
        }
    return {"grouped_by": label, "groups": out}


def main():
    parser_args = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser_args.add_argument("--db", default="roar.sqlite")
    parser_args.add_argument("--out-dir", default="analysis_out")
    parser_args.add_argument("--chemistry", default="chemistry.csv")
    parser_args.add_argument("--ecology", default="cluster_ecology.csv")
    parser_args.add_argument("--refs-fasta",
                             default="ROs_71_Clean/refs71.fasta")
    parser_args.add_argument("--refs-alignment",
                             default="ROs_71_Clean/refs71_hmmaln.sto")
    parser_args.add_argument("--alignment",
                             default="genomic_context/cand_aln.sto",
                             help="uye hizalamasi (varyant dagilimi icin)")
    parser_args.add_argument("--structures", default=STRUCTURE_CACHE_DIR,
                             help="yapi onbellegi dizini")
    parser_args.add_argument("--offline", action="store_true",
                             help="indirme yapma, sadece onbellegi kullan")
    args = parser_args.parse_args()

    try:
        from Bio.PDB import MMCIFParser
    except ImportError:
        print("[hata] biopython bulunamadi; Bio.PDB gerekli.", file=sys.stderr)
        return 2

    print("[okunuyor] kuratorlu tablolar")
    chemistry = {r["cluster"]: r
                 for r in csv.DictReader(open(args.chemistry))}
    ecology = {r["cluster"]: r
               for r in csv.DictReader(open(args.ecology))}
    references = read_fasta(args.refs_fasta)
    column_maps, matchcols, map_problems = build_column_maps(
        args.refs_alignment, references)

    # Kendi tam-hizalama okuyucumun ro_motif ile ayni match kolonlarini
    # verdigini SINA. Ayrisirsa butun kolon numaralari anlamsiz olur.
    library = read_stockholm_matchcols(args.refs_alignment)
    reader_agreement = (library == matchcols)
    if not reader_agreement:
        print("[hata] match kolonu cikarimi ro_motif ile uyusmuyor; durduruldu",
              file=sys.stderr)
        return 2
    print("[kontrol] match kolonu cikarimi ro_motif.read_stockholm_matchcols "
          "ile birebir ayni (%d kolon)" % len(next(iter(matchcols.values()))))
    for problem in map_problems:
        print("[uyari] %s" % problem)

    with_pdb = sorted((c, r["pdb"].strip().upper())
                      for c, r in chemistry.items() if r.get("pdb", "").strip())
    distinct_pdb = sorted({p for _, p in with_pdb})
    print("[bilgi] %d tipin %d'inde PDB kimligi var, %d ayri yapi"
          % (len(chemistry), len(with_pdb), len(distinct_pdb)))

    # --- yapilar ---
    os.makedirs(args.structures, exist_ok=True)
    fetch_status = {}
    for pdb_id in distinct_pdb:
        path, status = fetch_structure(pdb_id, args.structures, args.offline)
        fetch_status[pdb_id] = {"path": path, "status": status}
        print("[yapi] %s %s" % (pdb_id, status))

    aligner = make_aligner()
    cif_parser = MMCIFParser(QUIET=True)

    records, failures = [], []
    for cluster, pdb_id in with_pdb:
        reference_name, problem = reference_name_for_cluster(
            cluster, references.keys())
        if reference_name is None:
            failures.append({"type": cluster, "pdb_id": pdb_id,
                             "status": "undetermined", "problems": [problem]})
            continue
        if cluster not in column_maps.get(reference_name, {}) and \
                reference_name not in column_maps:
            failures.append({"type": cluster, "pdb_id": pdb_id,
                             "status": "undetermined",
                             "problems": ["no column map for reference %s"
                                          % reference_name]})
            continue
        entry = fetch_status.get(pdb_id, {})
        if not entry.get("path"):
            failures.append({"type": cluster, "pdb_id": pdb_id,
                             "status": "undetermined",
                             "problems": ["structure unavailable: %s"
                                          % entry.get("status", "unknown")]})
            continue
        record = analyse_structure(
            cluster, pdb_id, entry["path"], reference_name,
            references[reference_name], column_maps[reference_name],
            aligner, cif_parser)
        record["structure_source"] = entry["status"]
        if record["status"] == "undetermined":
            failures.append(record)
        else:
            records.append(record)
        print("[tip] %-14s %s  %s" % (cluster, pdb_id, record["status"]))

    # --- varyantlar ---
    variants = {}
    member_counts = {}
    sdp_overlap = {}
    if records:
        print("[okunuyor] uye hizalamasi %s" % args.alignment)
        aligned_members = read_stockholm_matchcols(args.alignment) \
            if os.path.exists(args.alignment) else {}
        print("[bilgi] hizalamada %d dizi" % len(aligned_members))
        connection = sqlite3.connect(args.db)
        confirmed = defaultdict(list)
        for candidate_id, cluster in connection.execute(
                "SELECT candidate_id, ro_cluster FROM ro WHERE is_confirmed=1"):
            confirmed[cluster].append(candidate_id)
        leaves_by_cluster = defaultdict(lambda: defaultdict(list))
        for candidate_id, cluster, leaf_id in connection.execute(
                "SELECT candidate_id, cluster, leaf_id FROM ro_leaf"):
            leaves_by_cluster[cluster][leaf_id].append(candidate_id)
        cluster_sdp = {}
        for cluster, columns in connection.execute(
                "SELECT cluster, columns FROM cluster_sdp"):
            try:
                cluster_sdp[cluster] = [m["column"] for m in json.loads(columns)]
            except (ValueError, KeyError, TypeError):
                continue
        connection.close()

        for record in records:
            cluster = record["type"]
            columns = record["pocket"]["%.1f" % PRIMARY_RADIUS_A][
                "alignment_columns"]
            members = [c for c in confirmed.get(cluster, [])
                       if c in aligned_members]
            member_counts[cluster] = {
                "confirmed_members_in_database": len(confirmed.get(cluster, [])),
                "members_present_in_alignment": len(members),
            }
            if not members:
                variants[cluster] = {
                    "status": "no_data",
                    "reason": ("this type has no confirmed member sequence in "
                               "the member alignment, so pocket variation "
                               "cannot be measured for it"),
                    "confirmed_members_in_database":
                        len(confirmed.get(cluster, [])),
                }
            else:
                sequences = [aligned_members[c] for c in members]
                per_column = {str(col): column_variation(col, sequences)
                              for col in columns}
                invariant = [c for c, v in per_column.items()
                             if v["invariant"] is True]
                variable = [c for c, v in per_column.items()
                            if v["invariant"] is False]
                leaves = {leaf_id: [c for c in ids if c in aligned_members]
                          for leaf_id, ids
                          in leaves_by_cluster.get(cluster, {}).items()}
                leaves = {k: v for k, v in leaves.items() if v}
                bins = {}
                for threshold in CONSERVATION_BINS:
                    bins["top_residue_fraction_at_least_%.2f" % threshold] = sum(
                        1 for v in per_column.values()
                        if v["top_fraction"] is not None
                        and v["top_fraction"] >= threshold)
                variants[cluster] = {
                    "status": "measured",
                    "source": ("aligned member sequences (match-state columns "
                               "of the motif profile)"),
                    "n_members": len(members),
                    "n_pocket_columns": len(columns),
                    "invariant_columns": sorted(int(c) for c in invariant),
                    "variable_columns": sorted(int(c) for c in variable),
                    "conservation_bins": bins,
                    "per_column": per_column,
                    "variant_level": variant_breakdown(columns, leaves,
                                                       aligned_members),
                }
            discriminating = set(cluster_sdp.get(cluster, []))
            sdp_overlap[cluster] = {
                "n_variant_discriminating_columns_from_cluster_sdp":
                    len(discriminating),
                "overlap_with_pocket_columns":
                    sorted(discriminating & set(columns)),
                "note": ("cluster_sdp columns were selected statistically by "
                         "variant_signature.py without any structural input; "
                         "an overlap means a statistically discriminating "
                         "position is also a structurally observed pocket "
                         "position"),
            }

    # --- kompozisyon ozetleri ---
    composition_rows = []
    for record in records:
        pocket = record["pocket"]["%.1f" % PRIMARY_RADIUS_A]
        composition_rows.append({
            "type": record["type"],
            "pdb_id": record["pdb_id"],
            "family": chemistry[record["type"]].get("family", ""),
            "reaction_class": chemistry[record["type"]].get("reaction_class", ""),
            "substrate_class": ecology.get(record["type"], {}).get(
                "substrate_class", ""),
            "substrate": chemistry[record["type"]].get("substrate_en", ""),
            "composition": pocket["composition"],
        })

    payload = {
        "module": "active_site.py",
        "generated_utc": datetime.now(timezone.utc).strftime(
            "%Y-%m-%dT%H:%M:%SZ"),
        "what_was_measured": (
            "For every curated reference type that carries a verified PDB "
            "identifier, the crystal structure was downloaded from RCSB, the "
            "catalytic mononuclear iron was told apart from the Rieske "
            "[2Fe-2S] irons on coordination-chemistry grounds, the residues "
            "lining the pocket around that iron were listed at two cut-off "
            "radii, each pocket residue was carried from PDB author numbering "
            "to the alignment match-state column of that type's reference "
            "sequence, and the residue distribution of the database members of "
            "that type was read out at every pocket column."),
        "parameters": {
            "pocket_radii_A": list(POCKET_RADII_A),
            "primary_radius_A": PRIMARY_RADIUS_A,
            "pocket_radius_justification": (
                "5.0 A is the first shell: the iron coordination sphere "
                "(about 2.2 A) plus one van der Waals contact, so it returns "
                "the metal ligands and the residues directly touching them. "
                "8.0 A is the substrate pocket: the substrates of these "
                "enzymes are extended molecules (naphthalene long axis about "
                "7 A, biphenyl about 9 A, steroid ring system longer still) "
                "and the far end of a bound substrate sits 6-9 A from the "
                "iron, so a smaller radius cannot enclose the residues that "
                "contact it. Both are reported so the reader can see how much "
                "of the answer depends on the radius instead of trusting one "
                "number."),
            "iron_ligand_distance_max_A": FE_LIGAND_MAX_A,
            "iron_iron_distance_max_A": FE_FE_MAX_A,
            "inorganic_sulfur_distance_max_A": INORGANIC_S_MAX_A,
            "iron_threshold_justification": (
                "Fe-N(His) is 2.0-2.2 A, Fe-O(carboxylate) 2.0-2.3 A and "
                "Fe-S(Cys) 2.2-2.3 A, so a 2.8 A ceiling admits all three "
                "bonds without reaching the second shell. Fe-Fe inside a "
                "[2Fe-2S] cluster measures 2.65-2.80 A, so a 3.5 A ceiling "
                "captures the pair while a mononuclear iron has no iron at "
                "all within it. Bridging Fe-S is 2.15-2.35 A, so 3.0 A "
                "tolerates crystallographic error."),
            "minimum_chain_identity_to_reference": MIN_CHAIN_IDENTITY,
            "minimum_chain_length": MIN_CHAIN_LENGTH,
            "minimum_members_per_variant": MIN_LEAF_MEMBERS,
            "alignment_scoring": {
                "substitution_matrix": "BLOSUM62",
                "gap_open": ALIGN_GAP_OPEN,
                "gap_extend": ALIGN_GAP_EXTEND,
                "end_gaps": "free",
            },
            "motif_columns_from_ro_motif": [
                {"column": c, "expected_residue": a, "role": role,
                 "label": label}
                for c, a, label, role in sorted(motif_columns())],
        },
        "method_notes": [
            "Structures are cached under the directory given by --structures "
            "(default 'structures/'), one request at a time with a %.1f s pause "
            "and a descriptive User-Agent, so a repeated run issues no network "
            "request at all." % DOWNLOAD_PAUSE_S,
            "Author numbering in a PDB entry is not the sequence index, so the "
            "structure's own observed sequence is aligned to the reference "
            "sequence and the column map is built from that alignment; nothing "
            "assumes the two numberings agree.",
            "The reference alignment ROs_71_Clean/refs71_hmmaln.sto was "
            "produced with hmmalign --trim, so the aligned sequence is a "
            "trimmed substring of the FASTA sequence. Counting aligned "
            "residues as FASTA indices would shift every column (for NahAc by "
            "4 residues, 449 against 433). The offset is located by searching "
            "for the trimmed sequence inside the FASTA sequence, and a "
            "reference whose offset cannot be located gets no column map.",
            "The match-state extraction used here was checked against "
            "ro_motif.read_stockholm_matchcols and agrees exactly; the module "
            "refuses to run if it does not.",
        ],
        "coverage": {
            "types_in_chemistry_csv": len(chemistry),
            "types_with_pdb_identifier": len(with_pdb),
            "distinct_pdb_entries": len(distinct_pdb),
            "structures_available": sum(1 for v in fetch_status.values()
                                        if v["path"]),
            "structures_downloaded": sum(1 for v in fetch_status.values()
                                         if v["status"] == "downloaded"),
            "structures_from_cache": sum(1 for v in fetch_status.values()
                                         if v["status"] == "cached"),
            "types_with_catalytic_iron_resolved": len(records),
            "types_with_catalytic_triad_verified_at_exact_columns": sum(
                1 for r in records
                if r["column_mapping_verification"][
                    "catalytic_triad_verified_at_exact_columns"]),
            "types_with_catalytic_triad_verified_within_window": sum(
                1 for r in records
                if r["column_mapping_verification"][
                    "catalytic_triad_verified_within_window"]),
            "types_with_a_column_shift": sorted(
                r["type"] for r in records
                if r["column_mapping_verification"]["column_shifts"]),
            "types_undetermined": len(failures),
            "types_with_member_variation_measured": sum(
                1 for v in variants.values() if v["status"] == "measured"),
            "independent_structure_sequence_pairs": len(
                {(r["pdb_id"], r["reference_sequence_sha1_12"])
                 for r in records}),
            "pseudo_replication_warning": (
                "The %d types resolve to only %d distinct PDB entries, because "
                "the curated reference set contains the same enzyme under more "
                "than one name (1NDO serves three types, 1Z03 and 1WW9 two "
                "each). Counts per type are therefore not independent "
                "observations; every aggregate below also reports the number "
                "of distinct (structure, reference sequence) pairs."
                % (len(records), len({r["pdb_id"] for r in records}))),
        },
        "structures": records,
        "undetermined": failures,
        "member_coverage": member_counts,
        "pocket_variation_across_members": variants,
        "overlap_with_variant_signatures": sdp_overlap,
        "electrostatics": {
            "what_this_is": "composition summary, NOT an electrostatic calculation",
            "explicitly_not_done": [
                "no Poisson-Boltzmann or any continuum electrostatics solver",
                "no pKa prediction or protonation-state assignment",
                "no electrostatic potential map, no field, no energy",
                "no solvent or counter-ion model",
            ],
            "why": (
                "A defensible electrostatic calculation needs specialist "
                "software, an explicit protonation assignment and a solvent "
                "model. What can be defended from the coordinates alone is the "
                "formal charge composition of the pocket at neutral pH, which "
                "is what is reported here. Histidine is counted separately "
                "rather than as a base, because the imidazole pKa of 6.0-7.0 "
                "leaves its protonation state undecided at neutral pH and the "
                "catalytic histidines are coordinated to iron anyway. Net "
                "charge is also given with the structurally verified iron "
                "ligands removed, because a carboxylate coordinated to iron "
                "should not be counted as a free negative charge."),
            "definitions": {
                "acidic": list(ACIDIC), "basic": list(BASIC),
                "histidine_counted_separately": True,
                "polar": list(POLAR), "hydrophobic": list(HYDROPHOBIC),
                "aromatic": list(AROMATIC),
                "net_formal_charge": "(K + R) - (D + E)",
            },
            "per_type": composition_rows,
            "what_the_composition_shows": None,
        },
        "limitations": [
            "Only %d of the 71 curated types carry a PDB identifier, and those "
            "%d types come from %d distinct structures, so nothing here "
            "generalises to the other types by itself."
            % (len(with_pdb), len(with_pdb), len(distinct_pdb)),
            "Pocket residues are read from one crystal conformation. Side "
            "chains move when a substrate binds, so a residue just outside the "
            "cut-off is not evidence that it never contacts substrate.",
            "The pocket is defined by distance to the iron, not by a bound "
            "substrate. Where a structure contains no substrate or inhibitor, "
            "the pocket is the cavity around the metal and not a measured "
            "binding site.",
            "Residue distributions at pocket columns describe the members "
            "assigned to that type by the profile library. Type assignment is "
            "a nearest-reference statement and about three quarters of members "
            "sit at distant or novel evidence tiers, so a varying pocket "
            "column may reflect that heterogeneity rather than a functional "
            "difference.",
            "Group comparisons of pocket composition have one to four types "
            "per group, after de-duplication even fewer. They are descriptive; "
            "no test is reported because none would be meaningful at that n.",
            "The acidic count of a pocket is dominated by residues that are "
            "conserved in the whole family rather than by anything "
            "substrate-specific: the iron carboxylate at column 355 is an "
            "aspartate in every resolved structure and column 354 carries an "
            "acidic or amide residue in most of them. Net pocket charge "
            "therefore mostly reports the conserved catalytic region, which is "
            "a further reason not to read it as a specificity signal.",
            "No binding energy, specificity prediction or substrate docking is "
            "attempted, and none can be supported by this data.",
        ],
    }
    if records:
        charges = [r["pocket"]["%.1f" % PRIMARY_RADIUS_A]["composition"][
            "net_formal_charge"] for r in records]
        aromatics = [r["pocket"]["%.1f" % PRIMARY_RADIUS_A]["composition"][
            "aromatic_F_W_Y"] for r in records]
        payload["electrostatics"]["what_the_composition_shows"] = {
            "net_formal_charge_range_over_resolved_types": [min(charges),
                                                            max(charges)],
            "every_pocket_is_net_negative": all(c < 0 for c in charges),
            "aromatic_count_range": [min(aromatics), max(aromatics)],
            "honest_reading": (
                "Every resolved pocket carries a net negative formal charge, "
                "and the spread is narrow (a few charge units over %d types). "
                "The group means by substrate class differ in the same "
                "direction one would guess from the chemistry, but the groups "
                "hold one to twelve types and only %d independent (structure, "
                "sequence) pairs in total, their ranges overlap, and the "
                "acidic count is largely set by the conserved carboxylates "
                "next to the iron. So the composition summary does not "
                "establish that pockets of different substrate classes differ "
                "electrostatically; it establishes that all of them are "
                "anionic and that this data cannot resolve finer differences."
                % (len(records),
                   len({(r["pdb_id"], r["reference_sequence_sha1_12"])
                        for r in records}))),
        }
        payload["composition_by_family"] = group_composition(
            records, lambda r: chemistry[r["type"]].get("family", ""), "family")
        payload["composition_by_substrate_class"] = group_composition(
            records,
            lambda r: ecology.get(r["type"], {}).get("substrate_class", ""),
            "substrate_class")
        payload["composition_by_reaction_class"] = group_composition(
            records,
            lambda r: chemistry[r["type"]].get("reaction_class", ""),
            "reaction_class")

        # Kolon bazinda birlesik gorunum: hangi kolon kac yapida cebi doseyor.
        for radius in POCKET_RADII_A:
            tally = defaultdict(lambda: {"types": [], "independent": set(),
                                         "residues": {}})
            for record in records:
                pocket = record["pocket"]["%.1f" % radius]
                for row in pocket["residues"]:
                    column = row["alignment_column"]
                    if column is None:
                        continue
                    bucket = tally[column]
                    bucket["types"].append(record["type"])
                    bucket["independent"].add(
                        (record["pdb_id"],
                         record["reference_sequence_sha1_12"]))
                    bucket["residues"][record["type"]] = row["residue_1letter"]
            ligand_columns = {c for c, _, _, _ in motif_columns()}
            payload.setdefault("pocket_columns", {})["%.1f" % radius] = {
                "radius_A": radius,
                "n_columns": len(tally),
                "columns": [{
                    "column": column,
                    "n_types": len(data["types"]),
                    "n_independent_structures": len(data["independent"]),
                    "types": sorted(data["types"]),
                    "residue_per_type": dict(sorted(data["residues"].items())),
                    "is_motif_model_column": column in ligand_columns,
                } for column, data in sorted(tally.items())],
                "columns_present_in_every_resolved_type": sorted(
                    c for c, d in tally.items() if len(d["types"]) == len(records)),
            }

    os.makedirs(args.out_dir, exist_ok=True)
    out_path = os.path.join(args.out_dir, "active_site.json")
    with open(out_path, "w") as handle:
        json.dump(payload, handle, indent=2, sort_keys=False)

    # ----------------------------------------------------------- stdout ozeti
    print()
    print("=" * 78)
    print("AKTIF BOLGE -- kapsam")
    print("=" * 78)
    cover = payload["coverage"]
    print("  PDB kimligi olan tip        : %d / %d"
          % (cover["types_with_pdb_identifier"], cover["types_in_chemistry_csv"]))
    print("  ayri yapi                   : %d" % cover["distinct_pdb_entries"])
    print("  indirilen / onbellekten     : %d / %d"
          % (cover["structures_downloaded"], cover["structures_from_cache"]))
    print("  katalitik demir cozuldu     : %d tip"
          % cover["types_with_catalytic_iron_resolved"])
    print("  triad TAM kolonda dogrulandi: %d tip"
          % cover["types_with_catalytic_triad_verified_at_exact_columns"])
    print("  triad +-%d pencerede dogrulandi: %d tip"
          % (COLUMN_WINDOW,
             cover["types_with_catalytic_triad_verified_within_window"]))
    if cover["types_with_a_column_shift"]:
        print("  kolon kaymasi olan tip      : %s"
              % ", ".join(cover["types_with_a_column_shift"]))
    print("  belirlenemedi               : %d tip" % cover["types_undetermined"])
    print("  bagimsiz (yapi,dizi) cifti  : %d"
          % cover["independent_structure_sequence_pairs"])
    for failure in failures:
        print("  [belirlenemedi] %s %s: %s"
              % (failure["type"], failure.get("pdb_id", "?"),
                 "; ".join(failure.get("problems", []))))

    print()
    print("=" * 78)
    print("TIP BASINA CEB (%s A birincil yaricap)" % PRIMARY_RADIUS_A)
    print("=" * 78)
    print("  %-14s %-5s %4s %4s %5s %6s %5s  %s"
          % ("tip", "pdb", "5A", "8A", "kolon", "netyuk", "arom", "triad"))
    for record in records:
        inner = record["pocket"]["5.0"]
        outer = record["pocket"]["%.1f" % PRIMARY_RADIUS_A]
        composition = outer["composition"]
        print("  %-14s %-5s %4d %4d %5d %6d %5d  %s"
              % (record["type"], record["pdb_id"],
                 inner["n_residues_same_chain"], outer["n_residues_same_chain"],
                 len(outer["alignment_columns"]),
                 composition["net_formal_charge"],
                 composition["aromatic_F_W_Y"],
                 ("OK" if record["column_mapping_verification"][
                     "catalytic_triad_verified_at_exact_columns"]
                  else "KAYMALI-OK" if record["column_mapping_verification"][
                     "catalytic_triad_verified_within_window"]
                  else "DOGRULANMADI")))

    print()
    print("=" * 78)
    print("DEMIR AYRIMI")
    print("=" * 78)
    for record in records:
        inventory = record["iron_inventory"]
        print("  %-14s %s: %d demir = %d Rieske (%d cift) + %d katalitik"
              % (record["type"], record["pdb_id"], inventory["n_iron_atoms"],
                 inventory["n_rieske_cluster_irons"],
                 len(inventory["rieske_pairs"]),
                 inventory["n_catalytic_mononuclear_irons"]))

    print()
    print("=" * 78)
    print("CAPRAZ KONTROL -- merkezler arasi mesafe ve kopru karboksilati")
    print("=" * 78)
    for record in records:
        separation = record.get("centre_separation", {})
        bridge = record.get("bridging_carboxylate", {})
        print("  %-14s en yakin Rieske Fe %5s A (zincir %s, ayni zincir %5s A)"
              "  kopru %s d(Fe)=%s A ikisini birden=%s"
              % (record["type"],
                 separation.get("nearest_rieske_iron_A"),
                 separation.get("nearest_rieske_iron_chain"),
                 separation.get("nearest_rieske_iron_within_same_chain_A"),
                 bridge.get("structure_residue"),
                 bridge.get("distance_from_catalytic_iron_A"),
                 bridge.get("bridges_both_centres")))

    if "pocket_columns" in payload:
        print()
        print("=" * 78)
        print("CEB KOLONLARI (%s A)" % PRIMARY_RADIUS_A)
        print("=" * 78)
        block = payload["pocket_columns"]["%.1f" % PRIMARY_RADIUS_A]
        print("  toplam %d kolon; %d kolon cozulen TUM tiplerde var"
              % (block["n_columns"],
                 len(block["columns_present_in_every_resolved_type"])))
        print("  her tipte olanlar: %s"
              % ", ".join(str(c) for c in
                          block["columns_present_in_every_resolved_type"]))

    measured = [(c, v) for c, v in variants.items() if v["status"] == "measured"]
    if measured:
        print()
        print("=" * 78)
        print("CEB KOLONLARINDA VARYANT DEGISIMI")
        print("=" * 78)
        print("  %-14s %6s %6s %9s %9s"
              % ("tip", "uye", "kolon", "degismez", "degisken"))
        for cluster, data in sorted(measured):
            print("  %-14s %6d %6d %9d %9d"
                  % (cluster, data["n_members"], data["n_pocket_columns"],
                     len(data["invariant_columns"]),
                     len(data["variable_columns"])))
        print()
        for cluster, data in sorted(measured):
            overlap = sdp_overlap.get(cluster, {}).get(
                "overlap_with_pocket_columns", [])
            if overlap:
                print("  [kesisim] %s: varyant imza kolonlari ile ortak ceb "
                      "kolonu %s" % (cluster, overlap))
    no_data = [c for c, v in variants.items() if v["status"] == "no_data"]
    if no_data:
        print("  [veri yok] %s -- bu tiplerin hizalamada uyesi yok"
              % ", ".join(sorted(no_data)))

    print()
    print("[not] JSON icindeki 'electrostatics' bolumu bir KOMPOZISYON "
          "ozetidir, elektrostatik hesap degildir.")
    print("[yazildi] %s" % out_path)
    return 0


if __name__ == "__main__":
    sys.exit(main())
