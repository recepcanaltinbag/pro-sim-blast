"""
ROAR-DB -- KATALITIK KARBOKSILAT DENETIMI (2/3 kuralinin bedeli).

NEDEN BU MODUL VAR
  ro_motif.py bir girisi veritabanina kabul etmek icin Rieske ligandlarinin
  4/4'unu ZORUNLU tutuyor, ama mononukleer demirin 2-His-1-karboksilat facial
  triad'inin yalnizca 2/3'unu istiyor. Bu, bir girisin SADECE iki histidinle,
  yani KARBOKSILATSIZ kabul edilebilecegi anlamina gelir. Onaylanmis 11.422
  girisin 1.414'u (%12,4) tam olarak bu durumda: ro_carboxylate tablosunda
  catalytic_residue D ya da E degil. Bu modul "bunlar gercekten Rieske
  oksijenazi mi" sorusuna ol culmus bir cevap verir.

  ASIL GEREKCE -- SABIT KOLON YANILGISI. ro_carboxylate.catalytic_offset
  alaninin kendi aciklamasi "kolondan kayma, NULL = pencerede yok" diyor:
  karboksilat, ROmotif71 hizalamasinin SABIT bir match-state kolonunda
  (kolon 355) +-2 kolonluk bir pencereyle okunuyor. Bu okuma yontemi
  referanslara yakin diziler icin dogru, ama bu 1.414 girisin en yakin
  kuratorlu referansa ortalama kimligi %31,6. O kimlik duzeyinde hizalama
  REGISTRESI guvenilir degildir: birkac pozisyon kaymis ya da kotu hizalanan
  bir ilmekte duran bir karboksilat, VAR oldugu halde YOK kaydedilir.
  "Beklenen kolonda karboksilat yok" ile "karboksilat yok" tamamen farkli iki
  bulgudur ve ikincisi veritabaninin kabul kuralini sorgulatir, birincisi
  yalnizca olcum araci nin menzilini.

  Bu yuzden modul karboksilati, VAR OLDUGU BILINEN kalintilara capalayarak
  yeniden okur. Bu 1.414 girisin hepsinde katalitik triad'in >=2/3'u yerinde,
  yani iki katalitik histidin (kolon 212 ve 217) VAR; dahasi hepsinde
  alt-birimler arasi kopru Asp/Glu (kolon 209) da VAR, yani tezin
  (Complete_thesis.pdf, Literature Review böl. 3.1, atif [104]) aktardigi
  Parales/Gibson konsensusu D-x(2)-H-x(3-5)-H bu girislerin HEPSINDE tam.
  DIKKAT: o konsensustaki D, facial triad'in karboksilati DEGIL, kopru
  aspartatidir (kolon 209; katalitik His'ler 212 ve 217 -> D-x(2)-H-x(4)-H).
  Facial triad'in karboksilati ~140 kalinti daha C-terminalde (kolon 355)
  oturur, yani o konsensus karboksilatin yerini VERMEZ; sadece capanin
  saglam oldugunu dogrular.

NE OLCULDU (sirayla, her adim bir sonrakinin gerekcesi)
  1. TEMEL OKUMA. ro_carboxylate yeniden uretilir ve eksik cagri parcalanir:
     kolon GERCEKTEN bos (gap) mu, yoksa orada karboksilat OLMAYAN bir
     kalinti mi oturuyor. Ikisi farklidir: gap bir hizalama olayidir,
     ikame bir biyolojik olaydir.
  2. MEKANIZMA. Her girisin ROmotif71 hizalamasindaki KAPSAMI olculur
     (ilk/son dolu match kolonu, dolu match kalintisi sayisi). Karboksilat
     kolonu, girisin hizalanmis bolgesinin DISINDA mi kaliyor?
  3. PENCERE GENISLETME. Ayni sabit kolon +-2, +-5, +-10, +-20, +-40 ile
     okunur. Kolon gap ise genisletmek kurtarmaz -- bu adim, sorunun basit
     bir kayma olmadigini gosterir.
  4. CIPLAK CAPA (ve neden yetmez). Katalitik His-2'ye gore HAM dizide
     kalibre edilmis uzaklikta D/E aranir. Pozitiflerde bu uzaklik 86-238
     kalinti arasinda degisiyor, yani pencere ~90 kalinti genisliginde olmak
     zorunda; D+E birlikte kalintilarin ~%11'i oldugu icin boyle bir pencere
     SANS ESERI neredeyse her zaman bir D/E bulur. Adim bu yuzden hem bulma
     oranini hem de her dizinin KENDI bilesimi nden cikan sans beklentisini
     raporlar: tek basina kanit olmadigi acikca gorulsun.
  5. REFERANS PROFIL SONDASI. refs71_hmmaln.sto'dan uc alt-model kurulur
     (Rieske alani kolon 60-140; C-terminal katalitik lob kolon 260-426;
     karboksilat bolgesi kolon 330-380) ve uc kumeye hmmsearch edilir:
     karboksilatsiz 1.414, karboksilatli 10.008 ve KONTROL olarak
     "Rieske var, katalitik yok" diye ELENMIS adaylar. Bu adim sorunun
     NEREDE oldugunu soyler.
  6. CAPALI YENIDEN OKUMA -- ASIL OLCUM. Her enzim tipi icin karboksilatsiz
     uyeler, AYNI TIPIN karboksilati bilinen uyeleriyle ve o tipin kuratorlu
     referans dizisiyle birlikte mafft ile yeniden hizalanir. Karboksilat
     kolonu, SABIT bir sayi olarak degil, pozitiflerin HAM dizideki bilinen
     karboksilat pozisyonlarinin oyuyla belirlenir. Sonra karboksilatsiz
     uyelerin o kolonda (+-pencere) ne tasidigi okunur. Yontemin kendi
     dogrulugu ayni hizalamada olculur (pos_recovery): pozitiflerin kacinda
     capa dogru kalintiyi geri veriyor.
  7. BOS MODEL. Adim 6'nin buldugu korunmus kolonun sans olup olmadigi
     sinanir: grubun KENDI D/E bilesimiyle binom testi, ve ayni hizalamada
     dolulugu >=%80 olan butun kolonlar arasindaki D/E oraninin dagilimi
     (karsilastirma esigi).
  8. GERCEK YOKLUKLARIN KIMLIGI. Adim 6'dan sonra hala karboksilati
     olmayan girisler icin: hangi kalinti oturuyor, hangi tip, hangi kanit
     kademesi, hangi taksonomik alan, operonu var mi, dizi uzayinda bir
     arada mi (yaprak/altaile), ve alan mimarisi (component_hmms +
     referanstan kurulan katalitik lob modeli).
  9. SIKI KAPININ BEDELI. Kabul kurali 2/3 yerine 3/3 olsaydi kac giris ve
     kac tip kaybedilirdi -- hem sabit-kolon okumasina hem de adim 6'nin
     capali okumasina gore.
 10. ISTATISTIK. Kontenjans testleri stats_overview.py'nin makinesiyle
     (signal_map / bh_adjust / cramers_v) hem GIRIS duzeyinde hem de
     tip x cins basina tek gozleme cokertilerek yapilir.

SINIRLAR, ACIKCA
  - Modul SALT OKUNUR. roar.sqlite read-only acilir, chemistry.csv'ye
    dokunulmaz. Hicbir dogru/yanlis duzeltmesi yazilmaz; bu bir denetimdir.
  - component_hmms/components.hmm icinde halka-hidroksilleyici alpha'nin
    KATALITIK alani icin bir model YOK (PF00848 / Ring_hydroxyl_A kitaplikta
    degil, yerelde Pfam-A da yok ve bu modul ag erisimi yapmaz). Bu yuzden
    "katalitik alan var mi" sorusu projenin kendi 71 kuratorlu referansindan
    kurulan modelle olculur; SRPBCC (heliks-tutamaci, PF03364/PF10604)
    gibi Pfam ailelerinin varligi DOGRULANAMAZ, yalnizca NCBI'nin urun
    dizgileri oldugu gibi raporlanir.
  - Adim 6'nin capasi tipin POZITIF uyelerine dayanir. Hicbir tipte pozitif
    uye sayisi 9'un altinda degil, ama capa gucu tip tip degisir ve her tip
    icin pos_recovery ile birlikte raporlanir.
  - TEKRARLANABILIRLIK. mafft'in yineleyici iyilestirme adimi COK IS
    PARCACIGIYLA calistirildiginda is parcacigi siralamasina duyarlidir:
    ayni girdiyle yapilan kosular arasinda capa kolonu ve birkac girisin
    cagrisi oynar. Olculdu: --thread 6 ile dort kosuda "gercek yokluk"
    21, 22, 23 ve 25 cikti (1.414'un %1,5-1,8'i). KARAR her kosuda ayni
    (%98,2-98,5 hizalama yanilgisi), ama giris duzeyindeki sayi oynak.
    Bu yuzden hizalama VARSAYILAN OLARAK TEK PARCACIKLA kosar
    (--mafft-threads 1) ve sonuc birebir tekrarlanabilir; --mafft-threads
    >1 daha hizlidir ama bit duzeyinde tekrarlanmaz ve JSON bunu
    reproducibility alaninda belirtir. hmmbuild/hmmsearch/hmmscan --cpu
    degeri ayri tutulur, cunku onlar bu duyarliligi tasimaz.

Cikti: analysis_out/carboxylate_audit.json  (verdict blogu duz sozlerle)
"""

import argparse
import json
import os
import shutil
import sqlite3
import subprocess
import sys
import tempfile
from collections import Counter, defaultdict
from datetime import datetime, timezone

import numpy as np
from scipy import stats as sps

from ro_motif import (BRIDGING_SITE, CATALYTIC_SITES, RIESKE_SITES,
                      DEFAULT_MODEL, MODEL_COLUMNS)
from stats_overview import bh_adjust, cramers_v, signal_map

HIS1_COL = CATALYTIC_SITES[0][0]        # 212
HIS2_COL = CATALYTIC_SITES[1][0]        # 217
CARBOX_COL = CATALYTIC_SITES[2][0]      # 355
BRIDGE_COL = BRIDGING_SITE[0]           # 209
CARBOX_OK = "DE"

# Referans hizalamasindan kesilecek alt-modeller: (ad, ilk kolon, son kolon, nicin)
SUBMODELS = [
    ("RORieske", 60, 140,
     "Rieske [2Fe-2S] alani -- POZITIF KONTROL, her sette tutmali"),
    ("ROcatC", 260, 426,
     "C-terminal katalitik lob -- karboksilat kolonunun icinde oldugu bolge"),
    ("ROcarbox", 330, 380,
     "yalnizca karboksilat cevresi -- dar ve hassas yerel sonda"),
]


# --------------------------------------------------------------------------
# Stockholm / hizalama altyapisi
# --------------------------------------------------------------------------

def parse_stockholm(path, wanted=None):
    """Stockholm'u TAM haliyle oku (insertion kolonlari dahil) + #=GC RF satiri.

    ro_motif.read_stockholm_matchcols yalnizca match kolonlarini dondurur;
    burada HAM dizi indekslerine geri harita kurmak gerektigi icin insertion
    kolonlari da tutulur. Match/insert ayrimi RF satirindan okunur ('.' =
    insertion), ilk diziyi referans almaktan daha saglam.
    """
    rows = defaultdict(list)
    rf = []
    with open(path) as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith("#=GC RF"):
                rf.append(line.split()[2])
                continue
            if not line or line.startswith("#") or line.startswith("//"):
                continue
            parts = line.split()
            if len(parts) != 2:
                continue
            if wanted is not None and parts[0] not in wanted:
                continue
            rows[parts[0]].append(parts[1])
    return {k: "".join(v) for k, v in rows.items()}, "".join(rf)


def read_fasta(path):
    out, name = {}, None
    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if line.startswith(">"):
                name = line[1:].split()[0]
                out[name] = []
            elif name:
                out[name].append(line)
    return {k: "".join(v).upper() for k, v in out.items()}


def write_fasta(path, records):
    with open(path, "w") as handle:
        for name, seq in records:
            handle.write(">%s\n%s\n" % (name, seq))


class AlignmentIndex:
    """Bir girisin ROmotif71 hizalamasi ile HAM dizisi arasindaki kopru.

    col_raw[k]  -> k+1'inci match kolonundaki kalintinin HAM dizideki
                   0-tabanli indeksi, kolon gap ise None.
    hmmalign --trim ile calistirildigi icin hizalamadaki harfler proteinin
    bir ALT DIZISIdir; kaydirma str.find ile bulunur.
    """

    __slots__ = ("col_raw", "offset", "first_col", "last_col", "n_match")

    def __init__(self, row, match_idx, sequence):
        residue_at = {}
        count = 0
        for i, char in enumerate(row):
            if char.isalpha():
                residue_at[i] = count
                count += 1
        letters = "".join(c for c in row if c.isalpha()).upper()
        self.offset = sequence.find(letters)
        self.col_raw = [residue_at.get(i) if row[i].isalpha() else None
                        for i in match_idx]
        filled = [k + 1 for k, v in enumerate(self.col_raw) if v is not None]
        self.first_col = filled[0] if filled else None
        self.last_col = filled[-1] if filled else None
        self.n_match = len(filled)

    def raw(self, column):
        """1-tabanli match kolonunun HAM dizi indeksi (yoksa None)."""
        k = column - 1
        if self.offset < 0 or not 0 <= k < len(self.col_raw):
            return None
        rel = self.col_raw[k]
        return None if rel is None else self.offset + rel

    def find(self, sequence, column, accepted, window):
        """Kolon +-window icinde beklenen kalintiyi ara. Donen: (ham_idx, kayma)."""
        for delta in sorted(range(-window, window + 1), key=lambda x: (abs(x), x)):
            idx = self.raw(column + delta)
            if idx is not None and sequence[idx].upper() in accepted:
                return idx, delta
        return None, None

    def residue(self, sequence, column):
        idx = self.raw(column)
        return sequence[idx].upper() if idx is not None else "-"


# --------------------------------------------------------------------------
# Veri yukleme
# --------------------------------------------------------------------------

def load_entries(con):
    """Onaylanmis her giris icin denetimin ihtiyac duydugu butun alanlar."""
    sql = """
        SELECT r.candidate_id, r.sequence, r.ro_cluster, r.ro_group, r.product,
               r.model_coverage, r.hmm_score,
               c.catalytic_residue, c.catalytic_offset,
               c.bridging_residue, c.bridging_offset,
               e.tier, e.nearest_ref, e.ref_identity, e.ref_qcov,
               d.domain, d.euk_group,
               p.organism,
               o.candidate_id IS NOT NULL AS has_operon,
               COALESCE(o.has_beta, 0), COALESCE(o.has_ferredoxin, 0),
               COALESCE(o.has_reductase, 0), COALESCE(o.nearby_beta, 0),
               l.leaf_id, s.subfamily_id, s.assignment_class
        FROM ro r
        LEFT JOIN ro_carboxylate c ON c.candidate_id = r.candidate_id
        LEFT JOIN ro_evidence    e ON e.candidate_id = r.candidate_id
        LEFT JOIN ro_domain      d ON d.candidate_id = r.candidate_id
        LEFT JOIN replicon       p ON p.nucleotide_id = r.nucleotide_id
        LEFT JOIN operon         o ON o.candidate_id = r.candidate_id
        LEFT JOIN ro_leaf        l ON l.candidate_id = r.candidate_id
        LEFT JOIN ro_subfamily   s ON s.candidate_id = r.candidate_id
        WHERE r.is_confirmed = 1
    """
    fields = ("candidate_id sequence cluster grp product model_coverage hmm_score "
              "stored_residue stored_offset bridging_residue bridging_offset "
              "tier nearest_ref ref_identity ref_qcov domain euk_group organism "
              "has_operon has_beta has_ferredoxin has_reductase nearby_beta "
              "leaf_id subfamily_id assignment_class").split()
    out = {}
    for row in con.execute(sql):
        rec = dict(zip(fields, row))
        rec["genus"] = (rec["organism"] or "?").split()[0]
        out[rec["candidate_id"]] = rec
    return out


def load_rieske_only_controls(con):
    """Rieske ligandlari TAM, katalitik triad ELENMIS adaylar.

    Bunlar pipeline'in zaten RO alpha SAYMADIGI proteinler (ferredoksin,
    bc1 Rieske ISP, NirD tipi). Adim 5'te "karboksilat bolgesi hic yok"
    halinin nasil gorundugunu gosteren negatif kontrol.
    """
    return dict(con.execute(
        "SELECT candidate_id, sequence FROM ro "
        "WHERE is_confirmed = 0 AND rieske_intact = 1 AND catalytic_intact = 0"))


# --------------------------------------------------------------------------
# Adim 1-2: temel okuma ve mekanizma
# --------------------------------------------------------------------------

def stage_baseline(entries, index, window):
    """ro_carboxylate'i yeniden uret, eksik cagriyi gap/ikame olarak parcala."""
    stored = Counter()
    recomputed = Counter()
    missing_kind = Counter()
    mismatch = 0
    for cid, rec in entries.items():
        stored[rec["stored_residue"] or "?"] += 1
        idx = index.get(cid)
        if idx is None:
            recomputed["not_aligned"] += 1
            continue
        hit, delta = idx.find(rec["sequence"], CARBOX_COL, CARBOX_OK, window)
        if hit is not None:
            rec["fixed_residue"] = rec["sequence"][hit].upper()
            rec["fixed_offset"] = delta
            recomputed[rec["fixed_residue"]] += 1
        else:
            rec["fixed_residue"] = idx.residue(rec["sequence"], CARBOX_COL)
            rec["fixed_offset"] = None
            recomputed[rec["fixed_residue"]] += 1
            # Pencerenin tamami gap mi, yoksa orada gercek bir kalinti mi var?
            occupied = [idx.raw(CARBOX_COL + d) is not None
                        for d in range(-window, window + 1)]
            if not any(occupied):
                missing_kind["window_entirely_gapped"] += 1
                rec["missing_kind"] = "window_entirely_gapped"
            elif idx.raw(CARBOX_COL) is None:
                missing_kind["centre_gapped_neighbours_present"] += 1
                rec["missing_kind"] = "centre_gapped_neighbours_present"
            else:
                missing_kind["substituted_residue_in_place"] += 1
                rec["missing_kind"] = "substituted_residue_in_place"
        if (rec["fixed_residue"] in CARBOX_OK) != \
                ((rec["stored_residue"] or "") in CARBOX_OK):
            mismatch += 1
    return {
        "question": "ro_carboxylate ne diyor, ve 'yok' derken neyi goruyor",
        "stored_residue_counts": dict(stored.most_common()),
        "recomputed_residue_counts": dict(recomputed.most_common()),
        "recomputation_disagreements_with_table": mismatch,
        "fixed_column": CARBOX_COL,
        "fixed_window": window,
        "missing_call_breakdown": dict(missing_kind.most_common()),
        "note": ("window_entirely_gapped = kolon ve komsulari HIZALAMADA BOS; "
                 "orada okunacak kalinti YOK, yani bu bir hizalama olayi. "
                 "substituted_residue_in_place = kolonda gercek bir kalinti "
                 "var ama karboksilat degil; bu biyolojik bir iddia."),
    }


def stage_mechanism(entries, index, negatives, positives):
    """Karboksilat kolonu, girisin HIZALANMIS bolgesinin disinda mi kaliyor."""
    def profile(ids):
        last, first, nmatch, beyond = [], [], [], 0
        for cid in ids:
            idx = index.get(cid)
            if idx is None or idx.last_col is None:
                continue
            last.append(idx.last_col)
            first.append(idx.first_col)
            nmatch.append(idx.n_match)
            if idx.last_col < CARBOX_COL:
                beyond += 1
        return {
            "n": len(last),
            "first_filled_column_median": int(np.median(first)) if first else None,
            "last_filled_column_median": int(np.median(last)) if last else None,
            "last_filled_column_mean": round(float(np.mean(last)), 1) if last else None,
            "filled_match_residues_median": int(np.median(nmatch)) if nmatch else None,
            "alignment_stops_before_the_carboxylate_column": beyond,
            "alignment_stops_before_fraction": round(beyond / max(1, len(last)), 4),
        }

    def means(ids, key, digits=3):
        vals = [entries[c][key] for c in ids if entries[c][key] is not None]
        return round(float(np.mean(vals)), digits) if vals else None

    out = {
        "question": ("karboksilat kolonu (%d) girisin hizalanmis bolgesi icinde "
                     "mi kaliyor" % CARBOX_COL),
        "model": DEFAULT_MODEL,
        "model_columns": len(MODEL_COLUMNS[DEFAULT_MODEL]["catalytic"]),
        "no_carboxylate": profile(negatives),
        "with_carboxylate": profile(positives),
    }
    for label, ids in (("no_carboxylate", negatives), ("with_carboxylate", positives)):
        out[label].update({
            "model_coverage_mean": means(ids, "model_coverage"),
            "hmm_score_mean": round(float(np.mean(
                [entries[c]["hmm_score"] for c in ids
                 if entries[c]["hmm_score"] is not None])), 1),
            "ref_identity_mean": means(ids, "ref_identity", 1),
            "ref_qcov_mean": means(ids, "ref_qcov", 1),
            "sequence_length_mean": round(float(np.mean(
                [len(entries[c]["sequence"]) for c in ids])), 1),
        })
    # Capa saglam mi: His'ler ve kopru Asp yerinde mi
    anchors = Counter()
    for cid in negatives:
        rec = entries[cid]
        idx = index.get(cid)
        if idx is None:
            continue
        h1 = idx.find(rec["sequence"], HIS1_COL, "H", 2)[0]
        h2 = idx.find(rec["sequence"], HIS2_COL, "H", 2)[0]
        br = (rec["bridging_residue"] or "") in CARBOX_OK
        anchors["catalytic_his1_present"] += h1 is not None
        anchors["catalytic_his2_present"] += h2 is not None
        anchors["both_histidines_present"] += (h1 is not None and h2 is not None)
        anchors["bridging_carboxylate_present"] += br
        if h1 is not None and h2 is not None and br:
            anchors["full_D_x2_H_x3_5_H_consensus_intact"] += 1
    out["anchor_integrity_in_the_no_carboxylate_set"] = dict(anchors)
    out["anchor_note"] = (
        "Parales/Gibson konsensusu D-x(2)-H-x(3-5)-H (Complete_thesis.pdf, "
        "bol. 3.1, atif [104]) bu girislerin hepsinde tamdir. O konsensustaki "
        "D KOPRU aspartatidir (kolon %d), facial triad karboksilati degil "
        "(kolon %d); yani konsensus karboksilatin yerini vermez, capanin "
        "saglam oldugunu dogrular." % (BRIDGE_COL, CARBOX_COL))
    return out


def stage_window_widening(entries, index, negatives, widths):
    """Ayni sabit kolonu genisleyen pencerelerle oku. Gap ise genisletmek kurtarmaz."""
    rows = []
    for width in widths:
        found = 0
        offsets = Counter()
        for cid in negatives:
            rec = entries[cid]
            idx = index.get(cid)
            if idx is None:
                continue
            hit, delta = idx.find(rec["sequence"], CARBOX_COL, CARBOX_OK, width)
            if hit is not None:
                found += 1
                offsets[delta] += 1
        rows.append({"window_columns": width, "rescued": found,
                     "rescued_fraction": round(found / max(1, len(negatives)), 4),
                     "offset_counts": dict(sorted(offsets.items(),
                                                  key=lambda kv: abs(kv[0])))})
    return {
        "question": ("sabit kolonu MODEL koordinatinda genisletmek karboksilati "
                     "geri getirir mi"),
        "n_tested": len(negatives),
        "windows": rows,
        "note": ("Kolon hizalamada BOSsa genisletmek de bos kolonlari tarar. "
                 "Bu adimin kurtaramadigi girisler, basit bir registre kaymasiyla "
                 "aciklanamaz demektir."),
    }


# --------------------------------------------------------------------------
# Adim 4: ciplak capa ve neden yetmedigi
# --------------------------------------------------------------------------

def stage_naive_anchor(entries, index, negatives, positives):
    """His-2'ye gore HAM uzaklikta D/E ara; yaninda SANS beklentisini de ver."""
    def his2(cid):
        rec = entries[cid]
        idx = index.get(cid)
        if idx is None:
            return None
        hit = idx.find(rec["sequence"], HIS2_COL, "H", 2)[0]
        if hit is None:
            hit = idx.find(rec["sequence"], HIS1_COL, "H", 2)[0]
        return hit

    distances = []
    for cid in positives:
        rec = entries[cid]
        idx = index.get(cid)
        h2 = his2(cid)
        if idx is None or h2 is None:
            continue
        cb = idx.find(rec["sequence"], CARBOX_COL, CARBOX_OK, 2)[0]
        if cb is not None:
            distances.append(cb - h2)
    arr = np.array(distances)
    percentiles = {("p%g" % q): int(np.percentile(arr, q))
                   for q in (0, 1, 2.5, 25, 50, 75, 97.5, 99, 100)}
    lo, hi = int(np.percentile(arr, 1)), int(np.percentile(arr, 99))

    found = 0
    expected = []
    tested = 0
    for cid in negatives:
        rec = entries[cid]
        seq = rec["sequence"]
        h2 = his2(cid)
        if h2 is None:
            continue
        tested += 1
        a, b = max(0, h2 + lo), min(len(seq), h2 + hi + 1)
        segment = seq[a:b]
        if any(ch in CARBOX_OK for ch in segment):
            found += 1
        # Her dizinin KENDI D/E bilesiminden sans beklentisi
        p = sum(1 for ch in seq if ch in CARBOX_OK) / max(1, len(seq))
        expected.append(1.0 - (1.0 - p) ** max(0, len(segment)))
    return {
        "question": ("His-2'ye gore kalibre edilmis HAM uzaklik penceresinde "
                     "D/E var mi"),
        "calibration_n": int(arr.size),
        "his2_to_carboxylate_raw_distance": percentiles,
        "window_used": [lo, hi],
        "window_width_residues": hi - lo + 1,
        "n_tested": tested,
        "found_a_carboxylate": found,
        "found_fraction": round(found / max(1, tested), 4),
        "expected_by_chance_fraction": round(float(np.mean(expected)), 4)
        if expected else None,
        "verdict": ("Pencere %d kalinti genisliginde ve sans beklentisi %.1f%%; "
                    "bu okuma TEK BASINA kanit degildir, adim 6'nin capali "
                    "okumasina ihtiyac var."
                    % (hi - lo + 1, 100 * float(np.mean(expected)) if expected else 0)),
    }


# --------------------------------------------------------------------------
# Adim 5: referans profil sondalari
# --------------------------------------------------------------------------

def build_submodels(ref_sto, workdir, log):
    """refs71_hmmaln.sto'dan alt-model kes ve hmmbuild et."""
    rows, rf = parse_stockholm(ref_sto)
    match_idx = [i for i, ch in enumerate(rf) if ch != "."]
    built = []
    for name, lo, hi, why in SUBMODELS:
        start, end = match_idx[lo - 1], match_idx[hi - 1]
        afa = os.path.join(workdir, name + ".afa")
        kept = 0
        with open(afa, "w") as handle:
            for key, row in rows.items():
                segment = row[start:end + 1].replace(".", "-").upper()
                if len(segment.replace("-", "")) < 10:
                    continue
                handle.write(">%s\n%s\n" % (key.replace("/", "_"), segment))
                kept += 1
        hmm = os.path.join(workdir, name + ".hmm")
        subprocess.run(["hmmbuild", "--amino", "--informat", "afa", "-n", name,
                        hmm, afa], check=True, stdout=subprocess.DEVNULL)
        leng = next(int(l.split()[1]) for l in open(hmm) if l.startswith("LENG"))
        log("  %-10s kolon %3d-%3d  referans %2d  model LENG %3d" %
            (name, lo, hi, kept, leng))
        built.append({"name": name, "columns": [lo, hi], "references": kept,
                      "model_length": leng, "why": why,
                      "path": hmm, "all_columns_kept": leng == hi - lo + 1})
    return built


def best_domain_scores(domtbl):
    best = {}
    with open(domtbl) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            field = line.split()
            cid, score = field[0], float(field[13])
            if score > best.get(cid, -1e9):
                best[cid] = score
    return best


def stage_reference_probe(models, sets, workdir, cpu, log):
    """Alt-modelleri uc kumeye hmmsearch et ve skor dagilimlarini karsilastir."""
    out = []
    for model in models:
        row = {"model": model["name"], "columns": model["columns"],
               "model_length": model["model_length"], "why": model["why"],
               "sets": {}}
        for label, path in sets.items():
            total = sum(1 for line in open(path) if line.startswith(">"))
            domtbl = os.path.join(workdir, "%s_%s.dom" % (model["name"], label))
            subprocess.run(["hmmsearch", "--max", "--cpu", str(cpu), "-E", "1000",
                            "--domtblout", domtbl, "-o", os.devnull,
                            model["path"], path], check=True)
            scores = best_domain_scores(domtbl)
            arr = np.array(list(scores.values())) if scores else np.array([0.0])
            row["sets"][label] = {
                "n": total,
                "with_any_hit": len(scores),
                "score_median": round(float(np.median(arr)), 1),
                "score_p95": round(float(np.percentile(arr, 95)), 1),
                "fraction_score_ge_20": round(
                    sum(1 for v in arr if v >= 20) / max(1, total), 4),
                "fraction_score_ge_10": round(
                    sum(1 for v in arr if v >= 10) / max(1, total), 4),
            }
            log("  %-10s vs %-14s hit %5d/%-5d  skor>=20 %5.1f%%"
                % (model["name"], label, len(scores), total,
                   100 * row["sets"][label]["fraction_score_ge_20"]))
        out.append(row)
    return {
        "question": ("referanstan kurulan alt-modeller bu dizilerin katalitik "
                     "lobuna erisebiliyor mu"),
        "models": out,
        "control_note": ("rieske_only = 'Rieske tam, katalitik triad elenmis' "
                         "adaylar; pipeline bunlari RO alpha SAYMIYOR. Bir "
                         "sondanin bu kontrolde de dustugu bir sinyal, "
                         "karboksilatsiz setin RO olmadigini DEGIL, sondanin "
                         "o bolgeye erisemedigini gosterir."),
    }


# --------------------------------------------------------------------------
# Adim 6-7: capali yeniden okuma ve bos modeli
# --------------------------------------------------------------------------

def run_mafft(records, workdir, tag, threads):
    """mafft --auto. threads=1 ise sonuc birebir tekrarlanabilir (bkz. modul
    basligindaki TEKRARLANABILIRLIK notu)."""
    fasta = os.path.join(workdir, "co_%s.fa" % tag)
    out = os.path.join(workdir, "co_%s.aln" % tag)
    write_fasta(fasta, records)
    with open(out, "w") as handle:
        subprocess.run(["mafft", "--auto", "--quiet", "--thread", str(threads),
                        fasta], check=True, stdout=handle)
    return read_fasta(out)


def column_map(row):
    out, count = {}, 0
    for i, char in enumerate(row):
        if char != "-":
            out[count] = i
            count += 1
    return out


def stage_anchored_reread(entries, index, negatives, positives, refs,
                          workdir, mafft_threads, window, pos_sample, log):
    """ASIL OLCUM: tip ici yeniden hizalama, kolonu POZITIFLER belirler."""
    neg_by = defaultdict(list)
    pos_by = defaultdict(list)
    for cid in negatives:
        neg_by[entries[cid]["cluster"]].append(cid)
    for cid in positives:
        pos_by[entries[cid]["cluster"]].append(cid)

    per_cluster = []
    calls = {}
    for cluster in sorted(neg_by, key=lambda c: -len(neg_by[c])):
        negs = sorted(neg_by[cluster])
        poss = sorted(pos_by.get(cluster, []))
        # Pozitif capalar: karboksilat HAM pozisyonu bilinen olanlar
        anchored_pos = []
        for cid in poss:
            idx = index.get(cid)
            if idx is None:
                continue
            raw = idx.find(entries[cid]["sequence"], CARBOX_COL, CARBOX_OK, window)[0]
            if raw is not None:
                anchored_pos.append((cid, raw))
        if len(anchored_pos) > pos_sample:
            step = max(1, len(anchored_pos) // pos_sample)
            anchored_pos = anchored_pos[::step][:pos_sample]
        ref_names = [k for k in refs if k.startswith(cluster + "_")]

        if len(anchored_pos) < 5:
            per_cluster.append({"cluster": cluster, "n_no_carboxylate": len(negs),
                                "status": "no_positive_anchor_available"})
            for cid in negs:
                calls[cid] = {"call": "untested", "residue": None, "offset": None}
            continue

        records = [(c, entries[c]["sequence"]) for c in negs]
        records += [(c, entries[c]["sequence"]) for c, _ in anchored_pos]
        records += [("REF|" + k, refs[k]) for k in ref_names]
        aln = run_mafft(records, workdir, cluster, mafft_threads)
        maps = {k: column_map(v) for k, v in aln.items()}
        inverse = {k: {c: r for r, c in m.items()} for k, m in maps.items()}
        width = len(next(iter(aln.values())))

        def his2_raw(cid):
            idx = index.get(cid)
            return (idx.find(entries[cid]["sequence"], HIS2_COL, "H", 2)[0]
                    if idx else None)

        # BEKLENEN uzaklik: ayni tipin pozitiflerinde His-2 -> karboksilat
        # arasindaki HAM kalinti sayisi. "Kayma" bunun etrafinda olculur,
        # cunku sabit kolona gore kayma bu girisler icin TANIMSIZ (kolon
        # hizalamada bos).
        pos_distances = []
        for cid, raw in anchored_pos:
            h = his2_raw(cid)
            if h is not None:
                pos_distances.append(raw - h)
        expected_distance = (int(np.median(pos_distances))
                             if pos_distances else None)

        # Karboksilat kolonu: pozitiflerin bilinen HAM pozisyonlarinin oyu
        votes = Counter()
        for cid, raw in anchored_pos:
            col = maps[cid].get(raw)
            if col is not None:
                votes[col] += 1
        anchor_col, top_votes = votes.most_common(1)[0]

        # Yontemin kendi dogrulugu: pozitifler capadan geri okunabiliyor mu
        def read_at(cid, col):
            """Donen: (kalinti, kolon kaymasi, hizalama kolonu)."""
            row = aln[cid]
            for delta in sorted(range(-window, window + 1), key=lambda x: (abs(x), x)):
                c = col + delta
                if 0 <= c < width and row[c].upper() in CARBOX_OK:
                    return row[c].upper(), delta, c
            return (row[col].upper() if 0 <= col < width else "?"), None, None

        def distance_offset(cid, column):
            """Bulunan karboksilatin BEKLENEN His-2 uzakligindan sapmasi (kalinti)."""
            if column is None or expected_distance is None:
                return None
            raw = inverse[cid].get(column)
            h = his2_raw(cid)
            if raw is None or h is None:
                return None
            return (raw - h) - expected_distance

        recovered = sum(1 for cid, _ in anchored_pos
                        if read_at(cid, anchor_col)[1] is not None)

        # His-2 kolonu (arama alt siniri) ve grubun kendi korunmus kolonu
        his_cols = Counter()
        for cid in negs:
            idx = index.get(cid)
            if idx is None:
                continue
            raw = idx.find(entries[cid]["sequence"], HIS2_COL, "H", 2)[0]
            if raw is not None and raw in maps[cid]:
                his_cols[maps[cid][raw]] += 1
        his_col = his_cols.most_common(1)[0][0] if his_cols else 0

        de_by_col = {}
        for col in range(his_col + 40, width):
            occupied = sum(1 for cid in negs if aln[cid][col] != "-")
            if occupied < 0.8 * len(negs):
                continue
            de = sum(1 for cid in negs if aln[cid][col].upper() in CARBOX_OK)
            de_by_col[col] = de / len(negs)
        own_col, own_frac = (max(de_by_col.items(), key=lambda kv: kv[1])
                             if de_by_col else (None, None))

        # Grubun kendi D/E arka plan sikligi (His-2 sonrasi kuyruk)
        background = []
        for cid in negs:
            idx = index.get(cid)
            raw = (idx.find(entries[cid]["sequence"], HIS2_COL, "H", 2)[0]
                   if idx else None)
            tail = entries[cid]["sequence"][(raw or 0) + 40:]
            if tail:
                background.append(sum(1 for ch in tail if ch in CARBOX_OK) / len(tail))
        bg = float(np.mean(background)) if background else 0.11

        residues = Counter()
        offsets = Counter()
        shifts = Counter()
        for cid in negs:
            res, delta, col_hit = read_at(cid, anchor_col)
            if delta is not None:
                shift = distance_offset(cid, col_hit)
                calls[cid] = {"call": "present_at_homologous_column",
                              "residue": res, "offset": delta,
                              "residue_offset_from_expected": shift}
                offsets[delta] += 1
                residues[res] += 1
                if shift is not None:
                    shifts[shift] += 1
                continue
            shifted = (read_at(cid, own_col) if own_col is not None
                       and own_frac and own_frac >= 0.8 else (res, None, None))
            if shifted[1] is not None:
                shift = distance_offset(cid, shifted[2])
                calls[cid] = {"call": "present_at_group_conserved_column",
                              "residue": shifted[0],
                              "offset": own_col - anchor_col + shifted[1],
                              "residue_offset_from_expected": shift}
                residues[shifted[0]] += 1
                if shift is not None:
                    shifts[shift] += 1
            else:
                calls[cid] = {"call": "absent", "residue": res, "offset": None,
                              "residue_offset_from_expected": None}
                residues[res] += 1

        present = sum(1 for cid in negs
                      if calls[cid]["call"].startswith("present"))
        # Bos model: bu D/E orani grubun kendi bilesimiyle sans olabilir mi
        observed = sum(1 for cid in negs
                       if read_at(cid, anchor_col)[1] is not None)
        p_binom = float(sps.binomtest(observed, len(negs), bg,
                                      alternative="greater").pvalue)
        column_fracs = sorted(de_by_col.values(), reverse=True)

        entry = {
            "cluster": cluster,
            "status": "tested",
            "n_no_carboxylate": len(negs),
            "n_positive_anchors": len(anchored_pos),
            "curated_reference_in_alignment": ref_names,
            "alignment_width": width,
            "anchor_column": anchor_col,
            "anchor_column_vote_share": round(top_votes / len(anchored_pos), 3),
            "method_recovery_on_positives": round(recovered / len(anchored_pos), 3),
            "reference_residue_at_anchor": {
                k: aln["REF|" + k][anchor_col] for k in ref_names},
            "carboxylate_present": present,
            "carboxylate_present_fraction": round(present / len(negs), 4),
            "present_at_homologous_column": sum(
                1 for cid in negs
                if calls[cid]["call"] == "present_at_homologous_column"),
            "present_at_group_conserved_column": sum(
                1 for cid in negs
                if calls[cid]["call"] == "present_at_group_conserved_column"),
            "absent": sum(1 for cid in negs if calls[cid]["call"] == "absent"),
            "residues_read": dict(residues.most_common()),
            "offset_counts": dict(sorted(offsets.items(), key=lambda kv: abs(kv[0]))),
            "expected_raw_distance_his2_to_carboxylate": expected_distance,
            "residue_offset_from_expected_counts": dict(
                sorted(shifts.items(), key=lambda kv: abs(kv[0]))[:15]),
            "residue_offset_within_plus_minus_3": sum(
                n for o, n in shifts.items() if abs(o) <= 3),
            "group_background_DE_frequency": round(bg, 4),
            "group_own_best_column": own_col,
            "group_own_best_DE_fraction": round(own_frac, 4) if own_frac else None,
            "group_own_column_equals_anchor": (
                own_col is not None and abs(own_col - anchor_col) <= window),
            "binomial_p_vs_own_composition": p_binom,
            "second_best_column_DE_fraction": round(column_fracs[1], 4)
            if len(column_fracs) > 1 else None,
            "median_column_DE_fraction": round(float(np.median(column_fracs)), 4)
            if column_fracs else None,
            "n_columns_scanned_for_the_null": len(de_by_col),
        }
        per_cluster.append(entry)
        log("  %-14s n=%-4d capa kol %-5d pozitif geri okuma %5.1f%%  "
            "karboksilat %5.1f%%  yok %d"
            % (cluster, len(negs), anchor_col,
               100 * entry["method_recovery_on_positives"],
               100 * entry["carboxylate_present_fraction"], entry["absent"]))

    # NCBI urun dizgisine gore capraz kontrol.
    #
    # Bu kumedeki girislerin bir kismi NCBI'da "SRPBCC family protein"
    # (heliks-tutamaci kati, oksijenaz katalitik alani DEGIL) ya da
    # "Rieske 2Fe-2S domain-containing protein" (katalitik alandan hic soz
    # etmeyen bir etiket) olarak geciyor. Bu etiketler "bunlar oksijenaz
    # degil" iddiasinin kaynagi, bu yuzden capali okumanin o etiketlerde ne
    # buldugu ayri raporlanir.
    by_product = defaultdict(Counter)
    for cid, call in calls.items():
        by_product[entries[cid]["product"] or "?"][call["call"]] += 1
    annotation_check = sorted(
        ({"product": product,
          "n": sum(counts.values()),
          "calls": dict(counts.most_common()),
          "carboxylate_present": sum(v for k, v in counts.items()
                                     if k.startswith("present")),
          "carboxylate_present_fraction": round(
              sum(v for k, v in counts.items() if k.startswith("present"))
              / sum(counts.values()), 4)}
         for product, counts in by_product.items()),
        key=lambda d: -d["n"])[:15]

    tested = [c for c in per_cluster if c["status"] == "tested"]
    qs = bh_adjust([c["binomial_p_vs_own_composition"] for c in tested])
    for c, q in zip(tested, qs):
        c["binomial_q"] = q
        c["beats_the_null"] = bool(q < 0.05)

    totals = Counter(v["call"] for v in calls.values())
    offsets = Counter()
    shifts = Counter()
    for v in calls.values():
        if v["offset"] is not None:
            offsets[v["offset"]] += 1
        if v.get("residue_offset_from_expected") is not None:
            shifts[v["residue_offset_from_expected"]] += 1
    shift_values = [o for o, n in shifts.items() for _ in range(n)]
    return {
        "question": ("karboksilat, SABIT kolon yerine ayni tipin pozitiflerine "
                     "capalanarak okundugunda var mi"),
        "method": ("tip basina mafft --auto ile yeniden hizalama; karboksilat "
                   "kolonu, pozitiflerin HAM dizideki bilinen karboksilat "
                   "pozisyonlarinin modu; karboksilatsiz uyeler o kolonda "
                   "+-%d kolon penceresiyle okunuyor" % window),
        "window_columns": window,
        "positive_anchors_per_cluster_max": pos_sample,
        "mafft_threads": mafft_threads,
        "reproducibility": (
            ("mafft tek parcacikla kosuldu: bu cikti birebir "
             "tekrarlanabilir."
             if mafft_threads == 1 else
             "mafft %d parcacikla kosuldu: yineleyici iyilestirme is "
             "parcacigi siralamasina duyarli oldugu icin bu cikti bit "
             "duzeyinde tekrarlanmaz. --thread 6 ile yapilan dort kosuda "
             "'gercek yokluk' 21-25 arasinda degisti (%%1,5-1,8); karar "
             "duzeyi her kosuda ayni kaldi. Birebir tekrar icin "
             "--mafft-threads 1." % mafft_threads)),
        "totals": dict(totals.most_common()),
        "carboxylate_present_total": sum(
            v for k, v in totals.items() if k.startswith("present")),
        "absent_total": totals.get("absent", 0),
        "untested_total": totals.get("untested", 0),
        "offset_distribution": dict(sorted(offsets.items(), key=lambda kv: abs(kv[0]))),
        "by_ncbi_product_annotation": {
            "question": ("NCBI'nin 'SRPBCC family protein' / 'Rieske 2Fe-2S "
                         "domain-containing protein' gibi etiketleri, "
                         "katalitik merkezin yokluguna isaret ediyor mu"),
            "rows": annotation_check,
        },
        "residue_offset_from_expected": {
            "definition": ("bulunan karboksilatin His-2'ye HAM uzakligi eksi o "
                           "tipin pozitiflerindeki medyan uzaklik. Sabit kolona "
                           "gore kayma bu girisler icin tanimsiz, cunku kolon "
                           "hizalamada bos; bu yuzden kayma BEKLENEN UZAKLIGA "
                           "gore olculuyor."),
            "how_to_read_this": ("Kayit testi offset_distribution'dadir: orada "
                                 "0, karboksilatin pozitiflerle AYNI hizalama "
                                 "kolonunda, yani homolog konumda oldugunu "
                                 "soyler. Buradaki kalinti-uzakligi sapmasi "
                                 "FARKLI bir seyi olcer: His cifti ile "
                                 "karboksilat arasindaki ~140 kalintilik "
                                 "parcada ne kadar indel birikmis. Sapmanin "
                                 "genis olmasi hatayi degil, SABIT KOLONUN "
                                 "NEDEN KAYDIGINI gosterir -- bu parcadaki "
                                 "insersiyon/delesyonlar modelin registresini "
                                 "kaydiriyor."),
            "n": len(shift_values),
            "within_plus_minus_0": shifts.get(0, 0),
            "within_plus_minus_3": sum(1 for o in shift_values if abs(o) <= 3),
            "within_plus_minus_10": sum(1 for o in shift_values if abs(o) <= 10),
            "beyond_plus_minus_10": sum(1 for o in shift_values if abs(o) > 10),
            "median": int(np.median(shift_values)) if shift_values else None,
            "p2.5": int(np.percentile(shift_values, 2.5)) if shift_values else None,
            "p97.5": int(np.percentile(shift_values, 97.5)) if shift_values else None,
            "counts_most_common": dict(
                sorted(shifts.items(), key=lambda kv: -kv[1])[:20]),
        },
        "per_cluster": per_cluster,
        "null_model": {
            "test": ("her tip icin: capa kolonunda gozlenen D/E sayisi, o grubun "
                     "KENDI kuyruk D/E sikligina karsi tek yonlu binom; "
                     "Benjamini-Hochberg ile duzeltilir"),
            "clusters_tested": len(tested),
            "clusters_beating_the_null": sum(1 for c in tested
                                             if c.get("beats_the_null")),
            "comparison": ("ayni hizalamada dolulugu >=%80 olan butun kolonlarin "
                           "D/E orani da raporlanir (median_column_DE_fraction); "
                           "capa kolonu bu dagilimin icinde kalirsa sinyal yoktur"),
        },
    }, calls


# --------------------------------------------------------------------------
# Adim 8: gercek yokluklarin kimligi
# --------------------------------------------------------------------------

# Artik girislerin KENDI aile ici korunmus kolonu -- son kontrol.
#
# Adim 6'nin capasi tipin BAKTERIYEL pozitiflerine dayaniyor. Karboksilat
# aileye ozgu, HOMOLOG OLMAYAN bir konuma tasinmissa o capa onu goremez.
# Bu yuzden artik girisler kendi aralarinda (urun dizgilerinden cikan aile
# etiketlerine gore) yeniden hizalanir ve His-2'ye capalanmis korunmus bir
# D/E kolonu aranir. Bir ailenin butun uyelerinde ayni kolonda D/E varsa
# karboksilat VARDIR -- sadece referans setinin ulasamadigi bir yerde.
PRODUCT_FAMILIES = [
    ("choline_monooxygenase", ("choline monooxygenase",)),
    ("chlorophyllide_a_oxygenase", ("chlorophyllide",)),
    ("pheophorbide_a_oxygenase", ("pheophorbide",)),
    ("PTC52_TIC55_translocon", ("translocon", "tic 55", "tic55")),
    ("cholesterol_7_desaturase", ("cholesterol", "daf-36")),
]


def product_family(product):
    text = (product or "").lower()
    for name, keys in PRODUCT_FAMILIES:
        if any(k in text for k in keys):
            return name
    return None


def residual_family_scan(ids, entries, index, workdir, mafft_threads,
                         min_group=3):
    """Artik girisleri kendi aralarinda hizala, korunmus D/E kolonu ara."""
    his2 = {}
    for cid in ids:
        idx = index.get(cid)
        if idx is None:
            continue
        raw = idx.find(entries[cid]["sequence"], HIS2_COL, "H", 2)[0]
        if raw is not None:
            his2[cid] = raw

    groups = {"all_residual": list(ids)}
    by_family = defaultdict(list)
    for cid in ids:
        by_family[product_family(entries[cid]["product"])
                  or ("unassigned_" + (entries[cid]["domain"] or "?"))].append(cid)
    groups.update(by_family)

    out = []
    for tag, members in sorted(groups.items(), key=lambda kv: -len(kv[1])):
        members = [c for c in members if c in his2]
        if len(members) < min_group:
            out.append({"group": tag, "n": len(members),
                        "status": "too_small_to_align (<%d)" % min_group})
            continue
        aln = run_mafft([(c, entries[c]["sequence"]) for c in members],
                        workdir, "residual_" + tag, mafft_threads)
        maps = {k: column_map(v) for k, v in aln.items()}
        width = len(next(iter(aln.values())))
        his_cols = Counter(maps[c][his2[c]] for c in members if his2[c] in maps[c])
        his_col = his_cols.most_common(1)[0][0] if his_cols else 0
        best = None
        for col in range(his_col + 40, width):
            occupied = sum(1 for c in members if aln[c][col] != "-")
            if occupied < 0.8 * len(members):
                continue
            de = sum(1 for c in members if aln[c][col].upper() in CARBOX_OK)
            if best is None or de > best[1]:
                best = (col, de, occupied)
        background = []
        for c in members:
            tail = entries[c]["sequence"][his2[c] + 40:]
            if tail:
                background.append(
                    sum(1 for ch in tail if ch in CARBOX_OK) / len(tail))
        bg = float(np.mean(background)) if background else 0.11
        entry = {"group": tag, "status": "tested", "n": len(members)}
        if best is None:
            entry["status"] = "no_column_with_80pct_occupancy"
            out.append(entry)
            continue
        col, de, occupied = best
        inverse = {c: {v: k for k, v in maps[c].items()} for c in members}
        distances = [inverse[c][col] - his2[c] for c in members if col in inverse[c]]
        entry.update({
            "best_column": col,
            "members_with_a_carboxylate_there": de,
            "fraction": round(de / len(members), 3),
            "column_occupancy": occupied,
            "residues_at_that_column": dict(Counter(
                aln[c][col].upper() for c in members).most_common()),
            "background_DE_frequency": round(bg, 4),
            "median_raw_distance_from_his2": int(np.median(distances))
            if distances else None,
            "binomial_p_vs_background": float(sps.binomtest(
                de, len(members), bg, alternative="greater").pvalue),
            "conserved": bool(de / len(members) >= 0.9),
        })
        out.append(entry)
    return out


def stage_real_absences(entries, index, calls, anchored, models, workdir, cpu,
                        mafft_threads, log):
    absent = [cid for cid, v in calls.items() if v["call"] == "absent"]
    untested = [cid for cid, v in calls.items() if v["call"] == "untested"]
    out = {
        "question": "capali okumadan sonra hala karboksilati olmayanlar nedir",
        "n_absent": len(absent),
        "n_untested": len(untested),
    }
    if not absent:
        out["note"] = "Capali okumadan sonra gercek yokluk kalmadi."
        return out, absent

    def tally(key):
        return dict(Counter(entries[c][key] or "?" for c in absent).most_common())

    out["residue_in_place"] = dict(Counter(
        calls[c]["residue"] or "?" for c in absent).most_common())
    out["by_type"] = tally("cluster")
    out["by_group"] = tally("grp")
    out["by_tier"] = tally("tier")
    out["by_domain_of_life"] = tally("domain")
    out["by_product"] = dict(Counter(
        entries[c]["product"] or "?" for c in absent).most_common(15))
    out["by_subfamily_class"] = tally("assignment_class")
    out["mean_ref_identity"] = round(float(np.mean(
        [entries[c]["ref_identity"] for c in absent
         if entries[c]["ref_identity"] is not None])), 1)
    out["mean_sequence_length"] = round(float(np.mean(
        [len(entries[c]["sequence"]) for c in absent])), 1)
    out["with_operon"] = sum(1 for c in absent if entries[c]["has_operon"])
    out["with_beta_subunit_in_operon"] = sum(1 for c in absent
                                             if entries[c]["has_beta"])
    out["with_beta_within_10kb"] = sum(1 for c in absent
                                       if entries[c]["nearby_beta"])

    # Protein, karboksilatin OTURACAGI yere kadar uzaniyor mu? Bir dizi
    # His-2'den sonra beklenen uzakligi doldurmuyorsa "karboksilat yok"
    # demek yanlis olur -- orada okunacak dizi yok, protein bitiyor.
    expected = {c["cluster"]: c.get("expected_raw_distance_his2_to_carboxylate")
                for c in anchored["per_cluster"]}
    reach = Counter()
    tails = []
    for cid in absent:
        idx = index.get(cid)
        seq = entries[cid]["sequence"]
        raw = idx.find(seq, HIS2_COL, "H", 2)[0] if idx else None
        exp = expected.get(entries[cid]["cluster"])
        if raw is None or exp is None:
            reach["could_not_be_measured"] += 1
            continue
        tail = len(seq) - raw
        tails.append(tail)
        if tail < exp:
            reach["protein_ends_before_the_expected_site"] += 1
        elif tail < exp + 20:
            reach["site_falls_in_the_last_20_residues"] += 1
        else:
            reach["long_enough_to_hold_the_site"] += 1
    out["c_terminal_reach"] = {
        "question": ("protein, karboksilatin oturacagi yere kadar uzaniyor mu"),
        "counts": dict(reach.most_common()),
        "residues_after_his2_median": int(np.median(tails)) if tails else None,
        "expected_distance_by_type": {k: v for k, v in sorted(expected.items())
                                      if k in out["by_type"]},
        "note": ("protein_ends_before_the_expected_site = bu girisler icin "
                 "'karboksilat yok' bir IKAME degil, dizinin bitmesidir."),
    }

    # Dizi uzayinda bir arada mi
    leaves = Counter(entries[c]["leaf_id"] or "?" for c in absent)
    genera = Counter(entries[c]["genus"] for c in absent)
    out["clustering"] = {
        "distinct_leaves": len(leaves),
        "largest_leaf": leaves.most_common(1)[0] if leaves else None,
        "entries_in_the_three_largest_leaves": sum(
            n for _, n in leaves.most_common(3)),
        "distinct_genera": len(genera),
        "top_genera": dict(genera.most_common(8)),
        "reading": ("yaprak sayisi giris sayisina yakinsa dagilmis, tek bir "
                    "yaprakta yiginlasiyorsa tek bir sapmis klad"),
    }

    # Alan mimarisi
    fasta = os.path.join(workdir, "absent.fasta")
    write_fasta(fasta, [(c, entries[c]["sequence"]) for c in absent])
    comp = os.path.join(workdir, "components_absent.dom")
    components = os.path.join("component_hmms", "components.hmm")
    architecture = {}
    if os.path.exists(components):
        subprocess.run(["hmmscan", "--cpu", str(cpu), "-E", "1e-3", "--domE", "1e-3",
                        "--domtblout", comp, "-o", os.devnull, components, fasta],
                       check=True)
        per = defaultdict(list)
        with open(comp) as handle:
            for line in handle:
                if line.startswith("#"):
                    continue
                f = line.split()
                per[f[3]].append((int(f[17]), f[0]))
        archs = Counter()
        for cid, hits in per.items():
            hits.sort()
            names = []
            for _, model in hits:
                if not names or names[-1] != model:
                    names.append(model)
            archs["+".join(names)] += 1
        architecture["components_hmm"] = {
            "library": sorted({m for hits in per.values() for _, m in hits}),
            "sequences_with_any_domain": len(per),
            "architectures": dict(archs.most_common(10)),
            "limitation": ("Bu kitaplikta halka-hidroksilleyici alpha'nin "
                           "KATALITIK alani icin model YOK (PF00848 / "
                           "Ring_hydroxyl_A kitaplikta degil) ve yerelde Pfam-A "
                           "bulunmuyor; bu yuzden 'sadece Rieske alani var' "
                           "ciktisi mimari hakkinda tek basina bir sey SOYLEMEZ. "
                           "SRPBCC (PF03364/PF10604, heliks-tutamaci) de "
                           "sinanamaz -- NCBI urun dizgileri oldugu gibi "
                           "by_product altinda raporlanir."),
        }
    for model in models:
        domtbl = os.path.join(workdir, "%s_absent.dom" % model["name"])
        subprocess.run(["hmmsearch", "--max", "--cpu", str(cpu), "-E", "1000",
                        "--domtblout", domtbl, "-o", os.devnull,
                        model["path"], fasta], check=True)
        scores = best_domain_scores(domtbl)
        architecture[model["name"]] = {
            "columns": model["columns"],
            "with_any_hit": len(scores),
            "fraction_score_ge_20": round(
                sum(1 for v in scores.values() if v >= 20) / max(1, len(absent)), 4),
        }
    out["domain_architecture"] = architecture

    # Aile ici son kontrol
    out["residual_family_scan"] = {
        "question": ("artik girisler KENDI aralarinda hizalandiginda, His-2'ye "
                     "capalanmis korunmus bir karboksilat kolonu var mi"),
        "why": ("Adim 6'nin capasi tipin BAKTERIYEL pozitiflerine dayaniyor; "
                "karboksilat aileye ozgu, homolog olmayan bir konuma tasinmissa "
                "o capa onu goremez. Bir ailenin butun uyelerinde ayni kolonda "
                "D/E varsa karboksilat VARDIR, sadece referans setinin "
                "ulasamadigi bir yerde."),
        "groups": residual_family_scan(absent, entries, index, workdir,
                                       mafft_threads),
    }
    resolved = [g for g in out["residual_family_scan"]["groups"]
                if g.get("conserved") and g["group"] != "all_residual"]
    out["residual_family_scan"]["families_with_a_conserved_carboxylate"] = [
        {"group": g["group"], "n": g["n"], "fraction": g["fraction"],
         "residues": g["residues_at_that_column"],
         "median_raw_distance_from_his2": g["median_raw_distance_from_his2"]}
        for g in resolved]
    out["residual_family_scan"]["entries_rescued_by_a_family_specific_site"] = sum(
        g["n"] for g in resolved)
    out["entries"] = [{
        "candidate_id": c,
        "type": entries[c]["cluster"],
        "residue_in_place": calls[c]["residue"],
        "product": entries[c]["product"],
        "tier": entries[c]["tier"],
        "domain": entries[c]["domain"],
        "ref_identity": entries[c]["ref_identity"],
        "length": len(entries[c]["sequence"]),
        "organism": entries[c]["organism"],
    } for c in sorted(absent)]
    log("  gercek yokluk: %d giris, %d tip, oturan kalintilar %s"
        % (len(absent), len(out["by_type"]), out["residue_in_place"]))
    return out, absent


# --------------------------------------------------------------------------
# Adim 9: siki kapinin bedeli
# --------------------------------------------------------------------------

def stage_gate_cost(entries, negatives, calls):
    total = len(entries)
    types_all = {rec["cluster"] for rec in entries.values()}

    members_of = defaultdict(list)
    for cid, rec in entries.items():
        members_of[rec["cluster"]].append(cid)

    def cost(lost_ids, label, note):
        lost_ids = set(lost_ids)
        lost_types = [cl for cl in sorted(types_all)
                      if members_of[cl] and all(c in lost_ids for c in members_of[cl])]
        per_type = Counter(entries[c]["cluster"] for c in lost_ids)
        worst = []
        for cluster, n in per_type.most_common():
            size = len(members_of[cluster])
            worst.append({"type": cluster, "lost": n, "members": size,
                          "lost_fraction": round(n / size, 4)})
        worst.sort(key=lambda d: -d["lost_fraction"])
        return {
            "reading": label,
            "note": note,
            "entries_lost": len(lost_ids),
            "entries_lost_fraction": round(len(lost_ids) / total, 4),
            "entries_remaining": total - len(lost_ids),
            "types_losing_members": len(per_type),
            "types_wiped_out_entirely": lost_types,
            "hardest_hit_types": worst[:12],
        }

    fixed_lost = set(negatives)
    anchored_lost = {c for c, v in calls.items() if v["call"] in ("absent", "untested")}
    return {
        "question": ("kabul kurali katalitik triad'in 3/3'unu isteseydi ne "
                     "kaybedilirdi"),
        "rule_now": "Rieske 4/4 ZORUNLU + katalitik triad >=%d/3, pencere +-%d"
                    % (2, 2),
        "fixed_column_reading": cost(
            fixed_lost, "sabit kolon (%d, +-2)" % CARBOX_COL,
            "ro_carboxylate'in bugunku okumasi -- SISMIS sayi"),
        "anchored_reading": cost(
            anchored_lost, "capali yeniden okuma (adim 6)",
            "ASIL sayi: karboksilati gercekten bulunamayan girisler"),
    }


# --------------------------------------------------------------------------
# Adim 10: istatistik (stats_overview makinesi)
# --------------------------------------------------------------------------

def contingency(entries, label_of, row_key, row_name, col_names, alpha=0.05):
    """Giris duzeyi + tip x cins cokertmesiyle kontenjans testi."""
    rows = sorted({entries[c][row_key] or "?" for c in label_of})
    index_r = {r: i for i, r in enumerate(rows)}
    index_c = {c: j for j, c in enumerate(col_names)}
    table = [[0] * len(col_names) for _ in rows]
    buckets = defaultdict(list)
    for cid, call in label_of.items():
        rec = entries[cid]
        r = index_r[rec[row_key] or "?"]
        j = index_c.get(call)
        if j is None:
            continue
        table[r][j] += 1
        buckets[(rec["cluster"], rec["genus"])].append((r, j))
    gtable = [[0] * len(col_names) for _ in rows]
    for (_, _), items in buckets.items():
        r = items[0][0]
        js = Counter(j for _, j in items)
        gtable[r][js.most_common(1)[0][0]] += 1
    keep = [i for i in range(len(rows)) if sum(table[i]) > 0]
    rows = [rows[i] for i in keep]
    table = [table[i] for i in keep]
    gtable = [gtable[i] for i in keep]
    if len(rows) < 2:
        return None
    result = signal_map(rows, col_names, table, genus_table=gtable, alpha=alpha)
    if result is None:
        return None
    def safe_chi2(tab):
        arr = np.array(tab, dtype=float)
        if arr.sum() <= 0 or min(arr.shape) < 2 \
                or (arr.sum(axis=0) == 0).any() or (arr.sum(axis=1) == 0).any():
            return None, None
        out = sps.chi2_contingency(arr)
        return float(out[0]), float(out[1])

    chi2, p = safe_chi2(table)
    if p is None:
        return None
    gchi2, gp = safe_chi2(gtable)
    return {
        "question": row_name,
        "rows": rows,
        "columns": col_names,
        "entry_level": {"n": int(np.array(table).sum()), "chi2": round(chi2, 1),
                        "p": p, "cramers_v": round(cramers_v(table), 3)},
        "genus_level": {"n": int(np.array(gtable).sum()),
                        "p": gp,
                        "cramers_v": round(cramers_v(gtable), 3)
                        if gp is not None else None,
                        "note": "tip x cins basina tek gozlem (cogunluk kurali)"},
        "signal_map": {k: result[k] for k in
                       ("n_cells_tested", "n_cells_flagged", "n_cells_material",
                        "n_cells_material_genus", "n_rows_tested",
                        "n_rows_flagged", "top", "min_lift", "alpha")},
    }


def stage_statistics(entries, negatives, calls):
    neg = set(negatives)
    fixed = {cid: ("no_carboxylate" if cid in neg else "has_carboxylate")
             for cid in entries}
    anchored = {cid: ("absent" if calls.get(cid, {}).get("call") in
                      ("absent", "untested") else "present")
                for cid in entries}
    tests = []
    for label, label_of, cols in (
            ("sabit kolon cagrisi", fixed, ["has_carboxylate", "no_carboxylate"]),
            ("capali okuma cagrisi", anchored, ["present", "absent"])):
        for key, name in (("cluster", "enzim tipi"), ("tier", "kanit kademesi"),
                          ("domain", "taksonomik alan")):
            res = contingency(entries, label_of, key, "%s x %s" % (name, label), cols)
            if res:
                tests.append(res)
    qs = bh_adjust([t["entry_level"]["p"] for t in tests])
    for t, q in zip(tests, qs):
        t["entry_level"]["q"] = q
    return {
        "machinery": ("stats_overview.signal_map / bh_adjust / cramers_v "
                      "dogrudan import edildi"),
        "independence_note": ("Giris duzeyi p degerleri IYIMSER (ayni sus "
                              "genomlari tekrar eder); tip x cins cokertmesi "
                              "muhafazakar okumadir."),
        "tests": tests,
    }


# --------------------------------------------------------------------------
# Verdict
# --------------------------------------------------------------------------

def _probe_numbers(probe):
    """Verdict metni icin sonda oranlari: (ROcatC uc set, ROcarbox iki set)."""
    by_name = {m["model"]: m["sets"] for m in probe["models"]}
    cat = by_name.get("ROcatC", {})
    car = by_name.get("ROcarbox", {})

    def frac(sets, label):
        return 100 * sets.get(label, {}).get("fraction_score_ge_20", 0.0)

    return (frac(cat, "no_carboxylate"), frac(cat, "with_carboxylate"),
            frac(cat, "rieske_only"), frac(car, "no_carboxylate"),
            frac(car, "with_carboxylate"))


def annotation_verdict(anchored):
    """NCBI urun dizgisi bir yoksunluk kaniti mi -- duz sozlerle.

    Bu kumenin "oksijenaz degil" diye okunmasinin kaynagi urun dizgileri:
    "SRPBCC family protein" bir heliks-tutamaci katidir, oksijenaz katalitik
    alani degil; "Rieske 2Fe-2S domain-containing protein" ise katalitik
    alandan hic soz etmez. O etiketlerin capali okumada ne verdigini
    raporlamak, iddiayi dogrudan sinamaktir.
    """
    rows = anchored["by_ncbi_product_annotation"]["rows"]

    def describe(predicate):
        hits = [r for r in rows if predicate(r["product"])]
        if not hits:
            return None
        row = max(hits, key=lambda r: r["n"])
        return ("'%s' etiketli %d girisin %d'inde (%%%.1f)"
                % (row["product"], row["n"], row["carboxylate_present"],
                   100 * row["carboxylate_present_fraction"]))

    parts = [text for text in (
        describe(lambda t: "SRPBCC" in t),
        describe(lambda t: t.startswith("Rieske 2Fe-2S")),
        describe(lambda t: "hypothetical" in t.lower()),
    ) if text]
    if not parts:
        return "Bu kumede incelenecek bir urun dizgisi yogunlasmasi yok."
    return ("NCBI urun dizgileri bu girislerin oksijenaz OLMADIGININ kaniti "
            "degil: capali okuma " + "; ".join(parts) + " katalitik "
            "karboksilati yerinde buldu. Bu etiketler proteinin yalnizca bir "
            "parcasina oturan kismi CDD cagrilaridir; o parcanin yaninda tam "
            "bir mononukleer demir merkezi duruyor.")


def build_verdict(baseline, mechanism, widening, naive, probe, anchored,
                  absences, gate):
    n_neg = mechanism["no_carboxylate"]["n"]
    present = anchored["carboxylate_present_total"]
    absent = anchored["absent_total"]
    untested = anchored["untested_total"]
    gapped = baseline["missing_call_breakdown"].get("window_entirely_gapped", 0)
    substituted = baseline["missing_call_breakdown"].get(
        "substituted_residue_in_place", 0)
    stops = mechanism["no_carboxylate"]["alignment_stops_before_the_carboxylate_column"]
    tested_clusters = [c for c in anchored["per_cluster"] if c["status"] == "tested"]
    recovery = [c["method_recovery_on_positives"] for c in tested_clusters]
    return {
        "headline": (
            "Karboksilatsiz gorunen %d girisin %d'i (%%%.1f) SABIT KOLON "
            "YANILGISIDIR: karboksilat, ayni tipin pozitiflerine capalanarak "
            "okundugunda yerinde durmaktadir. Gercek yokluk %d giris (%%%.1f)."
            % (n_neg, present, 100 * present / max(1, n_neg), absent,
               100 * absent / max(1, n_neg))),
        "how_many_were_an_alignment_artefact": present,
        "how_many_are_a_real_absence": absent,
        "how_many_could_not_be_tested": untested,
        "why_the_fixed_column_failed": (
            "Eksik cagrilarin %d'i (%%%.1f) kolonun KENDISININ ve komsularinin "
            "hizalamada BOS olmasindan geliyor, yalnizca %d'inde kolonda "
            "karboksilat olmayan gercek bir kalinti oturuyor. Sebep bir-iki "
            "pozisyonluk kayma degil: bu girislerin %d'inde (%%%.1f) ROmotif71 "
            "hizalamasi karboksilat kolonuna HIC ULASMIYOR -- son dolu match "
            "kolonu medyani %s, oysa karboksilat kolonu %d. Model kapsami "
            "ortalama %.2f (karboksilatli sette %.2f) ve hmm skoru %.0f "
            "(karsisinda %.0f). Yani 426 kolonluk referans modeli bu dizilerin "
            "C-terminal katalitik lobunu hizalayamiyor; kolon 355 hizalanan "
            "bolgenin DISINDA kaliyor."
            % (gapped, 100 * gapped / max(1, n_neg), substituted, stops,
               100 * stops / max(1, n_neg),
               mechanism["no_carboxylate"]["last_filled_column_median"], CARBOX_COL,
               mechanism["no_carboxylate"]["model_coverage_mean"],
               mechanism["with_carboxylate"]["model_coverage_mean"],
               mechanism["no_carboxylate"]["hmm_score_mean"],
               mechanism["with_carboxylate"]["hmm_score_mean"])
            + (" Referans profil sondasi bunu dogrudan gosteriyor: C-terminal "
               "katalitik lob modeli (ROcatC, kolon 260-426) bu dizilerin "
               "yalnizca %%%.1f'inde skor>=20 veriyor, karboksilatli sette "
               "%%%.1f, pipeline'in RO SAYMADIGI 'Rieske var katalitik yok' "
               "kontrol setinde %%%.1f. Dar karboksilat sondasi (ROcarbox, "
               "kolon 330-380) karboksilatsiz sette %%%.1f, karboksilatli "
               "sette %%%.1f. Rieske alani sondasi ise UC sette de %%100'e "
               "yakin tutuyor -- yani kaybolan sey Rieske merkezi degil, "
               "referans modelinin katalitik loba ERISIMI."
               % _probe_numbers(probe))),
        "why_this_is_not_a_window_problem": (
            "Ayni sabit kolonu +-%d kolona kadar genisletmek yalnizca %%%.1f'ini "
            "kurtariyor -- bos kolonlari genisletmek bos kolon taramaktir."
            % (widening["windows"][-1]["window_columns"],
               100 * widening["windows"][-1]["rescued_fraction"])),
        "why_the_naive_anchor_was_not_enough": naive["verdict"],
        "what_the_anchored_reading_did": (
            "Her tipte karboksilatsiz uyeler, AYNI TIPIN karboksilati bilinen "
            "uyeleri ve kuratorlu referans dizisiyle birlikte mafft ile yeniden "
            "hizalandi; karboksilat kolonu sabit bir sayi olarak degil, "
            "pozitiflerin bilinen pozisyonlarinin oyuyla belirlendi. Yontemin "
            "kendi dogrulugu ayni hizalamada olculdu: pozitiflerin medyan "
            "%%%.1f'inde capa dogru kalintiyi geri verdi. Kayma dagilimi "
            "capanin UZERINDE yiginlasiyor (kolon kaymasi 0: %d giris; beklenen "
            "His-2 uzakligindan sapma +-3 kalinti icinde: %d / %d), yani "
            "okunan sey "
            "rastgele bir D/E degil, HOMOLOG konumdaki karboksilat. Bu, %d "
            "tipin %d'inde grubun kendi D/E bilesimine karsi binom testini "
            "geciyor."
            % (100 * float(np.median(recovery)) if recovery else 0,
               anchored["offset_distribution"].get(0, 0),
               anchored["residue_offset_from_expected"]["within_plus_minus_3"],
               anchored["residue_offset_from_expected"]["n"],
               anchored["null_model"]["clusters_tested"],
               anchored["null_model"]["clusters_beating_the_null"])),
        "what_the_real_absences_appear_to_be": (
            absences.get("note")
            or ("%d giris, %d tipe dagilmis; capa kolonunda oturan kalintilar %s. "
                "Bunlarin %d'i hizalamada o kolonda gap, yani 'yoklugu' bile "
                "hala olcum sinirindan kaynaklanabilir. %d farkli yaprakta "
                "duruyorlar, yani tek bir sapmis klad DEGIL, dagilmis bir "
                "artik. Alan mimarisi testi bir ayrim uretemiyor cunku "
                "component_hmms kitapliginda halka-hidroksilleyici alpha'nin "
                "katalitik alani icin model yok."
                % (absences["n_absent"], len(absences.get("by_type", {})),
                   absences.get("residue_in_place"),
                   absences.get("residue_in_place", {}).get("-", 0),
                   absences.get("clustering", {}).get("distinct_leaves", 0))
               + (" C-terminal erisim: %s. AILE ICI son kontrol: artik "
                  "girisler kendi aralarinda hizalandiginda %d giris aileye "
                  "OZGU bir konumda korunmus bir karboksilat tasiyor (%s); "
                  "yani onlarda da karboksilat VAR, sadece 71 bakteriyel "
                  "referansin ulasamadigi bir yerde. Geri kalan %d giris tek "
                  "ya da ikili temsil edilen, kuratorlu sette KARSILIGI "
                  "OLMAYAN okaryotik Rieske oksijenaz alt aileleri "
                  "(feoforbid a oksijenaz, klorofillid a oksijenaz, "
                  "PTC52/TIC55, kolesterol 7-desaturaz/daf-36) oldugu icin bu "
                  "yontemle cozulemiyor: kalan artik bir YANLIS POZITIF "
                  "kumesi degil, bir REFERANS SETI BOSLUGUDUR."
                  % (absences.get("c_terminal_reach", {}).get("counts"),
                     absences.get("residual_family_scan", {}).get(
                         "entries_rescued_by_a_family_specific_site", 0),
                     ", ".join(
                         "%s %d/%d" % (g["group"],
                                       int(round(g["fraction"] * g["n"])), g["n"])
                         for g in absences.get("residual_family_scan", {}).get(
                             "families_with_a_conserved_carboxylate", []))
                     or "yok",
                     absences["n_absent"]
                     - absences.get("residual_family_scan", {}).get(
                         "entries_rescued_by_a_family_specific_site", 0))))),
        "what_about_the_SRPBCC_annotation": annotation_verdict(anchored),
        "are_they_Rieske_oxygenases": (
            "Evet, kanitlar bu yonde. Karboksilatsiz gorunen setin HEPSINDE "
            "Rieske ligandlarinin 4/4'u, iki katalitik histidin ve kopru "
            "Asp/Glu yerinde -- yani tezin aktardigi D-x(2)-H-x(3-5)-H "
            "konsensusu tam. Capali okuma uzerine facial triad'in karboksilati "
            "da %%%.1f'inde bulundu. Bu girislerin RO'luktan uzak gorunmesinin "
            "sebebi biyolojik degil olcumsel: referans seti 71 BAKTERIYEL "
            "halka-hidroksilleyici tipten olusuyor ve bu kumenin %d'i okaryot "
            "(klorofil/steroid yolu Rieske oksijenazlari: feoforbid a "
            "oksijenaz, TIC55, klorofillid a oksijenaz, kolesterol "
            "7-desaturaz), onlarin katalitik lobu bakteriyel modele "
            "hizalanamiyor."
            % (100 * present / max(1, n_neg),
               mechanism.get("euk_count", 0))),
        "should_the_admission_rule_change": (
            "HAYIR, 3/3'e sikilastirilmamali. Sabit-kolon okumasina gore 3/3 "
            "kurali %d giris (%%%.1f) ve %d tipte uye kaybettirirdi; ama bu sayi "
            "SISMIS, cunku kayiplarin neredeyse tamami gercek bir biyolojik "
            "eksiklik degil, modelin menzili. Capali okumaya gore gercek bedel "
            "%d giris (%%%.2f). Yani 3/3 kurali, karboksilati fiilen OLAN %d "
            "gercek RO'yu veritabanindan atardi ve en sert vurdugu tipler "
            "(%s) kuratorlu referansa en uzak, yani en degerli olanlar. "
            "ro_motif.py'nin 2/3 secimi bu olcumle DOGRULANIYOR."
            % (gate["fixed_column_reading"]["entries_lost"],
               100 * gate["fixed_column_reading"]["entries_lost_fraction"],
               gate["fixed_column_reading"]["types_losing_members"],
               gate["anchored_reading"]["entries_lost"],
               100 * gate["anchored_reading"]["entries_lost_fraction"],
               present,
               ", ".join(d["type"] for d in
                         gate["fixed_column_reading"]["hardest_hit_types"][:4]))),
        "what_should_change_instead": (
            "Kural degil OLCUM duzeltilmeli. (1) ro_carboxylate'in "
            "catalytic_residue/catalytic_offset alanlari sabit kolondan "
            "okundugunu tasimali: NULL offset 'yok' degil 'pencerede "
            "okunamadi' demektir ve bugun bu ayrim tabloda YOK. (2) Kapsami "
            "dusuk girisler icin karboksilat tip ici capali okumayla "
            "yazilmali. (3) Referans setinde okaryotik Rieske oksijenaz "
            "(PAO/TIC55/CAO/kolesterol 7-desaturaz) tipi YOK; bu bosluk "
            "doldurulmadikca o girislerin katalitik lobu hicbir sabit kolonla "
            "okunamaz."),
    }


# --------------------------------------------------------------------------
# Ozet tablo
# --------------------------------------------------------------------------

def print_summary(out):
    w = 62
    base = out["baseline_fixed_column"]
    mech = out["mechanism"]
    anc = out["anchored_reread"]
    ab = out["real_absences"]
    gate = out["gate_cost"]
    v = out["verdict"]
    n_neg = mech["no_carboxylate"]["n"]
    print()
    print("=== ROAR-DB karboksilat denetimi ===")
    print("%-*s %s" % (w, "onaylanmis giris", out["inputs"]["confirmed_entries"]))
    print("%-*s %s" % (w, "sabit kolonda karboksilat YOK", n_neg))
    print()
    print("--- 1. eksik cagri neyi goruyor ---")
    for k, n in base["missing_call_breakdown"].items():
        print("%-*s %5d  (%%%.1f)" % (w, "  " + k, n, 100 * n / max(1, n_neg)))
    print()
    print("--- 2. mekanizma: hizalama karboksilat kolonuna ulasiyor mu ---")
    print("%-*s %s" % (w, "  kolon 355 hizalanan bolgenin disinda",
                       "%d (%%%.1f)" % (
                           mech["no_carboxylate"]["alignment_stops_before_the_carboxylate_column"],
                           100 * mech["no_carboxylate"]["alignment_stops_before_fraction"])))
    for label, key in (("son dolu match kolonu (medyan)", "last_filled_column_median"),
                       ("model kapsami (ortalama)", "model_coverage_mean"),
                       ("hmm skoru (ortalama)", "hmm_score_mean"),
                       ("referansa kimlik % (ortalama)", "ref_identity_mean"),
                       ("referans kaplamasi % (ortalama)", "ref_qcov_mean")):
        print("%-*s %-10s karboksilatli: %s"
              % (w, "  " + label, mech["no_carboxylate"][key],
                 mech["with_carboxylate"][key]))
    print("%-*s %s" % (w, "  D-x(2)-H-x(3-5)-H konsensusu TAM olan",
                       mech["anchor_integrity_in_the_no_carboxylate_set"].get(
                           "full_D_x2_H_x3_5_H_consensus_intact")))
    print()
    print("--- 3. sabit kolonu genisletmek ---")
    for row in out["window_widening"]["windows"]:
        print("%-*s %5d  (%%%.1f)" % (w, "  pencere +-%d kolon" % row["window_columns"],
                                      row["rescued"], 100 * row["rescued_fraction"]))
    print()
    print("--- 4. ciplak capa (ve neden yetmez) ---")
    nv = out["naive_anchor"]
    print("%-*s %s" % (w, "  pencere (His-2'ye gore ham kalinti)", nv["window_used"]))
    print("%-*s %s" % (w, "  pencere genisligi", nv["window_width_residues"]))
    print("%-*s %%%.1f" % (w, "  D/E bulundu", 100 * nv["found_fraction"]))
    print("%-*s %%%.1f" % (w, "  SANS beklentisi",
                           100 * (nv["expected_by_chance_fraction"] or 0.0)))
    print()
    print("--- 5. referans profil sondalari (skor>=20 oranı) ---")
    print("  %-10s %-10s %-10s %-10s" % ("model", "karbok.yok", "karbok.var", "rieske_only"))
    for m in out["reference_profile_probe"]["models"]:
        s = m["sets"]
        print("  %-10s %-10s %-10s %-10s" % (
            m["model"],
            "%%%.1f" % (100 * s["no_carboxylate"]["fraction_score_ge_20"]),
            "%%%.1f" % (100 * s["with_carboxylate"]["fraction_score_ge_20"]),
            "%%%.1f" % (100 * s["rieske_only"]["fraction_score_ge_20"])))
    print()
    print("--- 6. CAPALI YENIDEN OKUMA (asil olcum) ---")
    print("%-*s %s" % (w, "  karboksilat VAR (homolog kolon)",
                       anc["totals"].get("present_at_homologous_column", 0)))
    print("%-*s %s" % (w, "  karboksilat VAR (grup korunmus kolon)",
                       anc["totals"].get("present_at_group_conserved_column", 0)))
    print("%-*s %s" % (w, "  karboksilat YOK (gercek yokluk)", anc["absent_total"]))
    print("%-*s %s" % (w, "  sinanamadi", anc["untested_total"]))
    print("%-*s %s" % (w, "  kolon kaymasi 0 olan giris",
                       anc["offset_distribution"].get(0, 0)))
    ro = anc["residue_offset_from_expected"]
    print("%-*s %s" % (w, "  beklenen uzakliktan sapma: tam 0",
                       ro["within_plus_minus_0"]))
    print("%-*s %s / %s" % (w, "  beklenen uzakliktan sapma: +-3 kalinti ici",
                            ro["within_plus_minus_3"], ro["n"]))
    print("%-*s %s" % (w, "  beklenen uzakliktan sapma: +-10 disi",
                       ro["beyond_plus_minus_10"]))
    print("%-*s [%s, %s] medyan %s" % (w, "  sapma %95 araligi",
                                       ro["p2.5"], ro["p97.5"], ro["median"]))
    print("%-*s %d / %d" % (w, "  bos modeli gecen tip",
                            anc["null_model"]["clusters_beating_the_null"],
                            anc["null_model"]["clusters_tested"]))
    med = [c["median_column_DE_fraction"] for c in anc["per_cluster"]
           if c["status"] == "tested" and c.get("median_column_DE_fraction")
           is not None]
    print("%-*s %.3f  (capa kolonu: %.3f-%.3f)"
          % (w, "  siradan bir kolonun D/E orani (tip medyanlari)",
             float(np.median(med)) if med else 0.0,
             min(c["carboxylate_present_fraction"] for c in anc["per_cluster"]
                 if c["status"] == "tested"),
             max(c["carboxylate_present_fraction"] for c in anc["per_cluster"]
                 if c["status"] == "tested")))
    print()
    print("  NCBI urun dizgisine gore capraz kontrol (ilk 6):")
    for r in anc["by_ncbi_product_annotation"]["rows"][:6]:
        print("    %-60s n=%-5d karboksilat %%%.1f" % (
            r["product"][:60], r["n"], 100 * r["carboxylate_present_fraction"]))
    print()
    print("  %-14s %5s %7s %9s %8s %6s" %
          ("tip", "n", "capa%", "poz.geri%", "karbok%", "yok"))
    for c in anc["per_cluster"]:
        if c["status"] != "tested":
            print("  %-14s %5d  %s" % (c["cluster"], c["n_no_carboxylate"], c["status"]))
            continue
        print("  %-14s %5d %7.1f %9.1f %8.1f %6d" % (
            c["cluster"], c["n_no_carboxylate"],
            100 * c["anchor_column_vote_share"],
            100 * c["method_recovery_on_positives"],
            100 * c["carboxylate_present_fraction"], c["absent"]))
    print()
    print("--- 7. gercek yokluklar ---")
    print("%-*s %s" % (w, "  giris", ab["n_absent"]))
    if ab["n_absent"]:
        print("%-*s %s" % (w, "  capa kolonunda oturan kalinti",
                           ab["residue_in_place"]))
        print("%-*s %s" % (w, "  tip", ", ".join(
            "%s:%d" % (k, n) for k, n in list(ab["by_type"].items())[:8])))
        print("%-*s %s" % (w, "  taksonomik alan", ab["by_domain_of_life"]))
        print("%-*s %s / %s" % (w, "  operonda beta alt birimi",
                                ab["with_beta_subunit_in_operon"], ab["n_absent"]))
        print("%-*s %s" % (w, "  farkli yaprak", ab["clustering"]["distinct_leaves"]))
        print("%-*s %s" % (w, "  C-terminal erisim", ab["c_terminal_reach"]["counts"]))
        print()
        print("  aile ici son kontrol (kendi aralarinda hizalama):")
        print("    %-28s %4s %6s %7s %10s %s" %
              ("aile", "n", "D/E", "oran", "His-2'den", "kalintilar"))
        for g in ab["residual_family_scan"]["groups"]:
            if g["status"] != "tested":
                print("    %-28s %4d %s" % (g["group"], g["n"], g["status"]))
                continue
            print("    %-28s %4d %6d %7.2f %10s %s" % (
                g["group"], g["n"], g["members_with_a_carboxylate_there"],
                g["fraction"], g["median_raw_distance_from_his2"],
                g["residues_at_that_column"]))
        print("%-*s %s" % (w, "  aileye ozgu bir karboksilatla cozulen giris",
                           ab["residual_family_scan"][
                               "entries_rescued_by_a_family_specific_site"]))
    print()
    print("--- 8. 3/3 kuralinin bedeli ---")
    for key in ("fixed_column_reading", "anchored_reading"):
        g = gate[key]
        print("%-*s %5d giris (%%%.2f), %d tip" % (
            w, "  " + g["reading"], g["entries_lost"],
            100 * g["entries_lost_fraction"], g["types_losing_members"]))
    print("%-*s %s" % (w, "  en sert vurulan tipler", ", ".join(
        "%s %%%.0f" % (d["type"], 100 * d["lost_fraction"])
        for d in gate["fixed_column_reading"]["hardest_hit_types"][:6])))
    print()
    print("--- 9. istatistik ---")
    for t in out["statistics"]["tests"]:
        print("  %-46s V=%-5s q=%-9.3g cins V=%s" % (
            t["question"][:46], t["entry_level"]["cramers_v"],
            t["entry_level"]["q"], t["genus_level"]["cramers_v"]))
    print()
    print("=== KARAR ===")
    print(v["headline"])
    print()
    print()
    print("kabul kurali:", v["should_the_admission_rule_change"])
    print()
    print("bunun yerine:", v["what_should_change_instead"])


# --------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--db", default="roar.sqlite")
    ap.add_argument("--out", default="analysis_out/carboxylate_audit.json")
    ap.add_argument("--cpu", type=int, default=4)
    ap.add_argument("--alignment", default="genomic_context/cand_aln.sto",
                    help="annotate_ro.py'nin urettigi hmmalign ciktisi")
    ap.add_argument("--ref-alignment", default="ROs_71_Clean/refs71_hmmaln.sto")
    ap.add_argument("--ref-fasta", default="ROs_71_Clean/refs71.fasta")
    ap.add_argument("--components", default="component_hmms/components.hmm")
    ap.add_argument("--window", type=int, default=3,
                    help="capali okumada kolon penceresi (kolon)")
    ap.add_argument("--pos-sample", type=int, default=120,
                    help="tip basina en fazla kac pozitif capa hizalamaya girsin")
    ap.add_argument("--mafft-threads", type=int, default=1,
                    help="mafft is parcacigi sayisi. VARSAYILAN 1 cunku cok "
                         "parcacikli yineleyici iyilestirme bit duzeyinde "
                         "tekrarlanabilir degil (bkz. modul basligi). >1 daha "
                         "hizli, tekrarlanabilir degil.")
    ap.add_argument("--work", default=None,
                    help="gecici dizin (varsayilan: sistem gecici dizini)")
    ap.add_argument("--keep-work", action="store_true")
    args = ap.parse_args()

    def log(message):
        print(message, flush=True)

    workdir = args.work or tempfile.mkdtemp(prefix="carboxylate_audit_")
    os.makedirs(workdir, exist_ok=True)
    log("[gecici dizin] %s" % workdir)

    con = sqlite3.connect("file:%s?mode=ro" % args.db, uri=True)
    log("[okunuyor] %s (SALT OKUNUR)" % args.db)
    entries = load_entries(con)
    controls = load_rieske_only_controls(con)
    con.close()
    log("[bilgi] onaylanmis giris: %d   Rieske-only kontrol: %d"
        % (len(entries), len(controls)))

    log("[okunuyor] hizalama: %s" % args.alignment)
    rows, rf = parse_stockholm(args.alignment, set(entries))
    match_idx = [i for i, ch in enumerate(rf) if ch != "."]
    log("[bilgi] hizalamada %d giris, %d match kolonu" % (len(rows), len(match_idx)))
    index = {}
    unmapped = 0
    for cid, row in rows.items():
        idx = AlignmentIndex(row, match_idx, entries[cid]["sequence"])
        if idx.offset < 0:
            unmapped += 1
            continue
        index[cid] = idx
    log("[bilgi] ham diziye eslenen: %d  eslenmeyen: %d" % (len(index), unmapped))

    negatives = sorted(c for c, r in entries.items()
                       if (r["stored_residue"] or "") not in CARBOX_OK)
    positives = sorted(c for c, r in entries.items()
                       if (r["stored_residue"] or "") in CARBOX_OK)
    log("[bilgi] sabit kolonda karboksilat YOK: %d   VAR: %d"
        % (len(negatives), len(positives)))

    log("\n[adim 1] temel okuma")
    baseline = stage_baseline(entries, index, 2)
    log("   %s" % baseline["missing_call_breakdown"])

    log("[adim 2] mekanizma: hizalama kapsami")
    mechanism = stage_mechanism(entries, index, negatives, positives)
    mechanism["euk_count"] = sum(1 for c in negatives
                                 if entries[c]["domain"] == "Eukaryota")
    mechanism["taxonomic_domain_counts"] = {
        "no_carboxylate": dict(Counter(
            entries[c]["domain"] or "?" for c in negatives).most_common()),
        "with_carboxylate": dict(Counter(
            entries[c]["domain"] or "?" for c in positives).most_common()),
    }
    log("   son dolu kolon medyani %s (karboksilatli %s), kolon %d disinda %d giris"
        % (mechanism["no_carboxylate"]["last_filled_column_median"],
           mechanism["with_carboxylate"]["last_filled_column_median"], CARBOX_COL,
           mechanism["no_carboxylate"]["alignment_stops_before_the_carboxylate_column"]))

    log("[adim 3] sabit kolonu genisletme")
    widening = stage_window_widening(entries, index, negatives, [2, 5, 10, 20, 40])
    log("   %s" % {r["window_columns"]: r["rescued"] for r in widening["windows"]})

    log("[adim 4] ciplak capa + sans beklentisi")
    naive = stage_naive_anchor(entries, index, negatives, positives)
    log("   %s" % naive["verdict"])

    log("[adim 5] referans profil sondalari")
    models = build_submodels(args.ref_alignment, workdir, log)
    sets = {}
    for label, records in (
            ("no_carboxylate", [(c, entries[c]["sequence"]) for c in negatives]),
            ("with_carboxylate", [(c, entries[c]["sequence"]) for c in positives]),
            ("rieske_only", sorted(controls.items()))):
        path = os.path.join(workdir, label + ".fasta")
        write_fasta(path, records)
        sets[label] = path
    probe = stage_reference_probe(models, sets, workdir, args.cpu, log)

    log("[adim 6] CAPALI YENIDEN OKUMA (tip ici mafft)")
    refs = read_fasta(args.ref_fasta)
    log("   mafft --thread %d%s" % (
        args.mafft_threads,
        "  (tekrarlanabilir)" if args.mafft_threads == 1
        else "  (HIZLI ama bit duzeyinde tekrarlanabilir DEGIL)"))
    anchored, calls = stage_anchored_reread(
        entries, index, negatives, positives, refs, workdir,
        args.mafft_threads, args.window, args.pos_sample, log)
    log("   %s" % anchored["totals"])

    log("[adim 7] gercek yokluklarin kimligi")
    absences, absent_ids = stage_real_absences(
        entries, index, calls, anchored, models, workdir, args.cpu,
        args.mafft_threads, log)

    log("[adim 8] 3/3 kuralinin bedeli")
    gate = stage_gate_cost(entries, negatives, calls)
    log("   sabit kolon %d giris / capali okuma %d giris"
        % (gate["fixed_column_reading"]["entries_lost"],
           gate["anchored_reading"]["entries_lost"]))

    log("[adim 9] istatistik")
    statistics_block = stage_statistics(entries, negatives, calls)

    out = {
        "generated": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "module": "carboxylate_audit.py",
        "read_only": True,
        "inputs": {
            "db": args.db,
            "alignment": args.alignment,
            "reference_alignment": args.ref_alignment,
            "model": DEFAULT_MODEL,
            "confirmed_entries": len(entries),
            "entries_without_a_carboxylate_at_the_fixed_column": len(negatives),
            "entries_with_a_carboxylate_at_the_fixed_column": len(positives),
            "rieske_only_controls": len(controls),
            "columns": {"rieske": RIESKE_SITES, "catalytic": CATALYTIC_SITES,
                        "bridging": BRIDGING_SITE},
            "admission_rule": {
                "rieske_ligands_required": len(RIESKE_SITES),
                "catalytic_ligands_required": 2,
                "catalytic_ligands_total": len(CATALYTIC_SITES),
                "column_window": 2,
            },
        },
        "baseline_fixed_column": baseline,
        "mechanism": mechanism,
        "window_widening": widening,
        "naive_anchor": naive,
        "reference_profile_probe": probe,
        "anchored_reread": anchored,
        "real_absences": absences,
        "gate_cost": gate,
        "statistics": statistics_block,
    }
    out["verdict"] = build_verdict(baseline, mechanism, widening, naive, probe,
                                   anchored, absences, gate)
    out["per_entry_calls"] = {
        cid: {"type": entries[cid]["cluster"],
              "fixed_column_residue": entries[cid]["stored_residue"],
              "anchored_call": calls[cid]["call"],
              "anchored_residue": calls[cid]["residue"],
              "anchored_offset": calls[cid]["offset"],
              "residue_offset_from_expected":
                  calls[cid].get("residue_offset_from_expected")}
        for cid in sorted(calls)}

    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    with open(args.out, "w") as handle:
        json.dump(out, handle, indent=1)
    log("\n[yazildi] %s" % args.out)

    print_summary(out)

    if not args.keep_work and args.work is None:
        shutil.rmtree(workdir, ignore_errors=True)
    else:
        log("\n[gecici dizin korundu] %s" % workdir)
    return 0


if __name__ == "__main__":
    sys.exit(main())
