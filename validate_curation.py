#!/usr/bin/env python3
"""
Kuratorlu veri ve veritabani dogrulayicisi -- run_all.sh icinde KAPI olarak kullanilir.

NEDEN: bu pipeline birkac SESSIZ hata yayinladi. cd-hit cokmesi yutuldu, bir CSV
virgul kacisi ('2\\,4') kolonlari kaydirdi ve TftA/DntAc satirlari yanlis yere
dustu, bir metrik en kotu hizalanmis ilmekleri "ayirt edici kolon" sandi,
hmmsearch hiti olmayan 945 aday is_confirmed=NULL kaldi. Hepsinin ortak yani:
hata ciktiyi BOZDU ama hicbir sey durmadi.

BU SCRIPT her kontrol icin tek satir PASS/FAIL ve bir sayi yazar, sonunda
basarisizlik varsa SIFIR-DISI cikar. Boylece run_all.sh'de bir adim olarak
durabilir ve bozuk veri HTML'e / web sitesine kadar gidemez.

Kullanim:
    python3 validate_curation.py
    python3 validate_curation.py --db roar.sqlite --ecology cluster_ecology.csv \\
        --chemistry chemistry.csv
    python3 validate_curation.py --warn-only   # rapor yaz ama her zaman 0 don

Cikis kodu: 0 = tum FAIL kontrolleri gecti, 1 = en az bir FAIL, 2 = kontrol
kosulamadi (veritabani/CSV yok).

SINIR: burada gecen her sey DOGRU demek degildir; bkz. README_PIPELINE.md
"Safeguards" bolumu -- bu dogrulayicinin yakalayamadiklari orada yaziyor.
"""

import argparse
import csv
import os
import re
import sqlite3
import sys

ROOT = os.path.dirname(os.path.abspath(__file__))

# Kuratorlu referans seti: pipeline'in TAMAMI bunun uzerine kuruluyor
# (tip atamasi, katalitik motif kolonlari, her kanit kademesi). Buradaki bir
# kusur tek bir sayiyi degil, butun siniflandirmayi kaydirir.
REFS_FASTA = os.path.join(ROOT, "ROs_71_Clean", "refs71.fasta")

# Iki referansin birbirine bu orandan fazla benzemesi RAPOR EDILIR ama hata
# sayilmaz: EdoA1/cumA1 %99,8 kimlikte FARKLI substratlara (etilbenzen /
# kumen) etki ediyor ve bu veritabaninin merkezi bulgusudur. Ayni sey
# DUPLIKASYON kusuruyla karistirilmamali -- o ayri ve sert bir kontrol.
NEAR_IDENTICAL_PCT = 99.0

# --- Izinli sozluk degerleri -------------------------------------------------
# add_reference.py yeni referans eklerken AYNI listeleri kullanir; ikisi
# ayrismasin diye asagida ayrica karsilastirilirlar.
VALID_SUBSTRATE_CLASS = {"xenobiotic", "natural_aromatic", "natural_specialized",
                         "unknown"}
VALID_CONFIDENCE = {"low", "medium", "high"}
VALID_REACTION_CLASS = {"cis_dihydroxylation", "angular_dioxygenation",
                        "dioxygenation_with_release", "O_demethylation",
                        "N_demethylation", "hydroxylation", "C_N_cleavage",
                        "unknown"}
VALID_FAMILY = {"alkylbenzenes", "pah", "biaryls_ethers", "nitroaromatics",
                "haloaromatics", "sulfoaromatics", "aromatic_acids", "anilines",
                "quaternary_amines", "alkaloids", "terpenoids_steroids", "unknown"}
VALID_SOURCE_KIND = {"paper", "structure", "curator", "paper+structure",
                     "structure+paper", "paper+curator", "general", "none"}

# --- Turkce sizintisi --------------------------------------------------------
# Web sitesini goren hicbir metin Turkce olmamali. Iki ayri tuzak var:
# (a) Turkce'ye ozgu harfler, (b) ASCII'ye sadelestirilmis Turkce kelimeler --
# bunlar daha sinsi, cunku kodlama kontrolleri onlari gormez.
# Turkce'ye ozgu harfler, KOD SAF ASCII kalsin diye kod noktasiyla yazildi:
# c-cedilla, g-breve, dotless-i, o-umlaut, s-cedilla, u-umlaut (kucuk+buyuk).
TURKISH_LETTERS = {chr(c) for c in (
    0x00E7, 0x011F, 0x0131, 0x00F6, 0x015F, 0x00FC,
    0x00C7, 0x011E, 0x0130, 0x00D6, 0x015E, 0x00DC,
)}
TURKISH_WORDS = ["agirlikli", "komsu", "kume", "yaprak", "bilinmiyor",
                 "dusuk", "orta", "yuksek", "degil"]
TURKISH_WORD_RE = re.compile(r"\b(" + "|".join(TURKISH_WORDS) + r")", re.IGNORECASE)

# Kelime arama ONEK esleseme yapar (Turkce eklemeli: komsu -> komsular,
# kume -> kumeler) ve bu bazen mesru metni yakalar. Bilinen yanlis pozitifler
# burada muaf tutulur; liste KISA kalmali, yoksa kontrolun degeri duser.
#   "Kumeu" : Yeni Zelanda'da bir yer adi, GenBank /geo_loc_name alanindan
#             geliyor ve 'kume' onekiyle eslesiyor.
LANGUAGE_FALSE_POSITIVES = ("Kumeu",)

# Kontrol edilecek kullanici-yuzlu DB alanlari: (tablo, anahtar kolon, metin kolonu)
USER_FACING_DB_FIELDS = [
    ("leaf_profile", "leaf_id", "label"),
    ("leaf", "leaf_id", "top_genera"),
    ("ro_etc", "candidate_id", "etc_profile"),
    ("operon", "candidate_id", "layout"),
    # Proje tarafindan URETILEN habitat etiketi (serbest metin degil, sozluk).
    # isolation_source.py'nin HABITAT_RULES sonucu; web'de gosteriliyor.
    ("replicon_source", "nucleotide_id", "habitat"),
]

# analysis_out/ icinde web arayuzunun DOGRUDAN yayinladigi metin dosyalari.
# Turkce sizintisi burada da yayina gider; bu yuzden ayrica taraniyor.
PUBLISHED_TEXT_FILES = [
    "cluster_ecology_stats.csv", "reference_pairs.csv", "representatives.csv",
    "etc_by_cluster.csv", "evidence_by_cluster.csv", "regulation_by_cluster.csv",
    "leaf_profiles.csv", "domain_by_cluster.csv", "null_model.csv",
    "cluster_variance.csv", "variant_signatures.csv", "leaves.csv",
    "novel_candidates.csv", "transposon_by_genus.csv", "sdp_positions.csv",
    "subfamilies.csv", "cluster_identity_matrix.csv", "ssn_nodes.csv",
    "habitat.json", "motif_stats.json", "stats.json", "redundancy.json",
    "cooccurrence.json", "operon_validation.json",
    "substrate_predictability.json",
]

# --- Uyesiz referans tipleri -------------------------------------------------
# cluster_ecology.csv ve chemistry.csv 71 kuratorlu referansin HEPSINI tutar,
# ama bu 10 tipin veritabaninda TEK bir adayi bile yok (olculdu: hicbir
# PF00355 proteini bu modellere en yuksek bit skorunu vermiyor; sadece
# dogrulanmamis degil, hic atanmamis). Yani "CSV'de var, ro'da yok" durumu
# bunlar icin BEKLENEN; listede olmayan yeni bir bosluk ise HATA'dir, cunku
# o, bir tipin uyelerini kaybettigi anlamina gelir.
#
# DIKKAT: ayni 10 kimlik check_reference_set() icinde de WARN olarak cikar,
# ama farkli bir soru icin. Burada sorulan "CSV ile DB ayristi mi"; orada
# sorulan "bu referans NEDEN uye toplamiyor". Ucu (3_309_NahAc, 3_315_NDO,
# 1_115_NdmC) referans setindeki duplikasyon/parca kusuru yuzunden bos;
# kusur duzeltilince bu listenin de kisalmasi gerekir.
REFERENCE_ONLY_CLUSTERS = {
    "1_103_CarAa", "1_115_NdmC", "2_206_IpBAa", "2_207_BpdC1", "2_210_XylC1",
    "3_309_NahAc", "3_310_DntAc", "3_312_NinA", "3_315_NDO", "4_407_PbaA1",
}

# SMILES'te mesru olan karakterler. Atom, bag, halka kapanisi, yuk, kiralite.
SMILES_LEGAL = set("ABCDEFGHIKLMNOPRSTUVWXYZabcdefghiklmnoprstuvy"
                   "0123456789[]()=#$:/\\@+-.%*")
SMILES_ATOM_RE = re.compile(r"\[[^\]]*\]|Cl|Br|[BCNOPSFIbcnopsfi]")


class Report:
    """PASS/FAIL satirlari ve cikis kodu."""

    def __init__(self):
        self.rows = []
        self.failed = 0

    def check(self, name, ok, count, detail="", severity="FAIL"):
        """severity='FAIL' kontrol kapiyi kapatir; 'WARN' sadece rapor eder."""
        status = "PASS" if ok else severity
        if not ok and severity == "FAIL":
            self.failed += 1
        self.rows.append((status, name, count, detail))
        return ok

    def emit(self):
        width = max(len(r[1]) for r in self.rows) if self.rows else 10
        for status, name, count, detail in self.rows:
            line = f"{status:4s}  {name:<{width}s}  n={count}"
            if detail:
                line += f"   {detail}"
            print(line)
        n_pass = sum(1 for r in self.rows if r[0] == "PASS")
        n_warn = sum(1 for r in self.rows if r[0] == "WARN")
        print("-" * (width + 24))
        print(f"      {len(self.rows)} checks: {n_pass} PASS, "
              f"{self.failed} FAIL, {n_warn} WARN")
        return 1 if self.failed else 0


def first(items, n=6):
    items = sorted(str(i) for i in items)
    out = ", ".join(items[:n])
    if len(items) > n:
        out += f", ... (+{len(items) - n})"
    return out


def turkish_hits(text):
    """Metinde Turkce harf veya Turkce kelime var mi -> bulgu listesi."""
    if text is None:
        return []
    text = str(text)
    for allowed in LANGUAGE_FALSE_POSITIVES:
        text = text.replace(allowed, "")
    hits = sorted({ch for ch in text if ch in TURKISH_LETTERS})
    hits += sorted({m.group(1).lower() for m in TURKISH_WORD_RE.finditer(text)})
    return hits


# ---------------------------------------------------------------------------
# SMILES makul mu
# ---------------------------------------------------------------------------

def smiles_problems(smiles):
    """Makul bir SMILES mi? Kimyasal dogruluk DEGIL, bicim dogrulugu.

    Kontroller: parantez ve koseli parantez dengesi, yalnizca mesru atom/bag
    karakterleri, en az bir atom, halka kapanis rakamlarinin eslesmesi.
    Burada gecmek molekulun DOGRU oldugu anlamina gelmez; yalnizca cizim
    kutuphanesine verildiginde kirilmayacagi anlamina gelir.
    """
    problems = []
    for opener, closer in (("(", ")"), ("[", "]")):
        depth = 0
        for ch in smiles:
            if ch == opener:
                depth += 1
            elif ch == closer:
                depth -= 1
                if depth < 0:
                    problems.append(f"'{closer}' before '{opener}'")
                    break
        if depth > 0:
            problems.append(f"{depth} unclosed '{opener}'")
    illegal = sorted({ch for ch in smiles if ch not in SMILES_LEGAL})
    if illegal:
        problems.append("illegal characters " + repr("".join(illegal)))
    if not SMILES_ATOM_RE.search(smiles):
        problems.append("no atom found")
    # Halka kapanis rakamlari: koseli parantez icindekiler (yuk, izotop) sayilmaz
    core = re.sub(r"\[[^\]]*\]", "A", smiles)
    digits = re.findall(r"%\d\d|\d", core)
    unpaired = sorted({d for d in set(digits) if digits.count(d) % 2})
    if unpaired:
        problems.append(f"unpaired ring-closure digits {unpaired}")
    return problems


# ---------------------------------------------------------------------------
# CSV okuma: kolon sayisi kaymasi (virgul kacisi hata sinifi)
# ---------------------------------------------------------------------------

def read_csv_strict(path, rep, label):
    """CSV'yi oku ve basliktan FARKLI alan sayisi olan satirlari bildir.

    2\\,4 kacisi hata sinifi tam olarak bu: satir sessizce fazla/eksik alanla
    ayristirilir, degerler bir kolon kayar ve hicbir sey uyarmaz.
    """
    with open(path, newline="", encoding="utf-8") as fh:
        rows = list(csv.reader(fh))
    if not rows:
        rep.check(f"{label}: file not empty", False, 0, "file has no rows")
        return [], []
    header = rows[0]
    bad = [(i + 1, len(r)) for i, r in enumerate(rows) if r and len(r) != len(header)]
    rep.check(f"{label}: field count matches header", not bad, len(bad),
              "" if not bad else f"{len(header)} header fields; offending "
                                 f"lines/counts: {first(bad)}")
    dicts = [dict(zip(header, r)) for r in rows[1:] if r]
    return header, dicts


# ---------------------------------------------------------------------------
# Kontroller
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# Kuratorlu referans setinin kendisi
# ---------------------------------------------------------------------------

def read_fasta(path):
    """[(header, sequence)] dondur. Dizi buyuk harfe cevrilir, bosluk atilir."""
    records, header, chunks = [], None, []
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(chunks)))
                header, chunks = line[1:].strip(), []
            elif line:
                chunks.append(line.upper().replace(" ", ""))
    if header is not None:
        records.append((header, "".join(chunks)))
    # Bazi kaynaklar dizinin sonuna '*' koyar; karsilastirmayi bozmasin.
    return [(h, s.rstrip("*")) for h, s in records]


def ref_id(header):
    """FASTA basligindan kume kimligi: evidence_tiers.py ile AYNI kural.

    Baslik '1_101_OxoO_monooxygenase_pro' bicimindedir; ilk uc alt-cizgi
    parcasi kume kimligidir. Kural ayrisirsa dogrulama baska bir sey sinar.
    """
    return "_".join(header.split("_")[:3])


def check_reference_set(rep, refs_path, pairs_path, eco, chem, con):
    """Referans setinin IC kusurlari.

    Burada aranan hatalar 2026-10 incelemesinde ELLE bulundu ve hepsi ayni
    sonucu dogurdu: RieskeDB71.hmm'e birbirinin AYNISI profiller girdi, bu
    yuzden tip atamasi ikizler arasinda keyfi boluldu (OxoO 2 / OMO 11,
    NDO 2 / NarAa 30) ve "en yuksek skorlu profil ile en yakin referans
    protein farkli tipte" orani (%29,3) buyuk olcude bundan sisti.
    """
    if not os.path.exists(refs_path):
        rep.check("refs71: reference FASTA present", False, 0,
                  f"{refs_path} not found; the whole reference-set check could "
                  f"not run")
        return

    records = read_fasta(refs_path)
    ids = [ref_id(h) for h, _ in records]
    seqs = [s for _, s in records]
    rep.check("refs71: FASTA is readable and non-empty", bool(records),
              len(records), f"{len(set(seqs))} distinct sequences")

    # 1. Birbirinin AYNISI diziler. Ayni enzimin iki adla girilmesi
    #    HMM kutuphanesine cift profil koyar.
    by_seq = {}
    for ident, seq in zip(ids, seqs):
        by_seq.setdefault(seq, []).append(ident)
    dup_groups = [g for g in by_seq.values() if len(g) > 1]
    detail = "; ".join(" == ".join(sorted(g)) for g in dup_groups)
    rep.check("refs71: no two reference sequences are identical",
              not dup_groups, len(dup_groups),
              f"{len(records)} records, {len(by_seq)} distinct sequences"
              if not dup_groups else
              f"{len(records)} records but only {len(by_seq)} distinct "
              f"sequences -> {detail}")

    # 2. Bir dizinin bir baskasinin ICINDE gecmesi: ayni protein, bir kopyasi
    #    kirpilmis. Parcadan kurulan profil sistematik olarak zayif olur.
    duplicated = {i for g in dup_groups for i in g}
    substrings = []
    for ident_a, seq_a in zip(ids, seqs):
        if not seq_a:
            continue
        for ident_b, seq_b in zip(ids, seqs):
            if ident_a == ident_b or seq_a == seq_b:
                continue
            if seq_a in seq_b:
                substrings.append(f"{ident_a} ({len(seq_a)} aa) is contained "
                                  f"in {ident_b} ({len(seq_b)} aa)")
    rep.check("refs71: no reference sequence is a substring of another",
              not substrings, len(substrings), first(substrings, 4))

    # 3. Her kimlik UC yerde de tam olarak bir kez: FASTA, cluster_ecology.csv,
    #    chemistry.csv. Biri eksikse web arayuzu o tipi yarim gosterir.
    fasta_counts = {}
    for ident in ids:
        fasta_counts[ident] = fasta_counts.get(ident, 0) + 1
    eco_counts, chem_counts = {}, {}
    for row in eco:
        key = row.get("cluster", "")
        eco_counts[key] = eco_counts.get(key, 0) + 1
    for row in chem:
        key = row.get("cluster", "")
        chem_counts[key] = chem_counts.get(key, 0) + 1
    every = set(fasta_counts) | set(eco_counts) | set(chem_counts)
    bad = []
    for ident in sorted(every):
        counts = (fasta_counts.get(ident, 0), eco_counts.get(ident, 0),
                  chem_counts.get(ident, 0))
        if counts != (1, 1, 1):
            bad.append(f"{ident}: fasta={counts[0]} ecology={counts[1]} "
                       f"chemistry={counts[2]}")
    rep.check("refs71: every reference id appears exactly once in FASTA, "
              "ecology and chemistry", not bad, len(bad),
              f"{len(every)} reference ids checked" if not bad
              else first(bad, 4))

    # 4. Neredeyse ayni ama AYNI OLMAYAN ciftler: WARN. Bu bir kusur DEGIL,
    #    bu veritabaninin bulgusu. Kimlik olcumu pipeline'in kendi diamond
    #    sonucundan (reference_pairs.csv) okunur, yeniden hesaplanmaz.
    chem_by = {r["cluster"]: r for r in chem}

    def substrate_of(ident):
        return (chem_by.get(ident, {}).get("substrate_en") or "unknown").strip()

    if not os.path.exists(pairs_path):
        rep.check(f"refs71: near-identical pairs above "
                  f"{NEAR_IDENTICAL_PCT:.0f}% identity", True, 0,
                  f"{os.path.basename(pairs_path)} not present; run "
                  f"evidence_tiers.py -- check skipped", severity="WARN")
    else:
        # Kimlik olcumu TURETILMIS bir dosyadan okunuyor. O dosya FASTA'dan
        # eskiyse olcum artik baska bir referans setini tarif ediyor; uyari
        # yaniltici olur, o yuzden acikca yazilir.
        stale_note = ""
        try:
            if os.path.getmtime(pairs_path) < os.path.getmtime(refs_path):
                stale_note = (" [WARNING: reference_pairs.csv is OLDER than "
                              "refs71.fasta, so these identities describe an "
                              "earlier reference set; re-run evidence_tiers.py]")
        except OSError:
            pass
        near = []
        with open(pairs_path, newline="", encoding="utf-8") as fh:
            for row in csv.DictReader(fh):
                try:
                    pid = float(row.get("identity", ""))
                except ValueError:
                    continue
                if pid <= NEAR_IDENTICAL_PCT:
                    continue
                a, b = row.get("ref_a", ""), row.get("ref_b", "")
                # Tam kopya / parca ciftleri yukarida SERT kontrol olarak
                # raporlandi; burada tekrar edilmemeleri gerekiyor, yoksa
                # gercek biyolojik bulgu duplikasyon gurultusunda kaybolur.
                if a in duplicated and b in duplicated:
                    continue
                if any(f"{a} (" in s and f"in {b} (" in s
                       or f"{b} (" in s and f"in {a} (" in s
                       for s in substrings):
                    continue
                near.append(f"{a} / {b} at {pid:.1f}% "
                            f"({substrate_of(a)} against {substrate_of(b)})")
        rep.check(f"refs71: near-identical pairs above "
                  f"{NEAR_IDENTICAL_PCT:.0f}% identity", not near, len(near),
                  "; ".join(near[:4]) + (f"; ... (+{len(near) - 4})"
                                         if len(near) > 4 else "") + stale_note,
                  severity="WARN")

    # 5. Hic uye toplamayan referanslar: WARN. Sifir uye hem ikizin hem de
    #    parca referansin imzasidir, o yuzden isimleriyle yazilir.
    members = {r[0]: r[1] for r in con.execute(
        "SELECT ro_cluster, COUNT(*) FROM ro WHERE is_confirmed=1 "
        "AND ro_cluster IS NOT NULL GROUP BY ro_cluster")}
    zero = [ident for ident in sorted(fasta_counts) if members.get(ident, 0) == 0]
    notes = []
    for ident in zero:
        tag = ""
        if ident in duplicated:
            twin = [i for g in dup_groups if ident in g for i in g if i != ident]
            tag = f" [identical twin of {','.join(twin)}]"
        elif any(s.startswith(f"{ident} (") for s in substrings):
            tag = " [fragment of another reference]"
        notes.append(ident + tag)
    rep.check("refs71: every reference recruits at least one confirmed member",
              not zero, len(zero),
              f"{len(fasta_counts) - len(zero)} of {len(fasta_counts)} "
              f"references have members; zero-member: " + first(notes, 12),
              severity="WARN")


def check_vocabulary_drift(rep):
    """add_reference.py yazarken BU dosyanin dogrularken kullandigi listeler ayni mi?

    Yeni referans ekleyen script ile dogrulayici ayrisirsa, eklenen satir
    dogrulamadan gecer ama web arayuzu onu taniyamaz.
    """
    try:
        sys.path.insert(0, ROOT)
        import add_reference as ar          # noqa: PLC0415
    except Exception:                       # noqa: BLE001
        rep.check("vocabulary: writer and validator agree", True, 0,
                  "add_reference.py not importable; comparison skipped",
                  severity="WARN")
        return
    pairs = [
        ("substrate_class", VALID_SUBSTRATE_CLASS, getattr(ar, "VALID_SUBSTRATE_CLASS", None)),
        ("confidence", VALID_CONFIDENCE, getattr(ar, "VALID_CONFIDENCE", None)),
        ("reaction_class", VALID_REACTION_CLASS, getattr(ar, "VALID_REACTION_CLASS", None)),
        ("family", VALID_FAMILY, getattr(ar, "VALID_FAMILY", None)),
        ("source_kind", VALID_SOURCE_KIND, getattr(ar, "VALID_SOURCE_KIND", None)),
    ]
    diffs = []
    for field, mine, theirs in pairs:
        if theirs is None:
            diffs.append(f"{field}: missing in add_reference.py")
        elif set(theirs) != mine:
            diffs.append(f"{field}: {sorted(set(theirs) ^ mine)}")
    rep.check("vocabulary: writer and validator agree", not diffs, len(diffs),
              first(diffs) if diffs else "")


def check_confirmed_ro(con, rep):
    n_conf = con.execute("SELECT COUNT(*) FROM ro WHERE is_confirmed=1").fetchone()[0]
    rep.check("ro: confirmed entries exist", n_conf > 0, n_conf)

    bad = con.execute("""
        SELECT candidate_id FROM ro WHERE is_confirmed=1
          AND (ro_cluster IS NULL OR TRIM(ro_cluster)='' OR ro_cluster='N/A')
        LIMIT 20""").fetchall()
    n = con.execute("""
        SELECT COUNT(*) FROM ro WHERE is_confirmed=1
          AND (ro_cluster IS NULL OR TRIM(ro_cluster)='' OR ro_cluster='N/A')
        """).fetchone()[0]
    rep.check("ro: confirmed has a real ro_cluster", n == 0, n,
              "" if n == 0 else f"e.g. {first(r[0] for r in bad)}")

    n = con.execute("""
        SELECT COUNT(*) FROM ro WHERE is_confirmed=1
          AND (sequence IS NULL OR TRIM(sequence)='')""").fetchone()[0]
    rep.check("ro: confirmed has a sequence", n == 0, n)

    n = con.execute("""
        SELECT COUNT(*) FROM ro WHERE is_confirmed=1
          AND (rieske_intact IS NULL OR rieske_intact!=1)""").fetchone()[0]
    rep.check("ro: confirmed has rieske_intact=1", n == 0, n)

    # is_confirmed'in NULL kalmasi daha once gercek bir hataydi (945 aday).
    n = con.execute("SELECT COUNT(*) FROM ro WHERE is_confirmed IS NULL").fetchone()[0]
    rep.check("ro: is_confirmed is never NULL", n == 0, n,
              "" if n == 0 else "candidates without an hmmsearch hit must be 0, "
                                "not NULL")
    return n_conf


def check_cluster_coverage(con, rep, eco, chem):
    db_clusters = {r[0] for r in con.execute(
        "SELECT DISTINCT ro_cluster FROM ro WHERE is_confirmed=1 "
        "AND ro_cluster IS NOT NULL")}
    eco_clusters = {r["cluster"] for r in eco}
    chem_clusters = {r["cluster"] for r in chem}

    missing_eco = db_clusters - eco_clusters
    rep.check("coverage: every ro cluster is in cluster_ecology.csv",
              not missing_eco, len(missing_eco), first(missing_eco))
    missing_chem = db_clusters - chem_clusters
    rep.check("coverage: every ro cluster is in chemistry.csv",
              not missing_chem, len(missing_chem), first(missing_chem))

    # Ters yon: CSV'de olup veritabaninda uyesi olmayan tipler. Bilinen 10
    # uyesiz referans muaf; baska bir bosluk cikarsa bir tip uyelerini
    # kaybetmis demektir.
    orphan_eco = eco_clusters - db_clusters - REFERENCE_ONLY_CLUSTERS
    rep.check("coverage: cluster_ecology.csv has no unexplained extra type",
              not orphan_eco, len(orphan_eco), first(orphan_eco))
    orphan_chem = chem_clusters - db_clusters - REFERENCE_ONLY_CLUSTERS
    rep.check("coverage: chemistry.csv has no unexplained extra type",
              not orphan_chem, len(orphan_chem), first(orphan_chem))

    # Muafiyet listesi bayatladi mi: listedeki bir tip uye kazanmissa
    # listeden cikarilmasi gerekir.
    now_populated = REFERENCE_ONLY_CLUSTERS & db_clusters
    rep.check("coverage: REFERENCE_ONLY_CLUSTERS list is still accurate",
              not now_populated, len(now_populated),
              "" if not now_populated else
              f"these now have members, remove them: {first(now_populated)}",
              severity="WARN")

    # Iki CSV birbirini kapsiyor mu
    only_eco = eco_clusters - chem_clusters
    only_chem = chem_clusters - eco_clusters
    rep.check("coverage: the two curated CSVs cover the same types",
              not only_eco and not only_chem, len(only_eco) + len(only_chem),
              "" if not (only_eco or only_chem) else
              f"ecology-only: {first(only_eco)} | chemistry-only: {first(only_chem)}")

    dup_eco = len(eco) - len(eco_clusters)
    dup_chem = len(chem) - len(chem_clusters)
    rep.check("coverage: no duplicate cluster rows in the CSVs",
              dup_eco == 0 and dup_chem == 0, dup_eco + dup_chem)
    return db_clusters


def check_ecology_csv(rep, eco):
    bad = [(r["cluster"], r["substrate_class"]) for r in eco
           if r.get("substrate_class") not in VALID_SUBSTRATE_CLASS]
    rep.check("cluster_ecology.csv: substrate_class in vocabulary",
              not bad, len(bad), first(f"{c}={v!r}" for c, v in bad))
    bad = [(r["cluster"], r["confidence"]) for r in eco
           if r.get("confidence") not in VALID_CONFIDENCE]
    rep.check("cluster_ecology.csv: confidence in vocabulary",
              not bad, len(bad), first(f"{c}={v!r}" for c, v in bad))
    bad = [r for r in eco if not (r.get("cluster") or "").strip()
           or not (r.get("substrate") or "").strip()]
    rep.check("cluster_ecology.csv: cluster and substrate are non-empty",
              not bad, len(bad))


def check_chemistry_csv(rep, chem):
    bad = [(r["cluster"], r["reaction_class"]) for r in chem
           if r.get("reaction_class") not in VALID_REACTION_CLASS]
    rep.check("chemistry.csv: reaction_class in vocabulary",
              not bad, len(bad), first(f"{c}={v!r}" for c, v in bad))
    bad = [(r["cluster"], r["family"]) for r in chem
           if r.get("family") not in VALID_FAMILY]
    rep.check("chemistry.csv: family in vocabulary",
              not bad, len(bad), first(f"{c}={v!r}" for c, v in bad))
    bad = [(r["cluster"], r["source_kind"]) for r in chem
           if r.get("source_kind") not in VALID_SOURCE_KIND]
    rep.check("chemistry.csv: source_kind in vocabulary",
              not bad, len(bad), first(f"{c}={v!r}" for c, v in bad))
    bad = [(r["cluster"], r.get("curation_confidence")) for r in chem
           if r.get("curation_confidence") not in VALID_CONFIDENCE]
    rep.check("chemistry.csv: curation_confidence in vocabulary",
              not bad, len(bad), first(f"{c}={v!r}" for c, v in bad))

    bad = []
    n_smiles = 0
    for r in chem:
        s = (r.get("substrate_smiles") or "").strip()
        if not s:
            continue
        n_smiles += 1
        problems = smiles_problems(s)
        if problems:
            bad.append(f"{r['cluster']}: {'; '.join(problems)}")
    rep.check("chemistry.csv: every SMILES is parseable",
              not bad, len(bad),
              f"{n_smiles} non-empty SMILES checked" if not bad else first(bad, 4))

    bad = [(r["cluster"], r["pdb"]) for r in chem
           if (r.get("pdb") or "").strip()
           and not re.fullmatch(r"[A-Za-z0-9]{4}", r["pdb"].strip())]
    n_pdb = sum(1 for r in chem if (r.get("pdb") or "").strip())
    rep.check("chemistry.csv: every pdb is four alphanumerics",
              not bad, len(bad),
              f"{n_pdb} non-empty pdb ids checked" if not bad
              else first(f"{c}={v!r}" for c, v in bad))

    # Kaynak beyan edilmeden substrat verilmesi: web sayfasi substrati
    # gosterir ama nereden geldigini soyleyemez.
    bad = [r["cluster"] for r in chem
           if (r.get("substrate_en") or "").strip() not in ("", "unknown")
           and not (r.get("source") or "").strip()]
    rep.check("chemistry.csv: a named substrate always cites a source",
              not bad, len(bad), first(bad))


def check_turkish(con, rep, eco_path, chem_path, eco, chem, out_dir):
    total = 0
    for table, key, col in USER_FACING_DB_FIELDS:
        if not con.execute("SELECT name FROM sqlite_master WHERE name=?",
                           (table,)).fetchone():
            rep.check(f"language: {table}.{col} is English", True, 0,
                      "table not in this database; check skipped",
                      severity="WARN")
            continue
        try:
            rows = con.execute(
                f'SELECT "{key}", "{col}" FROM "{table}" '
                f'WHERE "{col}" IS NOT NULL').fetchall()
        except sqlite3.DatabaseError as exc:
            rep.check(f"language: {table}.{col} is English", False, 0, str(exc))
            continue
        bad = []
        for ident, text in rows:
            hits = turkish_hits(text)
            if hits:
                bad.append(f"{ident}: {first(hits, 3)} in {str(text)[:40]!r}")
        total += len(bad)
        rep.check(f"language: {table}.{col} is English", not bad, len(bad),
                  f"{len(rows)} values checked" if not bad else first(bad, 3))

    for path, rows, label in ((eco_path, eco, "cluster_ecology.csv"),
                              (chem_path, chem, "chemistry.csv")):
        bad = []
        for i, row in enumerate(rows, start=2):
            for field, value in row.items():
                for hit in turkish_hits(value):
                    bad.append(f"line {i} {field}={hit}")
        total += len(bad)
        rep.check(f"language: {label} is English", not bad, len(bad),
                  f"{len(rows)} rows checked" if not bad else first(bad, 4))

    # analysis_out/ icinde web tarafindan yayinlanan metin dosyalari. Bu
    # dosyalar DB'den degil, uretildikleri ANDAKI CSV'den kopyalanir; CSV
    # Ingilizce'ye cevrilince eski ciktilar Turkce kalir.
    bad = []
    for name in PUBLISHED_TEXT_FILES:
        path = os.path.join(out_dir, name)
        if not os.path.exists(path):
            continue
        text = open(path, encoding="utf-8", errors="replace").read()
        hits = turkish_hits(text)
        if hits:
            bad.append(f"{name}: {first(hits, 4)}")
    rep.check("language: published analysis_out files are English",
              not bad, len(bad), first(bad, 4))
    return total


def check_referential(con, rep):
    def scalar(sql):
        return con.execute(sql).fetchone()[0]

    n = scalar("""SELECT COUNT(*) FROM neighbor n
                  LEFT JOIN ro r ON r.candidate_id=n.candidate_id
                  WHERE r.candidate_id IS NULL""")
    rep.check("refs: neighbor.candidate_id exists in ro", n == 0, n)

    n = scalar("""SELECT COUNT(*) FROM ro_leaf x
                  LEFT JOIN leaf l ON l.leaf_id=x.leaf_id
                  WHERE l.leaf_id IS NULL""")
    rep.check("refs: ro_leaf.leaf_id exists in leaf", n == 0, n)

    n = scalar("""SELECT COUNT(*) FROM leaf_sdp s
                  LEFT JOIN cluster_sdp c ON c.cluster=s.cluster
                  WHERE c.cluster IS NULL""")
    rep.check("refs: leaf_sdp.cluster exists in cluster_sdp", n == 0, n)

    n = scalar("""SELECT COUNT(*) FROM leaf_sdp s
                  LEFT JOIN leaf l ON l.leaf_id=s.leaf_id
                  WHERE l.leaf_id IS NULL""")
    rep.check("refs: leaf_sdp.leaf_id exists in leaf", n == 0, n)

    n = scalar("""SELECT COUNT(*) FROM gene_category g
                  LEFT JOIN neighbor n ON n.neighbor_id=g.neighbor_id
                  WHERE n.neighbor_id IS NULL""")
    rep.check("refs: gene_category.neighbor_id exists in neighbor", n == 0, n)

    n = scalar("""SELECT COUNT(*) FROM ro r
                  LEFT JOIN replicon p ON p.nucleotide_id=r.nucleotide_id
                  WHERE p.nucleotide_id IS NULL""")
    rep.check("refs: ro.nucleotide_id exists in replicon", n == 0, n)

    n = scalar("""SELECT COUNT(*) FROM ro_subfamily x
                  LEFT JOIN subfamily s ON s.subfamily_id=x.subfamily_id
                  WHERE s.subfamily_id IS NULL""")
    rep.check("refs: ro_subfamily.subfamily_id exists in subfamily", n == 0, n)

    if con.execute("SELECT name FROM sqlite_master "
                   "WHERE name='replicon_source'").fetchone():
        n = scalar("""SELECT COUNT(*) FROM replicon_source s
                      LEFT JOIN replicon p ON p.nucleotide_id=s.nucleotide_id
                      WHERE p.nucleotide_id IS NULL""")
        rep.check("refs: replicon_source.nucleotide_id exists in replicon",
                  n == 0, n)

    # leaf.size, o yapragin ro_leaf satir sayisina esit olmali. Bu alan web
    # sayfasinda "varyant boyutu" olarak gosteriliyor; ayrismasi sessizdir.
    rows = con.execute("""
        SELECT l.leaf_id, l.size,
               (SELECT COUNT(*) FROM ro_leaf r WHERE r.leaf_id=l.leaf_id)
        FROM leaf l""").fetchall()
    bad = [f"{lid}: size={sz} members={cnt}" for lid, sz, cnt in rows if sz != cnt]
    rep.check("refs: leaf.size equals its ro_leaf member count",
              not bad, len(bad),
              f"{len(rows)} leaves checked" if not bad else first(bad, 4))

    total_size = con.execute("SELECT COALESCE(SUM(size),0) FROM leaf").fetchone()[0]
    total_members = con.execute("SELECT COUNT(*) FROM ro_leaf").fetchone()[0]
    rep.check("refs: sum(leaf.size) equals the ro_leaf row count",
              total_size == total_members, abs(total_size - total_members),
              f"sum(leaf.size)={total_size} ro_leaf={total_members}")

    rows = con.execute("""
        SELECT s.subfamily_id, s.size,
               (SELECT COUNT(*) FROM ro_subfamily r
                WHERE r.subfamily_id=s.subfamily_id)
        FROM subfamily s""").fetchall()
    bad = [f"{sid}: size={sz} members={cnt}" for sid, sz, cnt in rows if sz != cnt]
    rep.check("refs: subfamily.size equals its member count",
              not bad, len(bad),
              f"{len(rows)} subfamilies checked" if not bad else first(bad, 4))


def check_arithmetic(con, rep, n_conf):
    """Dogrulanan RO sayisi her tabloda ayni mi?

    Her biri ayri bir script tarafindan dolduruluyor; biri yarida kalirsa
    sayilar ayrisir ve rapordaki yuzdeler sessizce yanlis olur.
    """
    tables = {
        "ro (is_confirmed=1)": "SELECT COUNT(*) FROM ro WHERE is_confirmed=1",
        "ro_search": "SELECT COUNT(*) FROM ro_search",
        "ro_evidence": "SELECT COUNT(*) FROM ro_evidence",
        "ro_domain": "SELECT COUNT(*) FROM ro_domain",
        "operon": "SELECT COUNT(*) FROM operon",
        "ro_leaf": "SELECT COUNT(*) FROM ro_leaf",
        "ro_etc": "SELECT COUNT(*) FROM ro_etc",
        "ro_regulation": "SELECT COUNT(*) FROM ro_regulation",
        "ro_subfamily": "SELECT COUNT(*) FROM ro_subfamily",
    }
    counts = {}
    for label, sql in tables.items():
        try:
            counts[label] = con.execute(sql).fetchone()[0]
        except sqlite3.DatabaseError:
            counts[label] = None
    distinct = {v for v in counts.values() if v is not None}
    ok = len(distinct) <= 1
    detail = f"all tables report {n_conf}" if ok else \
        "disagreement -> " + ", ".join(f"{k}={v}" for k, v in sorted(counts.items()))
    rep.check("arithmetic: confirmed-RO count agrees across tables",
              ok, len(distinct), detail)

    # Her dogrulanan RO her tabloda TEMSIL EDILIYOR mu (sayilar esit olsa bile
    # kumeler farkli olabilir).
    missing = {}
    for table in ("ro_search", "ro_evidence", "ro_domain", "operon", "ro_leaf",
                  "ro_etc", "ro_regulation", "ro_subfamily"):
        try:
            n = con.execute(f"""
                SELECT COUNT(*) FROM ro r WHERE r.is_confirmed=1
                  AND NOT EXISTS (SELECT 1 FROM "{table}" t
                                  WHERE t.candidate_id=r.candidate_id)""").fetchone()[0]
        except sqlite3.DatabaseError:
            continue
        if n:
            missing[table] = n
    rep.check("arithmetic: every confirmed RO is present in every derived table",
              not missing, sum(missing.values()),
              "" if not missing else ", ".join(f"{k} misses {v}"
                                               for k, v in missing.items()))

    # Ters yon: tureyen tabloda dogrulanmamis giris olmamali.
    extra = {}
    for table in ("ro_search", "ro_evidence", "ro_domain", "operon", "ro_leaf"):
        try:
            n = con.execute(f"""
                SELECT COUNT(*) FROM "{table}" t
                  LEFT JOIN ro r ON r.candidate_id=t.candidate_id
                  WHERE r.candidate_id IS NULL OR r.is_confirmed!=1""").fetchone()[0]
        except sqlite3.DatabaseError:
            continue
        if n:
            extra[table] = n
    rep.check("arithmetic: derived tables hold only confirmed ROs",
              not extra, sum(extra.values()),
              "" if not extra else ", ".join(f"{k} has {v}" for k, v in extra.items()))


def check_search_index_fresh(con, rep, eco, chem):
    """ro_search, kuratorlu CSV'lerin SU ANKI halini mi yansitiyor?

    build_search_index.py substrat/aile/reaksiyon alanlarini CSV'lerden
    KOPYALAR. CSV'ler sonradan guncellenirse tablo bayatlar ve web sitesi eski
    degeri gosterir. SQLite tablo basina mtime tutmadigi icin bunu zaman
    damgasiyla yakalamak mumkun degil; tek yol icerigi karsilastirmak.
    """
    chem_by = {r["cluster"]: r for r in chem}
    eco_by = {r["cluster"]: r for r in eco}
    rows = con.execute("""
        SELECT cluster, substrate, family, reaction, COUNT(*)
        FROM ro_search GROUP BY cluster, substrate, family, reaction""").fetchall()
    bad = []
    entries = 0
    for cluster, substrate, family, reaction, n in rows:
        ch = chem_by.get(cluster, {})
        want_sub = (ch.get("substrate_en") or "").strip() \
            or (eco_by.get(cluster, {}).get("substrate") or "").strip()
        want_fam = (ch.get("family") or "").strip()
        want_rxn = (ch.get("reaction_class") or "").strip()
        got = ((substrate or "").strip(), (family or "").strip(),
               (reaction or "").strip())
        if got != (want_sub, want_fam, want_rxn):
            entries += n
            bad.append(f"{cluster} ({n} entries): index has "
                       f"{got} but the CSVs say "
                       f"{(want_sub, want_fam, want_rxn)}")
    rep.check("freshness: ro_search matches the curated CSVs",
              not bad, len(bad),
              f"{len(rows)} cluster groups checked" if not bad else
              f"{entries} indexed entries affected; " + first(bad, 3))


def check_cross_file_substrate(rep, eco, chem):
    """Iki kuratorlu CSV ayni tip hakkinda celisiyor mu?

    Ozellikle: chemistry.csv substrati ADLANDIRMIS ama cluster_ecology.csv
    hala 'unknown' diyorsa, ekoloji hipotez testi (substrate_class) o tipi
    disarida birakir ve web sayfasi iki farkli sey soyler.
    """
    eco_by = {r["cluster"]: r for r in eco}
    bad = []
    for r in chem:
        e = eco_by.get(r["cluster"])
        if not e:
            continue
        named = (r.get("substrate_en") or "").strip().lower() not in ("", "unknown")
        eco_unknown = (e.get("substrate") or "").strip().lower() in ("", "unknown") \
            or (e.get("substrate_class") or "").strip() == "unknown"
        if named and eco_unknown:
            bad.append(f"{r['cluster']}: chemistry says "
                       f"{r['substrate_en'][:34]!r} but ecology says "
                       f"substrate={e.get('substrate')!r} "
                       f"class={e.get('substrate_class')!r}")
    rep.check("consistency: the two CSVs agree on whether a substrate is known",
              not bad, len(bad), first(bad, 3))


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--db", default=os.path.join(ROOT, "roar.sqlite"))
    ap.add_argument("--ecology", default=os.path.join(ROOT, "cluster_ecology.csv"))
    ap.add_argument("--chemistry", default=os.path.join(ROOT, "chemistry.csv"))
    ap.add_argument("--out-dir", default=os.path.join(ROOT, "analysis_out"))
    ap.add_argument("--refs", default=REFS_FASTA,
                    help="curated reference FASTA the whole pipeline rests on")
    ap.add_argument("--warn-only", action="store_true",
                    help="report everything but always exit 0")
    args = ap.parse_args()

    for path in (args.db, args.ecology, args.chemistry):
        if not os.path.exists(path):
            print(f"CANNOT RUN: {path} not found", file=sys.stderr)
            return 2

    rep = Report()
    print(f"validate_curation.py -- db={os.path.basename(args.db)} "
          f"ecology={os.path.basename(args.ecology)} "
          f"chemistry={os.path.basename(args.chemistry)}")
    print()

    check_vocabulary_drift(rep)
    _, eco = read_csv_strict(args.ecology, rep, "cluster_ecology.csv")
    _, chem = read_csv_strict(args.chemistry, rep, "chemistry.csv")

    con = sqlite3.connect(f"file:{args.db}?mode=ro", uri=True)
    con.text_factory = str
    try:
        check_reference_set(rep, args.refs,
                            os.path.join(args.out_dir, "reference_pairs.csv"),
                            eco, chem, con)
        n_conf = check_confirmed_ro(con, rep)
        check_cluster_coverage(con, rep, eco, chem)
        check_ecology_csv(rep, eco)
        check_chemistry_csv(rep, chem)
        check_cross_file_substrate(rep, eco, chem)
        check_turkish(con, rep, args.ecology, args.chemistry, eco, chem,
                      args.out_dir)
        check_referential(con, rep)
        check_arithmetic(con, rep, n_conf)
        check_search_index_fresh(con, rep, eco, chem)
    finally:
        con.close()

    code = rep.emit()
    if code and not args.warn_only:
        print("\nVALIDATION FAILED -- do not publish this build.", file=sys.stderr)
    return 0 if args.warn_only else code


if __name__ == "__main__":
    sys.exit(main())
