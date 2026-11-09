"""
Komsu genlerin HAREKETLILIK ve DUZENLEYICI ailesi -- urun adindan degil, modelden.

NEDEN. `build_db.py` komsu genlere kategoriyi URUN ADINDAN duzenli ifadeyle
koyar (`gene_category.method = 'regex_v1'`). Iki olculmus zaafi var:

  1. "transposon" etiketi urun adinda transpozaz GECMEYEN genleri de topluyor;
     `control_elements.py` bunu zaten yaziyor: 117 kez "Holliday junction
     resolvase RuvX" transpozon sayiliyor. Rekombinasyon makinesi ile
     hareketli eleman ayni kovaya giriyor.
  2. Duzenleyici ailesi (`regulator_family_v1`) ancak ad aileyi SOYLUYORSA
     doluyor. "transcriptional regulator" diyen bir ad aile vermez; o yuzden
     40.524 duzenleyicinin buyuk kismi `regulator_unclassified` kaliyor ve
     `ro_regulation.upstream_family` cogunlukla bos.

BU SCRIPT ad yerine Pfam modeline bakar. Iki AYRI kutuphane kullanilir
(`family_hmms/mobility.hmm`, `family_hmms/regulator.hmm`) ve ikisi de
`component_hmms/components.hmm` ile KARISTIRILMAZ: bilesen taramasindaki her
model bir operon bileseni etiketine baglidir, oraya baska soru sokmak o
taramanin anlamini degistirir.

NE YAZILIR
  - `neighbor_family_hit`: protein basina ham model hitleri (skor, E, kaplama).
  - `neighbor_family`: UC yeni yontem; `gene_category` ile AYNI sekilde
    (neighbor_id, category, method). AYRI TABLO olmasinin sebebi su:
    `gene_category` erken yazilan ve neredeyse her adimin girdisi sayilan
    paylasilan bir tablodur; oraya satir eklemek yeniden kurma grafiginde
    63 adimi birden bayat gosterir, oysa hicbirinin cikti sayisi degismez.
    Paylasilan kaynagi kirletmek yerine turetilmis tablo kurulur -- tip
    basina alan icin `chemistry.csv`de uygulanan kuralin aynisi.
      method='mobility_hmm_v1'          -> category = IS3-like | Tn3 | IS110 | ...
      method='recombinase_hmm_v1'       -> category = tyrosine_recombinase |
                                           serine_recombinase | holliday_resolvase
      method='regulator_family_hmm_v1'  -> category = LysR | TetR | ... | AraC
    Rekombinaz AYRI yontemdir cunku hareketli eleman DEGILDIR; regex ikisini
    ayni 'transposon' kovasina koyuyor ve ayirmak icin once ayri olculmeleri
    gerekir. Yontemler yan yana durdugu icin UYUSMAZLIK dogrudan
    sorgulanabilir; yayimlanan sey de zaten uyusmazligin kendisidir.

NE YAZILMAZ. `ro_regulation.upstream_family` ELLENMEZ. O sutunu
`analyze_regulation.py` yazar ve kaynagi regex'tir; HMM sonucu oraya sessizce
sizarsa sutunun ne anlattigi belirsizlesir. HMM okumasi kendi yonteminde
durur, karsilastirma turetilmis tabloda yapilir.

KAPSAM UYARISI. `neighbor_protein` yalnizca DOGRULANMIS RO'larin komsularini
tutar; regex ise tum komsulara uygulanmistir. Bu yuzden HMM kapsamı regex
kapsamindan kucuktur ve ciktida iki sayi da ayri verilir.

Cikti: neighbor_family_hit + neighbor_family + ro_upstream_family_hmm tablolari
       analysis_out/family_scan.json
"""

import argparse
import hashlib
import json
import os
import sqlite3
import subprocess
from collections import Counter, defaultdict

from ro_filter import parse_domtblout
from stale_intermediate import stale

HMM_DIR = "family_hmms"
WORK_DIR = "family_work"
LIBRARIES = ("mobility", "recombinase", "regulator")

# Bilesen taramasiyla AYNI esik (build_operons.DEFAULT_EVALUE); iki tarama
# arasinda farkli esik kullanmak karsilastirmayi bozar.
DEFAULT_EVALUE = 1e-5

# Aile atamasi icin modelin en az yarisi kaplanmali. HTH alanlari kisa
# (47-80 aa) oldugu icin tam kaplama beklemek duyarliligi gereksiz dusurur;
# yarim kaplamanin altindaki hit "model gorunmus" saymaz.
MIN_COVERAGE = 0.5

# Model -> aile. Bir proteinde birden fazla model varsa aile EN IYI SKORLU
# modelden gelir; ad ciftleri (HTH + efektor alani) zaten ayni aileyi verir.
MOBILITY_FAMILY = {
    "DDE_Tnp_1": "IS4/IS5-like",
    "DDE_Tnp_1_7": "IS4/IS5-like",
    "rve": "IS3-like",
    "HTH_Tnp_1": "IS3-like",
    "Transposase_mut": "MULE",
    "DDE_Tnp_IS66": "IS66",
    "TnpB_IS66": "IS66",
    "DDE_Tnp_Tn3": "Tn3",
    "Transposase_20": "IS110",
    "Y1_Tnp": "IS200/IS605",
    "Transposase_1": "IS630/Tc1",
    "Zn_ribbon_InsA": "IS1",
}
RECOMBINASE_FAMILY = {
    "Phage_integrase": "tyrosine_recombinase",
    "Resolvase": "serine_recombinase",
    "HTH_7": "serine_recombinase",
    "Recombinase": "serine_recombinase",
    "Zn_ribbon_recom": "serine_recombinase",
    "RuvC": "holliday_resolvase",
    "RuvX": "holliday_resolvase",
}
REGULATOR_FAMILY = {
    "TetR_N": "TetR",
    "GntR": "GntR",
    "HTH_1": "LysR",
    "LysR_substrate": "LysR",
    "MarR": "MarR",
    "AraC_binding": "AraC",
    "HTH_18": "AraC",
    "HTH_IclR": "IclR",
    "IclR_C": "IclR",
}
FAMILY_MAP = {"mobility": MOBILITY_FAMILY,
              "recombinase": RECOMBINASE_FAMILY,
              "regulator": REGULATOR_FAMILY}
METHOD = {"mobility": "mobility_hmm_v1",
          "recombinase": "recombinase_hmm_v1",
          "regulator": "regulator_family_hmm_v1"}


def write_fasta(con, path):
    """Icerik degismediyse dosyaya DOKUNMA -- mtime tetikler, ozet karar verir."""
    tmp = path + ".new"
    digest = hashlib.md5()
    count = 0
    with open(tmp, "w") as handle:
        for key, seq in con.execute(
                "SELECT protein_key, translation FROM neighbor_protein "
                "ORDER BY protein_key"):
            record = f">{key}\n{seq}\n"
            handle.write(record)
            digest.update(record.encode())
            count += 1
    new = digest.hexdigest()
    old = None
    if os.path.exists(path):
        h = hashlib.md5()
        with open(path, "rb") as fh:
            for chunk in iter(lambda: fh.read(1 << 20), b""):
                h.update(chunk)
        old = h.hexdigest()
    if old == new:
        os.remove(tmp)
        print(f"[korundu] {path} degismedi ({count:,} protein)")
    else:
        os.replace(tmp, path)
        print(f"[yazildi] {path} ({count:,} protein)")
    return count


def search(library, fasta, cpu, force):
    hmm = os.path.join(HMM_DIR, f"{library}.hmm")
    if not os.path.exists(hmm):
        raise SystemExit(f"[HATA] {hmm} yok -- once: "
                         f"python3 fetch_family_hmms.py --library {library}")
    domtbl = os.path.join(WORK_DIR, f"{library}.domtbl")
    why = "force" if force else stale(domtbl, fasta, hmm)
    if why:
        print(f"[uretiliyor] {domtbl} ({why})")
        subprocess.run(["hmmsearch", "--domtblout", domtbl, "--noali",
                        "--cpu", str(cpu), "-E", str(DEFAULT_EVALUE), hmm, fasta],
                       check=True, stdout=subprocess.DEVNULL)
    else:
        print(f"[atlandi] {domtbl} guncel")
    return domtbl


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--db", default="roar.sqlite")
    ap.add_argument("--out", default="analysis_out/family_scan.json")
    ap.add_argument("--work-dir", default=WORK_DIR)
    ap.add_argument("--cpu", type=int, default=8)
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()

    os.makedirs(args.work_dir, exist_ok=True)
    con = sqlite3.connect(args.db)
    con.executescript("""
        CREATE TABLE IF NOT EXISTS neighbor_family_hit (
            protein_key TEXT REFERENCES neighbor_protein(protein_key),
            library     TEXT,
            model       TEXT,
            score       REAL,
            evalue      REAL,
            coverage    REAL
        );
        CREATE INDEX IF NOT EXISTS idx_nfh_key ON neighbor_family_hit(protein_key);
        CREATE TABLE IF NOT EXISTS neighbor_family (
            neighbor_id INTEGER REFERENCES neighbor(neighbor_id),
            category    TEXT,
            method      TEXT
        );
        CREATE INDEX IF NOT EXISTS idx_nf_method ON neighbor_family(method);
    """)

    fasta = os.path.join(args.work_dir, "neighbor_proteins.fasta")
    n_proteins = write_fasta(con, fasta)

    # protein_key -> [neighbor_id, ...]  (ayni gen birden cok RO'nun komsusu olabilir)
    key_to_neighbors = defaultdict(list)
    for nid, nuc, start, end, strand in con.execute(
            "SELECT neighbor_id, nucleotide_id, start, end, strand FROM neighbor"):
        key_to_neighbors[f"{nuc}:{start}-{end}:{strand}"].append(nid)

    report = {"n_neighbour_proteins": n_proteins,
              "evalue": DEFAULT_EVALUE, "min_coverage": MIN_COVERAGE,
              "libraries": {}}
    for library in LIBRARIES:
        domtbl = search(library, fasta, args.cpu, args.force)
        hits = parse_domtblout(domtbl)
        mapping = FAMILY_MAP[library]

        rows, assigned, model_counts = [], {}, Counter()
        for key, hs in hits.items():
            best = None
            for h in hs:
                rows.append((key, library, h["query_name"], h["score"],
                             h["e_value"], round(h["model_coverage"], 3)))
                if h["model_coverage"] < MIN_COVERAGE:
                    continue
                model_counts[h["query_name"]] += 1
                if best is None or h["score"] > best["score"]:
                    best = h
            if best is not None:
                assigned[key] = mapping.get(best["query_name"], "unassigned_model")

        con.execute("DELETE FROM neighbor_family_hit WHERE library = ?", (library,))
        con.executemany("INSERT INTO neighbor_family_hit VALUES (?,?,?,?,?,?)", rows)

        method = METHOD[library]
        cat_rows = [(nid, fam, method)
                    for key, fam in assigned.items()
                    for nid in key_to_neighbors.get(key, ())]
        con.execute("DELETE FROM neighbor_family WHERE method = ?", (method,))
        con.executemany(
            "INSERT INTO neighbor_family (neighbor_id, category, method) VALUES (?,?,?)",
            cat_rows)
        con.commit()

        report["libraries"][library] = {
            "method": method,
            "models": json.load(open(os.path.join(HMM_DIR,
                                                  f"{library}.manifest.json")))["models"],
            "proteins_with_any_hit": len(hits),
            "proteins_assigned": len(assigned),
            "neighbour_rows_assigned": len(cat_rows),
            "model_hits_above_coverage": dict(model_counts.most_common()),
            "families": dict(Counter(assigned.values()).most_common()),
        }
        print(f"[{library}] {len(hits):,} proteinde hit, {len(assigned):,} proteine aile, "
              f"{len(cat_rows):,} komsu satiri")

    # ------------------------------------------------ regex ile karsilastirma
    regex_cat, regex_fam = {}, {}
    for nid, cat, method in con.execute(
            "SELECT neighbor_id, category, method FROM gene_category "
            "WHERE method IN ('regex_v1','regulator_family_v1')"):
        (regex_cat if method == "regex_v1" else regex_fam)[nid] = cat

    product = dict(con.execute("SELECT neighbor_id, product FROM neighbor"))

    # 1. Hareketlilik: regex "transposon" der; model ne der?
    #
    # UC KUME, IKI AYRI SORU. Regex etiketi tek bir kovadir; model iki ayri
    # sey gorur:
    #   - transpozaz cekirdegi  -> gercekten hareketli eleman,
    #   - rekombinaz/integraz   -> DNA'yi kesip baglar ama HAREKETLI ELEMAN
    #                              DEGILDIR (XerC, RuvX, cozucu HTH'ler).
    # Ikisi ayri olculmezse "regex fazla sayiyor" denemez: fazlaligin ne
    # kadari yanlis sinif, ne kadari eksik model belli olmaz.
    mob_nids = {nid for nid, in con.execute(
        "SELECT neighbor_id FROM neighbor_family WHERE method = 'mobility_hmm_v1'")}
    rec_nids = {nid for nid, in con.execute(
        "SELECT neighbor_id FROM neighbor_family WHERE method = 'recombinase_hmm_v1'")}
    # TARANAN KUME: neighbor_protein'i olan her komsu satiri. Regex tum
    # komsulara uygulandi, HMM yalnizca bunlara; karsilastirma ancak ayni
    # kume uzerinde anlamli.
    scanned = set()
    for key, in con.execute("SELECT protein_key FROM neighbor_protein"):
        scanned.update(key_to_neighbors.get(key, ()))

    regex_mob = {nid for nid in scanned if regex_cat.get(nid) == "transposon"}
    misclassified = (regex_mob & rec_nids) - mob_nids
    unexplained = regex_mob - mob_nids - rec_nids
    report["mobility_vs_regex"] = {
        "scanned_neighbour_rows": len(scanned),
        "regex_transposon": len(regex_mob),
        "hmm_transposase": len(mob_nids),
        "hmm_recombinase": len(rec_nids),
        "transposase_confirmed": len(regex_mob & mob_nids),
        "recombinase_not_mobile": len(misclassified),
        "no_model_hit": len(unexplained),
        "transposase_missed_by_regex": len(mob_nids - regex_mob),
        "recombinase_outside_regex_transposon": len(rec_nids - regex_mob),
        "recombinase_not_mobile_top_products":
            dict(Counter(product.get(n) or "" for n in misclassified).most_common(15)),
        "no_model_hit_top_products":
            dict(Counter(product.get(n) or "" for n in unexplained).most_common(15)),
        "transposase_missed_by_regex_top_products":
            dict(Counter(product.get(n) or ""
                         for n in mob_nids - regex_mob).most_common(10)),
    }

    # 2. Duzenleyici ailesi: regex aile vermemisse model veriyor mu?
    hmm_fam = {}
    for nid, cat in con.execute(
            "SELECT neighbor_id, category FROM neighbor_family "
            "WHERE method = 'regulator_family_hmm_v1'"):
        hmm_fam[nid] = cat
    regulators = {nid for nid in scanned if regex_cat.get(nid) == "regulator"}
    unclassified = {nid for nid in regulators
                    if regex_fam.get(nid) in (None, "regulator_unclassified")}
    agree = Counter()
    for nid in regulators:
        r, h = regex_fam.get(nid), hmm_fam.get(nid)
        if r in (None, "regulator_unclassified"):
            agree["regex_blank_hmm_named" if h else "both_blank"] += 1
        elif h is None:
            agree["regex_named_hmm_blank"] += 1
        else:
            agree["agree" if r == h else "disagree"] += 1
    conflicts = Counter((regex_fam[nid], hmm_fam[nid]) for nid in regulators
                        if regex_fam.get(nid) not in (None, "regulator_unclassified")
                        and hmm_fam.get(nid) and regex_fam[nid] != hmm_fam[nid])
    report["regulator_vs_regex"] = {
        "regex_regulator_rows": len(regulators),
        "regex_family_blank": len(unclassified),
        "rescued_by_hmm": agree["regex_blank_hmm_named"],
        "agreement": dict(agree),
        "disagreements": {f"{a}->{b}": n for (a, b), n in conflicts.most_common(15)},
        "hmm_only_families":
            dict(Counter(hmm_fam[nid] for nid in unclassified
                         if nid in hmm_fam).most_common()),
        "hmm_outside_regex_regulator":
            len(set(hmm_fam) - regulators),
    }

    # ------------------------------- 3. operon basina yukari akis ailesi
    # NEDEN AYRI TABLO. `ro_regulation.upstream_family` regex'in yazdigi
    # sutundur ve sahibi `analyze_regulation.py`dir. HMM okumasi oraya
    # yazilirsa sutunun ne anlattigi belirsizlesir (iki yontem tek ad).
    # Bu yuzden birlesik okuma TURETILMIS tabloda durur ve hangi satirin
    # hangi yontemden geldigi `source` ile acikca yazilir.
    con.executescript("""
        CREATE TABLE IF NOT EXISTS ro_upstream_family_hmm (
            candidate_id  TEXT PRIMARY KEY,
            family_regex  TEXT,
            family_hmm    TEXT,
            family_best   TEXT,   -- regex varsa o, yoksa HMM
            source        TEXT    -- regex | hmm | none | no_upstream_gene
        );
    """)
    reg_rows, source_counts = [], Counter()
    for cid, nid, fam in con.execute(
            "SELECT candidate_id, upstream_gene_id, upstream_family FROM ro_regulation"):
        named = fam if fam and fam != "regulator_unclassified" else None
        hmm = hmm_fam.get(nid)
        if named:
            best, source = named, "regex"
        elif hmm:
            best, source = hmm, "hmm"
        elif nid is None:
            best, source = None, "no_upstream_gene"
        else:
            best, source = None, "none"
        source_counts[source] += 1
        reg_rows.append((cid, named, hmm, best, source))
    con.execute("DELETE FROM ro_upstream_family_hmm")
    con.executemany("INSERT INTO ro_upstream_family_hmm VALUES (?,?,?,?,?)", reg_rows)
    con.commit()
    report["ro_upstream_family"] = {
        "n_ro": len(reg_rows),
        "by_source": dict(source_counts),
        "named_before": source_counts["regex"],
        "named_after": source_counts["regex"] + source_counts["hmm"],
        "families_added_by_hmm":
            dict(Counter(h for _c, r, h, _b, s in reg_rows
                         if s == "hmm").most_common()),
    }
    u = report["ro_upstream_family"]
    print(f"[yukari akis ailesi] {u['named_before']:,} -> {u['named_after']:,} "
          f"adlandirilmis RO ({len(reg_rows):,} icinde)")

    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    with open(args.out, "w") as fh:
        json.dump(report, fh, indent=1)
    print(f"\n[yazildi] {args.out}")
    m = report["mobility_vs_regex"]
    print(f"[hareketlilik] regex 'transposon' {m['regex_transposon']:,} satir -> "
          f"transpozaz dogrulandi {m['transposase_confirmed']:,}, "
          f"rekombinaz (hareketli degil) {m['recombinase_not_mobile']:,}, "
          f"hicbir model tutmadi {m['no_model_hit']:,}")
    print(f"[hareketlilik] regex'in kacirdigi transpozaz: "
          f"{m['transposase_missed_by_regex']:,}")
    r = report["regulator_vs_regex"]
    print(f"[duzenleyici] regex duzenleyici {r['regex_regulator_rows']:,}, "
          f"ailesi bos {r['regex_family_blank']:,}, modelin adlandirdigi "
          f"{r['rescued_by_hmm']:,}")
    con.close()


if __name__ == "__main__":
    main()
