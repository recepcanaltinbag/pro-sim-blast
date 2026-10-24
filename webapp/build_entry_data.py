"""Giris sayfalarinin verisini parcali JSON olarak uret -- statik sitenin agirligi icin.

NEDEN BU SCRIPT VAR
  Statik ihracat her dogrulanmis giris icin bir HTML sayfasi yaziyordu:
  11.423 sayfa x 22,2 KB = 265 MB, ve giris sayisi 19.109'a cikinca ayni
  tasarimla ~405 MB olacakti (site toplami 474 MB, 19 bin dosya). Olculdu:
  sayfanin %66'si komsuluk SVG'si (6,9 KB) ve komsu genler tablosu (7,5 KB) --
  yani AYNI VERININ iki render'i, HTML'e gomulu oldugu icin 19 bin kez
  tekrarlanmis halde. Girisin kendi bilgisi ise 2,2 KB.

  Burada veri BIR KEZ yayimlanir, sayfa tarayicida cizilir. Silinen sey veri
  degil, verinin tekrarlanan HTML kopyasi.

IKI SIKISTIRMA
  1. Komsu satirlari ANAHTARSIZ DIZI + metin sozlugu olarak yazilir.
     474.724 komsu satirinda yalnizca 13.437 ayri `product` metni var; naif
     JSON'da satir 311 B, boyle 49 B.
  2. Operon-regulator SVG'si SUNUCUDA hesaplanir ve (layout, architecture)
     basina BIR KEZ saklanir: 19.109 giris, 4.657 ayri kombinasyon. Boylece
     atlas.operon_regulator_svg (157 satir) JS'e port EDILMEZ -- port edilmeyen
     kod ayrisamaz.

  Komsuluk SVG'si (app.neighborhood_svg, ~50 satir) tarayicida cizilir, cunku
  her girisin koordinatlari farkli. O port verify_entry_render.py ile
  TUM girislerde birebir dizgi karsilastirmasiyla sinanir.

PARCALAMA
  Parca anahtari md5(candidate_id)[:2] -- 256 parca, kararli ve deterministik.
  Bir girisi acan tarayici yalnizca o parcayi indirir.
"""
import argparse
import hashlib
import json
import os
import sqlite3
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

SHARD_HEX = 2            # md5'in ilk 2 hex hanesi -> 256 parca

# Operon SVG'sinde girise ozgu metinleri yer tutucuya cevirmek icin kullanilan
# sentinel'ler. bp int olmak zorunda (fonksiyon int(bp) yapiyor).
BP_SENTINEL = -987654321
BP_SENTINEL_STR = str(BP_SENTINEL)
PROD_SENTINEL = "\u0001ROARPROD\u0001"

# Komsu satirinin dizi duzeni. SIRA ONEMLI: entry.js bu sirayi okur.
NB_FIELDS = ["gene_offset", "distance", "strand", "same_strand", "component",
             "evidence", "category", "product", "gene", "protein_id",
             "start", "end", "spans_origin", "in_operon"]
# Sozluge alinan alanlar: tekrar orani yuksek olanlar.
NB_INTERNED = {"component", "evidence", "category", "product", "gene"}

# `sequence` BU LISTEDE YOK: ayri parcaya gider (asagida seqs).
RO_FIELDS = ["candidate_id", "product", "organism", "taxonomy", "description",
             "nucleotide_id", "replicon_length", "cds_count", "is_circular",
             "is_plasmid", "protein_id", "locus_tag", "gene", "start", "end",
             "strand", "ro_cluster", "ro_group", "leaf_id", "leaf_size",
             "hmm_score", "hmm_evalue", "model_coverage", "rieske_intact",
             "catalytic_intact", "assignment_class", "core_identity",
             "domain", "euk_group", "is_confirmed"]

OPERON_FIELDS = ["n_genes", "start", "end", "has_beta", "has_ferredoxin",
                 "has_reductase", "nearby_beta", "nearby_ferredoxin",
                 "nearby_reductase", "layout"]
OPGENE_FIELDS = ["position", "component", "category", "product", "gap_bp"]
# Operon genleri komsularin ALT KUMESI, yani ayni metinler. Ayni sozluge
# alinirlar ve bedelsiz tekrar kazanci saglanir (olculdu: opg 10,9 MB -> ~3 MB).
OPGENE_INTERNED = {"component", "category", "product"}
EV_FIELDS = ["tier", "nearest_ref", "ref_identity", "ref_qcov", "same_as_hmm"]
ETC_FIELDS = ["etc_profile", "in_operon"]
REG_FIELDS = ["architecture", "operon_5p", "intergenic_bp", "upstream_product",
              "upstream_family", "upstream_divergent"]


def shard_of(candidate_id):
    return hashlib.md5(candidate_id.encode()).hexdigest()[:SHARD_HEX]


class Interner:
    """Metin -> indeks. Ayni metin bir kez yazilir."""

    def __init__(self):
        self.table = []
        self.index = {}

    def __call__(self, value):
        if value is None:
            return None
        key = str(value)
        got = self.index.get(key)
        if got is None:
            got = len(self.table)
            self.index[key] = got
            self.table.append(key)
        return got


def pack(row, fields, intern=None, interned=frozenset()):
    """Satiri ANAHTARSIZ diziye cevir.

    Alan ADLARI her giris icin tekrar yazilmasin: 31 anahtarli `ro` blogunda
    bu tek basina giris basina ~394 B, 19.109 girisle 7,2 MB demekti. Sira
    meta.json'daki *_fields listelerinden okunur; entry.js ayni sirayi kullanir.
    """
    if row is None:
        return None
    keys = row.keys()
    out = []
    for f in fields:
        v = row[f] if f in keys else None
        out.append(intern(v) if (intern and f in interned) else v)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--db", default=os.path.join(os.path.dirname(
        os.path.dirname(os.path.abspath(__file__))), "roar.sqlite"))
    ap.add_argument("--out-dir", required=True,
                    help="statik site koku; data/ altina yazilir")
    ap.add_argument("--limit", type=int, help="gelistirme icin ilk N giris")
    a = ap.parse_args()

    import app as A                                   # noqa: PLC0415
    import atlas as AT                                # noqa: PLC0415

    con = sqlite3.connect(f"file:{a.db}?mode=ro", uri=True)
    con.row_factory = sqlite3.Row

    q = "SELECT candidate_id FROM ro WHERE is_confirmed=1 ORDER BY candidate_id"
    if a.limit:
        q += f" LIMIT {a.limit}"
    ids = [r[0] for r in con.execute(q)]
    print(f"[giris] {len(ids)}")

    shards = {}                 # shard -> {"e": {...}, "s": [sozluk]}
    interners = {}
    opsvg = {}                  # svg_key -> svg  (parcalanmadan once)
    seqs = {}                   # shard -> {candidate_id: dizi}
    opsvg_of = {}               # candidate_id -> svg_key

    for n, cid in enumerate(ids, 1):
        d = A.ro_detail(con, cid)
        if not d:
            continue
        ro, nbs = d["ro"], d["neighbors"]
        operon, opgenes = d["operon"], d["operon_genes"]
        operon_ids = d["operon_ids"]
        x = A.ro_extras(con, cid)
        reg = x.get("ro_regulation")

        sh = shard_of(cid)
        bucket = shards.setdefault(sh, {"e": {}})
        intern = interners.setdefault(sh, Interner())

        rows = []
        for nb in nbs:
            rec = []
            for f in NB_FIELDS:
                if f == "in_operon":
                    rec.append(1 if nb["neighbor_id"] in operon_ids else 0)
                else:
                    v = nb[f] if f in nb.keys() else None
                    rec.append(intern(v) if f in NB_INTERNED else v)
            rows.append(rec)

        # Operon-regulator SVG: SEKIL paylasilir, METIN istemcide yerlesir.
        #
        # Ilk denemede anahtara butun regulasyon alanlari konmustu ve dedup
        # calismadi: 5.000 giris icin 4.987 ayri SVG. Sebep, intergenic_bp ve
        # upstream_product'in SVG icinde METIN olarak gecmesi, yani anahtari
        # neredeyse benzersiz yapmasi.
        #
        # atlas.operon_regulator_svg okundu: bp GEOMETRIYI yalnizca var/yok
        # olarak etkiliyor ((show_reg or bp) kosulu), degeri "N bp" metni;
        # product ise sadece <title> icinde. Dolayisiyla SEKIL yalnizca
        # (layout, architecture, upstream_family, bool(bp)) ile belirlenir.
        #
        # Fonksiyon SENTINEL degerlerle cagrilir ve ciktidaki sentinel'ler
        # yer tutucuya cevrilir. Boylece atlas.py'ye DOKUNULMAZ -- geometrinin
        # tek kaynagi o kalir, JS yalnizca iki dizgi degistirir.
        key = None
        if operon is not None:
            arch = reg["architecture"] if reg and "architecture" in reg.keys() else None
            fam = reg["upstream_family"] if reg and "upstream_family" in reg.keys() else None
            bp = reg["intergenic_bp"] if reg and "intergenic_bp" in reg.keys() else None
            prod = reg["upstream_product"] if reg and "upstream_product" in reg.keys() else None
            shape = json.dumps([operon["layout"], arch, fam, bp is not None], sort_keys=True)
            key = hashlib.md5(shape.encode()).hexdigest()[:12]
            if key not in opsvg:
                svg = AT.operon_regulator_svg(operon["layout"], {
                    "architecture": arch, "upstream_family": fam,
                    "intergenic_bp": BP_SENTINEL if bp is not None else None,
                    "upstream_product": PROD_SENTINEL})
                svg = svg.replace(f"{BP_SENTINEL} bp", "__BP__")
                svg = svg.replace(PROD_SENTINEL, "__PROD__")
                assert BP_SENTINEL_STR not in svg, "bp sentinel ciktida kaldi"
                opsvg[key] = svg
            # Girise ozgu metinler ayrica tasinir.
            bucket.setdefault("t", {})[cid] = [
                (str(int(bp)) + " bp") if bp is not None else "promoter",
                str(prod if prod is not None else "regulator")]
            opsvg_of[cid] = key

        # Dizi AYRI parcada: giris basina ~401 B, 19.109 girisle 7,3 MB. Sayfanin
        # ust kismi onu beklemeden cizilsin diye ayri dosyadan, istendiginde gelir.
        seq = ro["sequence"]
        seqs.setdefault(sh, {})[cid] = seq
        bucket["e"][cid] = {
            "ro": pack(ro, RO_FIELDS),
            "ev": pack(x.get("ro_evidence"), EV_FIELDS),
            "etc": pack(x.get("ro_etc"), ETC_FIELDS),
            "reg": pack(reg, REG_FIELDS),
            "op": pack(operon, OPERON_FIELDS),
            "opg": [pack(g, OPGENE_FIELDS, intern, OPGENE_INTERNED) for g in opgenes],
            "nb": rows,
            "opsvg": key,
        }
        if n % 2000 == 0:
            print(f"   {n}/{len(ids)}")

    data_dir = os.path.join(a.out_dir, "data")
    os.makedirs(os.path.join(data_dir, "entry"), exist_ok=True)
    os.makedirs(os.path.join(data_dir, "opsvg"), exist_ok=True)

    total = 0
    for sh, bucket in sorted(shards.items()):
        bucket["s"] = interners[sh].table
        p = os.path.join(data_dir, "entry", f"{sh}.json")
        with open(p, "w", encoding="utf-8") as fh:
            json.dump(bucket, fh, separators=(",", ":"), ensure_ascii=False)
        total += os.path.getsize(p)
    print(f"[yazildi] data/entry/  {len(shards)} parca, {total/1048576:.1f} MB")

    # Operon SVG'leri de parcalanir; anahtarin ilk 2 hanesi parca.
    buckets = {}
    for key, svg in opsvg.items():
        buckets.setdefault(key[:SHARD_HEX], {})[key] = svg
    tot2 = 0
    for sh, b in sorted(buckets.items()):
        p = os.path.join(data_dir, "opsvg", f"{sh}.json")
        with open(p, "w", encoding="utf-8") as fh:
            json.dump(b, fh, separators=(",", ":"), ensure_ascii=False)
        tot2 += os.path.getsize(p)
    print(f"[yazildi] data/opsvg/  {len(opsvg)} ayri SVG "
          f"({len(opsvg_of)} giris icin), {len(buckets)} parca, {tot2/1048576:.1f} MB")

    os.makedirs(os.path.join(data_dir, "seq"), exist_ok=True)
    tot3 = 0
    for sh, b in sorted(seqs.items()):
        p = os.path.join(data_dir, "seq", f"{sh}.json")
        with open(p, "w", encoding="utf-8") as fh:
            json.dump(b, fh, separators=(",", ":"), ensure_ascii=False)
        tot3 += os.path.getsize(p)
    print(f"[yazildi] data/seq/    {len(seqs)} parca, {tot3/1048576:.1f} MB")

    meta = {
        "shard_hex": SHARD_HEX,
        "nb_fields": NB_FIELDS,
        "nb_interned": sorted(NB_INTERNED),
        "ro_fields": RO_FIELDS,
        "ev_fields": EV_FIELDS,
        "etc_fields": ETC_FIELDS,
        "reg_fields": REG_FIELDS,
        "operon_fields": OPERON_FIELDS,
        "opgene_fields": OPGENE_FIELDS,
        "opgene_interned": sorted(OPGENE_INTERNED),
        "component_colors": A.COMPONENT_COLORS,
        "neighbor_window": A.NEIGHBOR_WINDOW,
        "tier_colors": AT.TIER_COLORS,
        "tier_meaning": AT.TIER_MEANING,
        "cluster_gene": {c: A.cluster_parts(c).get("gene")
                         for c in {r[0] for r in con.execute(
                             "SELECT DISTINCT ro_cluster FROM ro WHERE is_confirmed=1") if r[0]}},
        "cluster_substrate": {k: v.get("substrate") for k, v in A.ECOLOGY.items()},
    }
    p = os.path.join(data_dir, "meta.json")
    with open(p, "w", encoding="utf-8") as fh:
        json.dump(meta, fh, separators=(",", ":"), ensure_ascii=False)
    print(f"[yazildi] data/meta.json  {os.path.getsize(p)/1024:.0f} KB")
    print(f"[toplam]  {(total+tot2+tot3+os.path.getsize(p))/1048576:.1f} MB")
    con.close()


if __name__ == "__main__":
    main()
