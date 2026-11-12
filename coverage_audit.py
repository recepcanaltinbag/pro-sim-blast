"""Kapsama acigini OLCER: dis evrenin hangi taksonlari bizde yok?

SORU. "Veri tabanindaki veriler eksik geliyor" bir payda sorusudur, ve
olculebilir: bu DB'nin dogrulanmis uyeleri, ayni domain'i tasiyan BILINEN
protein evreninin ne kadarini kapsiyor? Dis evren olarak Pfam PF00848
(halka hidroksilleyen oksijenazin alfa alt birimi, katalitik domain)
alinir; InterPro onu UniProtKB uzerinde sayar.

NEDEN DIZI DEGIL TAKSON. Iki tarafin DIZI sayisi ayni birim DEGIL:
InterPro'nun birimi bir UniProtKB kaydi (ayni proteini tasiyan yuz sus
yuz kayit uretir), bizimki genomdan turetilmis bir giris. Dolayisiyla
"%71 dizi kapsamasi" rakami bir kapsama olcusu degil, iki ayri sayma
kuralinin orani. TAKSON ekseninde bu uyusmazlik kayboluyor: bir taxid
her iki tarafta ayni seydir. Bu yuzden karar ekseni takson.

KAPSAM AYRIMI -- BU BETIGIN EN ONEMLI ADIMI. PF00848 bakteriyel degildir:
ayni domain bitki/fungus kolin monooksijenazinda da var. Bu DB ise
bakteriyel genomlardan kurulu. Kapsam ayrilmazsa "eksik takson" sayisi,
bu DB'nin hic hedeflemedigi okaryotlarla sisirilir -- yani olculen sey
kapsama acigi degil, kapsam tanimi olur. O yuzden dis evren superkingdom
rank'indan siniflandirilir ve HEM tam acik HEM kapsam-ici acik ayri ayri
yazilir; kapsam, bizim tarafimizda ne OLDUGU sayilarak beyan edilir.

CINS KURALI. Iki tarafin taxid'leri `taxonomy_ids.genus_from_ranks` ile
cinse cevrilir -- bu projenin cins kararinin tek sahibi o fonksiyondur.
Dis evrenin taxid'leri bu DB'nin hicbir tablosunda olmadigi icin burada
ikinci bir cozucu yazmak cazipti; yazilmadi, cunku ayni olcumun ikinci
uygulamasi birincisiyle sessizce ayrisir.

BU SAYI NE DEGIL
  * Bir cinsin PF00848'de OLMAMASI, o cinste enzimin olmadigini
    GOSTERMEZ. UniProtKB'de temsil edilmemek bir kanit yoklugudur.
  * Bizde olup PF00848'de olmayan bir cins HATA DEGIL: bizim isabetlerimiz
    projenin kendi modelinden gelir, Pfam'in kendi esiginden degil.
  * Eksik cins listesi bir INDIRME LISTESI onerisidir, bir hata listesi
    degil; sirasi "kac protein kazanilir"a gore.

AG. `--fetch` olmadan betik aga CIKMAZ; onbellek yetkilidir. InterPro'nun
sayisi surukleniyor (2026-11 icinde 22.620 -> 22.553), bu yuzden onbellek
tarihiyle ve kendi sayisiyla saklanir ve rapor ikisini de yazar.

Cikti: analysis_out/coverage_audit.json  (rapor)
       analysis_out/pf00848_members.json.gz  (onbellek, --fetch ile)
"""

import argparse
import gzip
import json
import os
import sqlite3
import sys
import time
import urllib.error
import urllib.request
from collections import Counter, defaultdict

import ro_taxonomy
import taxonomy_ids

PFAM_ENTRY = "PF00848"
PFAM_LABEL = ("ring-hydroxylating oxygenase alpha subunit, catalytic domain "
              "(Rieske non-haem iron oxygenase family)")
#  SAYFA BOYU 100, 200 DEGIL. Olculdu: 200'luk sayfalarla InterPro 2.600
#  ve 7.600 kayitta ust uste 408 (Request Timeout) verdi -- sunucu sorguyu
#  kendi zaman siniri icinde bitiremiyor. Kucuk sayfa daha fazla istek
#  demek, ama bitmeyen bir cekimden iyidir.
API = ("https://www.ebi.ac.uk/interpro/api/protein/UniProt/entry/pfam/"
       "{acc}/?page_size=100")
CACHE = "analysis_out/pf00848_members.json.gz"
#  Yarida kalan cekim buraya yazilir; depoda degil, cunku tamamlanmamis
#  bir liste bir artefakt degil.
CHECKPOINT = "/media/lin-bio/back2/roar_taxonomy/pf00848_fetch.partial.json"
OUT = "analysis_out/coverage_audit.json"

#  Rapora yazilan ornek sayisi. Liste bir indirme onerisi; her cins icin
#  birkac ornek accession teshisi mumkun kilar, tam dokum dosyayi sisirir.
EXAMPLES = 3

#  YAYIMLANAN ETIKET, DAHILI ANAHTAR DEGIL. `ro_taxonomy.UNKNOWN_KEY`
#  bir gruplama anahtari ve Turkce ("taxid yok"); yayimlanan metin
#  Ingilizce olmak zorunda, ve anahtari degistirmek butun projede
#  gruplamayi degistirirdi. O yuzden anahtar oldugu gibi kalir ve
#  yalnizca RAPORA yazilirken etiketlenir. Site denetimi bunu yakaladi:
#  dahili bir anahtar yayimlanan bir dosyaya sizdiginda, okuyan kisi
#  onu bir veri degeri sanar.
PUBLISHED_LABEL = {ro_taxonomy.UNKNOWN_KEY: "no genus in the NCBI lineage"}


def label(key):
    return PUBLISHED_LABEL.get(key, key)


def fetch_members(entry, checkpoint, pause=1.0):
    """InterPro'nun PF00848 protein listesini accession+taxid olarak ceker.

    TLS DOGRULAMASI ACIK. Bu depodaki eski `download_pfam.py`
    `ssl._create_unverified_context()` kullaniyor; buraya kopyalanmadi.
    Sertifika dogrulamasini kapatmak, indirilen listenin KAYNAGINI
    dogrulanamaz kilar -- ve bu liste bundan sonra hangi genomlarin
    indirilecegine karar verecek.

    DEVAM EDEBILIR. Olculdu: ilk deneme 7.600/22.553 kayitta InterPro'nun
    408'ine takildi ve cekilen her sey COPE GITTI. 113 sayfalik bir cekim
    tek gecici hataya bagli olamaz; o yuzden her sayfadan sonra kayitlar
    VE siradaki cursor bir ara dosyaya yazilir, ve betik yeniden
    kosturuldugunda kaldigi yerden devam eder. Ara dosya bir onbellek
    DEGIL: tamamlanmamis bir liste kapsama iddiasina girmemeli, bu yuzden
    onbellek yalnizca cekim bittiginde yazilir.
    """
    state = {"url": API.format(acc=entry), "records": [], "total": None,
             "pages": 0}
    if os.path.exists(checkpoint):
        with open(checkpoint, encoding="utf-8") as handle:
            state = json.load(handle)
        sys.stderr.write(f"[devam] {len(state['records'])} kayittan "
                         f"{state['pages']}. sayfadan sonra\n")

    while state["url"]:
        payload = None
        for attempt in range(8):
            try:
                req = urllib.request.Request(
                    state["url"], headers={"Accept": "application/json"})
                with urllib.request.urlopen(req, timeout=180) as res:
                    payload = json.loads(res.read().decode())
                break
            except urllib.error.HTTPError as err:
                #  408 InterPro'nun kendi onerdigi bekleme; digerleri de
                #  gecici olabilir, ama sessizce atlanmaz.
                if attempt == 7:
                    raise
                #  Retry-After varsa sunucunun kendi suresi kullanilir.
                wait = err.headers.get("Retry-After") if err.headers else None
                time.sleep(int(wait) if wait and wait.isdigit()
                           else (61 if err.code == 408 else 10 * (attempt + 1)))
            except (urllib.error.URLError, TimeoutError):
                if attempt == 7:
                    raise
                time.sleep(10 * (attempt + 1))
        if state["total"] is None:
            state["total"] = payload["count"]
        for item in payload["results"]:
            meta = item["metadata"]
            org = meta.get("source_organism") or {}
            state["records"].append([meta.get("accession"), org.get("taxId")])
        state["pages"] += 1
        state["url"] = payload.get("next")
        with open(checkpoint, "w", encoding="utf-8") as handle:
            json.dump(state, handle)
        sys.stderr.write(f"\r[fetch] {len(state['records'])}/{state['total']}")
        sys.stderr.flush()
        if state["url"]:
            time.sleep(pause)
    sys.stderr.write("\n")

    #  Cursor sayfalamasi ayni kaydi iki kere verebilir (kume cekim
    #  sirasinda degisirse) ve devam dosyasi bir sayfayi tekrarlayabilir.
    #  Accession'a gore tekillestirilir; KAC TEKRAR oldugu de yazilir,
    #  cunku sessizce dusurulen bir tekrar, kayip bir kayitla ayni gorunur.
    seen, unique = set(), []
    for accession, taxid in state["records"]:
        if accession in seen:
            continue
        seen.add(accession)
        unique.append([accession, taxid])

    #  InterPro'nun BEYAN ETTIGI sayi ile CEKILEN sayi ayri ayri saklanir:
    #  esit olmadiklari an bu bir sayfalama kaybidir ve raporda gorunur.
    return {"entry": entry,
            "fetched_at": time.strftime("%Y-%m-%d"),
            "interpro_declared_count": state["total"],
            "records_fetched": len(unique),
            "duplicate_records_dropped": len(state["records"]) - len(unique),
            "pages": state["pages"],
            "records": unique}


def read_cache(path):
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        return json.load(handle)


def write_cache(path, data):
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        json.dump(data, handle)


def classify(taxids, nodes, names, merged, deleted):
    """taxid kumesini cinse ve superkingdom'a cevirir.

    Dondurur: taxid -> (genus_key, genus_name, superkingdom). `genus_key`
    `ro_taxonomy` ile AYNI anahtar bicimindedir (vekil rank'li anahtar
    gercek bir cinsle ayni kovaya dusmesin diye).
    """
    out = {}
    for taxid in taxids:
        current, status = taxonomy_ids.resolve(taxid, nodes, merged, deleted)
        if current is None:
            out[taxid] = (f"unresolved:{status}", None, None)
            continue
        ranks = taxonomy_ids.nearest_ranks(current, nodes, names)
        state, gtaxid, gname, proxy_rank, proxy_name = \
            taxonomy_ids.genus_from_ranks(ranks)
        if state == "ranked":
            key = gname
        elif proxy_rank:
            key = f"rank:{proxy_rank}:{proxy_name}"
        else:
            key = ro_taxonomy.UNKNOWN_KEY
        kingdom = ranks.get("superkingdom") or ranks.get("domain")
        out[taxid] = (key, gname, kingdom[1] if kingdom else None)
    return out


def our_taxids(con):
    """DOGRULANMIS uyesi olan GENOM kayitlarinin guncel taxid'leri.

    Payda `ro.is_confirmed=1`: aday satirlar bu DB'nin bir iddiasi degil,
    taranmis bir kumedir, ve kapsama iddiasi yalnizca dogrulanmis uyeler
    uzerinden kurulabilir.

    VE `record_kind='genome'`. Olculdu: bu DB'nin okaryot tarafi 958
    kayittan olusuyor ve 958'inin TAMAMI `XM_`/`NM_` bir mRNA transkripti
    (930 + 28); tek bir okaryot genomu yok. Transkriptler payda icinde
    birakilirsa "Eukaryota" bu DB'nin KAPSAMI sayilir, ve eksik takson
    listesinin basina bir fungus (Clonostachys, 67 protein) cikar -- yani
    olculen sey kapsama acigi olmaktan cikip bir kapsam hatasi olur.
    Kayit turu `replicon_kind` tablosundan okunur (ROADMAP 135.7); kaydin
    kendi LOCUS beyanina dayanir, accession onegine degil.
    """
    rows = con.execute("""
        SELECT DISTINCT t.taxid_current
          FROM replicon_taxon t
          JOIN ro r ON r.nucleotide_id = t.nucleotide_id
          JOIN replicon_kind k ON k.nucleotide_id = t.nucleotide_id
         WHERE r.is_confirmed = 1 AND t.taxid_current IS NOT NULL
           AND k.record_kind = 'genome'""").fetchall()
    excluded = con.execute("""
        SELECT COUNT(DISTINCT t.taxid_current)
          FROM replicon_taxon t
          JOIN ro r ON r.nucleotide_id = t.nucleotide_id
          JOIN replicon_kind k ON k.nucleotide_id = t.nucleotide_id
         WHERE r.is_confirmed = 1 AND t.taxid_current IS NOT NULL
           AND k.record_kind <> 'genome'""").fetchone()[0]
    return {int(r[0]) for r in rows}, excluded


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--db", default="roar.sqlite")
    parser.add_argument("--dump", default="/media/lin-bio/back2/roar_taxonomy")
    parser.add_argument("--cache", default=CACHE)
    parser.add_argument("--out", default=OUT)
    parser.add_argument("--checkpoint", default=CHECKPOINT,
                        help="yarida kalan cekimin devam dosyasi")
    parser.add_argument("--fetch", action="store_true",
                        help="InterPro'ya cikip onbellegi yeniler")
    args = parser.parse_args()

    if args.fetch:
        write_cache(args.cache, fetch_members(PFAM_ENTRY, args.checkpoint))
        os.remove(args.checkpoint)
    if not os.path.exists(args.cache):
        sys.exit(f"onbellek yok: {args.cache} (once --fetch ile kos)")
    cache = read_cache(args.cache)
    print(f"[dis]   {cache['records_fetched']} kayit, "
          f"{cache['fetched_at']} tarihli onbellek")

    nodes, names, merged, deleted = taxonomy_ids.read_taxdump(args.dump)
    con = sqlite3.connect(args.db)

    mine, excluded_taxids = our_taxids(con)
    outside = defaultdict(list)                   # taxid -> [accession, ...]
    no_taxid = 0
    for accession, taxid in cache["records"]:
        if taxid is None:
            no_taxid += 1
            continue
        outside[int(taxid)].append(accession)

    tax = classify(set(outside) | mine, nodes, names, merged, deleted)

    #  KAPSAM, bizim tarafimizda ne oldugunu SAYARAK beyan edilir.
    my_kingdoms = Counter(tax[t][2] for t in mine)
    in_scope = {k for k, n in my_kingdoms.items() if k and n >= 1}

    my_genera = {tax[t][0] for t in mine}
    ext_genus = defaultdict(lambda: {"proteins": 0, "taxids": set(),
                                     "examples": [], "kingdoms": set(),
                                     "genus_taxid": None})
    for taxid, accessions in outside.items():
        key, gname, kingdom = tax[taxid]
        rec = ext_genus[key]
        rec["proteins"] += len(accessions)
        rec["taxids"].add(taxid)
        rec["kingdoms"].add(kingdom)
        rec["genus_name"] = gname
        for acc in accessions[:EXAMPLES]:
            if len(rec["examples"]) < EXAMPLES:
                rec["examples"].append(acc)

    def rows_for(keys):
        out = []
        for key in sorted(keys, key=lambda k: -ext_genus[k]["proteins"]):
            rec = ext_genus[key]
            out.append({"genus_key": label(key),
                        "genus_name": rec.get("genus_name"),
                        "superkingdom": sorted(
                            k for k in rec["kingdoms"] if k) or None,
                        "pfam_proteins": rec["proteins"],
                        "distinct_taxids": len(rec["taxids"]),
                        "example_accessions": rec["examples"]})
        return out

    def scoped(keys):
        return [k for k in keys
                if any(x in in_scope for x in ext_genus[k]["kingdoms"])]

    #  CINS OLMAYAN anahtarlar listeden cikarilir. `rank:...` bir vekil
    #  (NCBI o soyagacinda cins rank'li dugum tutmuyor), `taxid yok` ise
    #  cozulemeyen taxid. Ikisi de bir INDIRME onerisi olamaz: "eksik
    #  cins" diye sunulan bir satirin karsiliginda indirilecek bir cins
    #  yok. Sayilari dusurulmez, ayri blokta yazilir.
    def is_genus(key):
        return not (key.startswith("rank:") or key.startswith("unresolved:")
                    or key == "taxid yok")

    unresolved = {k: ext_genus[k]["proteins"]
                  for k in ext_genus if not is_genus(k)}
    missing = [k for k in ext_genus if k not in my_genera and is_genus(k)]
    shared = [k for k in ext_genus if k in my_genera and is_genus(k)]
    missing_in_scope = scoped(missing)

    #  Acigin NERESINDE oldugu: kapsanan cinslerde duran protein sayisi ile
    #  hic kapsanmayan cinslerde duran protein sayisi ayri ayri. Bu ikisi
    #  ayri isler: birincisi derinlik, ikincisi genislik.
    proteins_in_shared = sum(ext_genus[k]["proteins"] for k in shared)
    proteins_missing = sum(ext_genus[k]["proteins"] for k in missing)
    proteins_missing_scope = sum(ext_genus[k]["proteins"]
                                 for k in missing_in_scope)

    report = {
        "question": ("Which taxa carrying this domain in the known protein "
                     "universe are absent from this database?"),
        "external_universe": {
            "entry": PFAM_ENTRY,
            "label": PFAM_LABEL,
            "source": "InterPro API, UniProtKB proteins matching the Pfam entry",
            "cache_date": cache["fetched_at"],
            "interpro_declared_count": cache["interpro_declared_count"],
            "records_fetched": cache["records_fetched"],
            "duplicate_records_dropped": cache.get(
                "duplicate_records_dropped", 0),
            "records_without_taxid": no_taxid,
            #  YON ONEMLI. Tehlikeli yon EKSIK cekmektir: yarim bir dis
            #  liste paydayi kucultur ve kapsama oranini oldugundan YUKSEK
            #  gosterir. FAZLA cekmek bunun tersi -- ve olculdu: InterPro
            #  ilk sayfada 22.553 beyan ediyor, cursor'u sonuna kadar
            #  yurumek 22.620 AYRI accession veriyor; yani `count` bir
            #  onbellek tahmini, yurunen liste ise tukenmis.
            "walk_not_truncated": (cache["records_fetched"]
                                   >= cache["interpro_declared_count"]),
            "declared_minus_fetched": (cache["interpro_declared_count"]
                                       - cache["records_fetched"]),
        },
        "method": {
            "decision_axis": "taxon, not sequence",
            "why_not_sequence": (
                "the two sides do not count the same unit: InterPro counts one "
                "UniProtKB record per protein (a hundred strains carrying the "
                "same protein give a hundred records), while this database "
                "counts one genome-derived entry; a ratio of two different "
                "counting rules is not a coverage measure. A taxid, by "
                "contrast, is the same thing on both sides."),
            "genus_rule": ("the nearest ancestor whose NCBI rank is exactly "
                           "'genus'; decided by taxonomy_ids.genus_from_ranks, "
                           "the single owner of the genus question"),
            "our_denominator": ("replicons carrying at least one confirmed "
                                "member (ro.is_confirmed=1)"),
            "scope": {
                "in_scope_superkingdoms": sorted(in_scope),
                "why": ("PF00848 is not bacterial: the same domain occurs in "
                        "plant and fungal choline monooxygenase, while this "
                        "database is built from bacterial genomes. Without the "
                        "split, the missing-taxon count is inflated by lineages "
                        "this database never targeted, so what is measured "
                        "stops being a coverage gap and becomes a scope "
                        "definition."),
                "transcript_only_superkingdom_excluded": (
                    "Eukaryota enters this database only through 958 XM_/NM_ "
                    "mRNA records and no eukaryotic genome; counting it as "
                    "scope would put a fungus at the top of the missing list"),
                "our_side_by_superkingdom": dict(
                    sorted(((str(k), n) for k, n in my_kingdoms.items()),
                           key=lambda kv: -kv[1])),
            },
        },
        "what_this_is_not": [
            "A genus absent from PF00848 is not a genus without the enzyme: "
            "not being represented in UniProtKB is an absence of evidence.",
            "A genus present here but absent from PF00848 is not an error: our "
            "hits come from the project's own model at its own threshold, not "
            "from the Pfam gathering threshold.",
            "The missing-genus list is a proposed download list, not an error "
            "list; it is ordered by how many proteins would be gained.",
        ],
        "denominators": {
            "our_taxids": len(mine),
            "our_taxids_excluded_as_transcript_only": excluded_taxids,
            "our_genera": len(my_genera),
            "external_taxids": len(outside),
            "external_genera": len([k for k in ext_genus if is_genus(k)]),
            "external_keys_that_are_not_a_genus": len(unresolved),
            "external_proteins_under_those_keys": sum(unresolved.values()),
        },
        "gap": {
            "genera_missing": len(missing),
            "genera_missing_in_scope": len(missing_in_scope),
            "genera_shared": len(shared),
            "genera_ours_only": len(my_genera - set(ext_genus)),
            "pfam_proteins_in_shared_genera": proteins_in_shared,
            "pfam_proteins_in_missing_genera": proteins_missing,
            "pfam_proteins_in_missing_genera_in_scope": proteins_missing_scope,
            "where_the_gap_is": (
                "depth (more proteins inside genera already covered)"
                if proteins_in_shared >= proteins_missing
                else "breadth (proteins in genera not covered at all)"),
        },
        "external_taxa_with_no_genus": {
            label(k): v for k, v in
            sorted(unresolved.items(), key=lambda kv: -kv[1])},
        "missing_genera_in_scope": rows_for(missing_in_scope),
        "missing_genera_out_of_scope": rows_for(
            [k for k in missing if k not in set(missing_in_scope)]),
    }

    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    with open(args.out, "w", encoding="utf-8") as handle:
        json.dump(report, handle, indent=2, ensure_ascii=False)
        handle.write("\n")

    print(f"[bizde] {len(mine)} taxid, {len(my_genera)} cins")
    print(f"[dis]   {len(outside)} taxid, {len(ext_genus)} cins")
    print(f"[acik]  {len(missing)} cins eksik, "
          f"{len(missing_in_scope)} tanesi kapsam ici "
          f"({proteins_missing_scope} protein)")
    print(f"[yon]   {report['gap']['where_the_gap_is']}")
    print(f"[yazildi] {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
