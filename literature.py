"""Her tip icin LITERATUR tanimlayicilarini (DOI, PubMed) toplar.

SORU. Tip sayfalari enzimin substratini, reaksiyonunu ve kuratorluk
gerekcesini yaziyor, ama okuyucu "bunu nereden biliyoruz" diye sorup
makaleye gitmek istediginde elinde bir baglanti yok. Gerekce `source`
alaninda duz yazi olarak duruyor ve icindeki DOI/PMID cumlenin ortasinda
kaliyor.

UC KAYNAK, UC AYRI IDDIA. Bu betigin en onemli karari, uc kaynagi ayni
kefeye KOYMAMAK. Her kayit hangi yolla geldigini tasir, cunku uclu ayni
seyi iddia etmiyor:

  curated       Bizim kendi kuratorluk metnimizde YAZILI olan tanimlayici.
                En guclusu: bu projenin atama gerekcesi olarak gosterdigi
                kaynak. Kuratorun cumlesi de saklanir.
  pdb_primary   Tipin PDB kaydinin BIRINCIL atifi. Yapidan gelen atamalarda
                (source_kind'i structure olanlar) bu, atamanin dayandigi
                makaledir. Yapisi olup atamasi baska yerden gelen tiplerde
                ise "bu yapiyi cozen makale"dir, baskasi degil.
  uniprot       Referans proteinin UniProt kaydinin atif listesi. Bu
                "bu enzimi karakterize eden makale" DEMEK DEGILDIR: giris
                icin gosterilen her atif listeye girer. UniProt bir atifi
                FUNCTION ya da CATALYTIC ACTIVITY aciklamasina KANIT olarak
                bagladiysa, o kayit `function_evidence` ile isaretlenir --
                ancak o isaret bir islev iddiasinin kaynagidir.

Uc kaynagi birlestirip "tipin makalesi" diye tek liste yayinlamak kolaydi;
yayinlanmadi, cunku o liste olcumun kendisinden daha fazlasini iddia
ederdi.

ALAN ADLARI KAYNAGINDAN OKUNUR. RCSB'nin atif alanlari BUYUK harfli:
`pdbx_database_id_DOI`, `pdbx_database_id_PubMed`. Ilk denemede kucuk
harfle sorulduklari icin uc yapinin ucunde de `None` dondu ve sonuc
"bu yapilarin DOI'si yok" gibi okundu. Aracin urettigi yokluk, verinin
yoklugundan ayirt edilemiyordu; bu yuzden alan adlari burada sabit ve
beklenen tipleriyle birlikte yaziliyor.

DOI COZULUR. Yayinlanan her DOI doi.org'a HEAD ile sorulur ve cozulup
cozulmedigi kaydedilir. Cozulmeyen bir DOI yanlis bir atiftir, ve yanlis
atif eksik atiftan kotudur: okuyucu onu dogru sanip pesine gider.
Cozulmeyenler ciktida KALIR ama `resolved: false` ile, ve kapi onlari
sayar.

BU SAYI NE DEGIL
  * Bir tipin DOI'si olmamasi, literaturu olmadigini GOSTERMEZ; bizim
    kuratorluk metnimizde tanimlayici yazmamis olmasini gosterir.
  * `uniprot` kaynakli bir atif, o makalenin bu tipin substrat atamasini
    destekledigini GOSTERMEZ. Kaynak alani okunmadan hicbir kayit
    "bu enzimin makalesi" diye okunmamalidir.

AG. `--fetch` olmadan betik aga CIKMAZ; onbellek yetkilidir. Onbellek
istek basina bir dosyadir, adi sorulan kimligin kendisinden gelir (sira
degil), cunku sira ile adlandirilan onbellek is listesi degistiginde
sessizce yanlis dosyayi okur.

Cikti: analysis_out/literature.json
       analysis_out/literature_cache/*.json  (onbellek, --fetch ile)
"""

import argparse
import csv
import glob
import json
import os
import re
import time
import urllib.error
import urllib.request

import cross_references

UA = "ROAR-DB/1.0 (literature.py; mailto:noreply@example.org)"

#  Kuratorluk metninden tanimlayici cikarma. DOI deseni kasten genis:
#  kuratorun cumlesi noktalama ile bitebiliyor, o yuzden sagdan temizlenir.
DOI_RE = re.compile(r"\b10\.\d{4,9}/[^\s,;)\]\"']+", re.I)
PMID_RE = re.compile(r"\bPMID:?\s*(\d{6,8})\b", re.I)
UNIPROT_RE = re.compile(r"UniProt[:\s]+([A-Z0-9]{6,10})")
PDB_RE = re.compile(r"\bPDB\s+([0-9][A-Za-z0-9]{3})")

#  RCSB alan adlari. Kucuk harfli varyantlari sessizce None donduruyor.
RCSB_DOI = "pdbx_database_id_DOI"
RCSB_PMID = "pdbx_database_id_PubMed"


def trim_doi(doi):
    """Kuratorun cumlesinden gelen noktalamayi atar."""
    return doi.rstrip(".,;:)]}'\"")


def valid_doi(doi):
    return bool(re.fullmatch(r"10\.\d{4,9}/\S+", doi or ""))


class Cache:
    """Istek basina bir dosya, adi SORULAN KIMLIKTEN turetilir."""

    def __init__(self, path, fetch):
        self.path = path
        self.fetch = fetch
        os.makedirs(path, exist_ok=True)
        self.hits = 0
        self.misses = 0
        self.errors = {}

    def _file(self, key):
        safe = re.sub(r"[^A-Za-z0-9._-]", "_", key)
        return os.path.join(self.path, safe + ".json")

    def get(self, key, url, pause=0.34):
        f = self._file(key)
        if os.path.exists(f):
            self.hits += 1
            with open(f) as fh:
                return json.load(fh)
        if not self.fetch:
            self.misses += 1
            return None
        req = urllib.request.Request(url, headers={"User-Agent": UA})
        try:
            with urllib.request.urlopen(req, timeout=60) as resp:
                data = json.load(resp)
        except (urllib.error.URLError, json.JSONDecodeError, OSError) as err:
            #  Hata ONBELLEGE YAZILMAZ. Yazilsa bir kerelik 503, kalici
            #  "bu kaydin atifi yok" haline gelirdi.
            self.errors[key] = str(err)
            return None
        tmp = f + ".partial"
        with open(tmp, "w") as fh:
            json.dump(data, fh)
        os.replace(tmp, f)
        time.sleep(pause)
        return data


def anchors(root):
    """Her tip icin: kuratorluk metni, UniProt erisimleri, PDB kimlikleri."""
    chem_path = os.path.join(root, "chemistry.csv")
    out = {}
    with open(chem_path) as fh:
        for row in csv.DictReader(fh):
            cluster = row["cluster"]
            text = " ".join(filter(None, (row.get("source"), row.get("notes"))))
            yaml_path = os.path.join(root, "types", cluster + ".yaml")
            if os.path.exists(yaml_path):
                with open(yaml_path) as yh:
                    text += " " + yh.read()
            ups = {m for m in UNIPROT_RE.findall(text)
                   if cross_references.PATTERNS["uniprot"].match(m)}
            pdbs = set(PDB_RE.findall(text))
            if (row.get("pdb") or "").strip():
                pdbs.add(row["pdb"].strip())
            pdbs = {p.upper() for p in pdbs
                    if re.fullmatch(r"[0-9][A-Za-z0-9]{3}", p)}
            out[cluster] = {
                "text": text,
                "source_kind": row.get("source_kind") or "",
                "uniprot": sorted(ups),
                "pdb": sorted(pdbs),
            }
    return out


def curated_records(cluster, info):
    """Kendi metnimizde yazili tanimlayicilar, kuratorun cumlesiyle."""
    text = info["text"]
    found = {}
    for match in DOI_RE.finditer(text):
        doi = trim_doi(match.group(0))
        if not valid_doi(doi):
            continue
        found.setdefault(("doi", doi.lower()),
                         (match.start(), match.start() + len(doi)))
    for match in PMID_RE.finditer(text):
        found.setdefault(("pmid", match.group(1)), (match.start(), match.end()))

    #  CUMLE SINIRI. Ilk yazim ilk noktaya kadar kesiyordu ve nokta
    #  DOI'nin ICINDE kaliyordu: "doi:10.1016/j.checat." gibi bir baglam
    #  okuyucuya yardim etmek yerine yanlis bir tanimlayici gosteriyordu.
    #  Cumle sonu artik "nokta + bosluk ya da metin sonu" olarak okunuyor,
    #  yani tanimlayicilarin icindeki noktalar sinir sayilmiyor.
    boundary = re.compile(r"\.(?=\s|$)")
    stops = [m.start() for m in boundary.finditer(text)]

    records = []
    for (kind, value), span in sorted(found.items(), key=lambda kv: kv[1]):
        at, at_end = span
        #  Tanimlayicinin GECTIGI cumle saklanir: bir DOI kuratorluk
        #  metninde bir DUZELTMENIN gerekcesi olarak da gecebiliyor, ve o
        #  zaman "bu tipin makalesi" degildir. Cumleyi yayinlamak, kaydi
        #  okuyanin bunu kendisi gormesini saglar.
        before = [i for i in stops if i < at]
        after = [i for i in stops if i >= at_end]
        left = (before[-1] + 1) if before else 0
        right = after[0] if after else len(text) - 1
        sentence = text[left:right + 1].strip()
        records.append({
            "source": "curated",
            "doi": value if kind == "doi" else None,
            "pmid": value if kind == "pmid" else None,
            "context": re.sub(r"\s+", " ", sentence)[:400],
        })
    return records


def pdb_records(cluster, info, cache):
    records = []
    for pdb in info["pdb"]:
        data = cache.get("pdb_" + pdb,
                         "https://data.rcsb.org/rest/v1/core/entry/" + pdb)
        if not data:
            continue
        for cit in (data.get("citation") or []):
            if cit.get("id") != "primary" and cit.get("rcsb_is_primary") != "Y":
                continue
            doi = cit.get(RCSB_DOI)
            pmid = cit.get(RCSB_PMID)
            records.append({
                "source": "pdb_primary",
                "anchor": pdb,
                "doi": (doi or "").lower() or None,
                "pmid": str(pmid) if pmid else None,
                "title": cit.get("title"),
                "journal": cit.get("rcsb_journal_abbrev") or cit.get("journal_abbrev"),
                "year": cit.get("year"),
                "authors": (cit.get("rcsb_authors") or [])[:6],
                #  Atamanin yapidan geldigi tiplerde bu atif gerekcenin
                #  kendisidir; digerlerinde yalnizca yapiyi cozen makale.
                "is_assignment_basis": "structure" in info["source_kind"],
            })
    return records


def uniprot_records(cluster, info, cache):
    records = []
    for acc in info["uniprot"]:
        data = cache.get("uniprot_" + acc,
                         "https://rest.uniprot.org/uniprotkb/%s.json" % acc)
        if not data:
            continue
        #  Hangi atif bir ISLEV iddiasina kanit olarak bagli?
        evidence = set()
        for comment in (data.get("comments") or []):
            if comment.get("commentType") not in ("FUNCTION", "CATALYTIC ACTIVITY"):
                continue
            for text in (comment.get("texts") or []):
                for ev in (text.get("evidences") or []):
                    if ev.get("source") == "PubMed":
                        evidence.add(str(ev.get("id")))
            for ev in ((comment.get("reaction") or {}).get("evidences") or []):
                if ev.get("source") == "PubMed":
                    evidence.add(str(ev.get("id")))

        for ref in (data.get("references") or []):
            cit = ref.get("citation") or {}
            xrefs = {x.get("database"): x.get("id")
                     for x in (cit.get("citationCrossReferences") or [])}
            pmid = xrefs.get("PubMed")
            doi = (xrefs.get("DOI") or "").lower() or None
            records.append({
                "source": "uniprot",
                "anchor": acc,
                "doi": doi,
                "pmid": str(pmid) if pmid else None,
                "title": cit.get("title"),
                "journal": cit.get("journal"),
                "year": (cit.get("publicationDate") or "")[:4] or None,
                "authors": (cit.get("authors") or [])[:6],
                "function_evidence": bool(pmid and str(pmid) in evidence),
            })
    return records


def resolve_dois(records, cache):
    """Yayinlanan her DOI SICILDE kayitli mi?

    AYIRT EDICI DOGRULANDI. Ilk yazim `HEAD https://doi.org/<doi>` ile
    yayinciya gidiyordu ve 125 DOI'nin 79'unu "cozulmedi" sayiyordu. Bu
    rakam inandirici degildi: UniProt ve RCSB'den gelen yayinlanmis
    DOI'lerin ucte ikisi olu olamaz. Olculdu -- besinin besi de
    yayincidan HTTP 403 donuyor ama DOI handle sicilinden
    `responseCode: 1` donuyor, yani hepsi kayitli. Wiley, Elsevier ve ACS
    betikten gelen istegi engelliyor; olculen sey "bu DOI gercek mi"
    degil, "bu yayinci bota cevap veriyor mu"ydu.

    Soru artik SAHIBINE soruluyor: handle sicili. responseCode 1 kayitli,
    100 handle bulunamadi demektir. Yayincinin sayfasina hic gidilmiyor.
    """
    seen = {}
    for rec in records:
        doi = rec.get("doi")
        if not doi:
            continue
        if doi in seen:
            rec["resolved"] = seen[doi]
            continue
        #  Onbellek dosyasi DOGRUDAN okunuyor: bu istek bir HEAD, yani
        #  Cache.get'in JSON govdesi bekleyen yolu buraya uymuyor.
        f = cache._file("resolve_" + doi)
        if os.path.exists(f):
            with open(f) as fh:
                seen[doi] = json.load(fh).get("resolved")
            rec["resolved"] = seen[doi]
            continue
        if not cache.fetch:
            rec["resolved"] = None
            continue
        req = urllib.request.Request("https://doi.org/api/handles/" + doi,
                                     headers={"User-Agent": UA})
        ok = None
        try:
            with urllib.request.urlopen(req, timeout=45) as resp:
                ok = json.load(resp).get("responseCode") == 1
        except urllib.error.HTTPError as err:
            #  404 ve 100: handle yok, yani DOI kayitli DEGIL. Digerleri
            #  sicilin kendi hatasi, ve "kayitli degil" diye okunamaz.
            ok = False if err.code == 404 else None
        except (urllib.error.URLError, OSError, json.JSONDecodeError):
            ok = None
        if ok is not None:
            tmp = f + ".partial"
            with open(tmp, "w") as fh:
                json.dump({"doi": doi, "resolved": ok}, fh)
            os.replace(tmp, f)
            time.sleep(0.25)
        seen[doi] = ok
        rec["resolved"] = ok
    return seen


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--root", default=".")
    ap.add_argument("--out", default="analysis_out/literature.json")
    ap.add_argument("--cache", default="analysis_out/literature_cache")
    ap.add_argument("--fetch", action="store_true",
                    help="aga cik; verilmezse yalnizca onbellek okunur")
    args = ap.parse_args()

    info = anchors(args.root)
    cache = Cache(args.cache, args.fetch)

    types = {}
    flat = []
    for cluster in sorted(info):
        recs = curated_records(cluster, info[cluster])
        recs += pdb_records(cluster, info[cluster], cache)
        recs += uniprot_records(cluster, info[cluster], cache)
        #  Ayni DOI birden fazla kaynaktan gelebilir. BIRLESTIRILMEZ:
        #  hangi yolla geldigi iddianin parcasi. Yalnizca ayni kaynak +
        #  ayni tanimlayici tekrari atilir.
        uniq = {}
        for rec in recs:
            key = (rec["source"], rec.get("anchor"), rec.get("doi"), rec.get("pmid"))
            uniq.setdefault(key, rec)
        recs = list(uniq.values())
        for rec in recs:
            rec["cluster"] = cluster
        types[cluster] = recs
        flat += recs

    resolved = resolve_dois(flat, cache)

    by_source = {}
    for rec in flat:
        by_source[rec["source"]] = by_source.get(rec["source"], 0) + 1
    with_any = sorted(c for c in types if types[c])
    with_doi = sorted(c for c in types
                      if any(r.get("doi") for r in types[c]))
    unresolved = sorted({r["doi"] for r in flat
                         if r.get("doi") and r.get("resolved") is False})
    unchecked = sorted({r["doi"] for r in flat
                        if r.get("doi") and r.get("resolved") is None})

    report = {
        "module": "literature.py",
        "generated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "what_was_measured":
            "Literature identifiers for each curated type, from THREE "
            "SEPARATE sources with the source recorded ON EACH RECORD: "
            "identifiers written in this project's own curation (curated), "
            "the primary citation of the type's PDB entry (pdb_primary), "
            "and the reference list of the reference protein's UniProt "
            "entry (uniprot).",
        "what_this_is_not":
            "A record appearing here does NOT mean that paper supports this "
            "type's substrate assignment. Only `curated` records are the "
            "basis this project cites; `uniprot` records are the entry's "
            "reference list, and only those marked `function_evidence` are "
            "attached by UniProt as evidence for a statement about what the "
            "enzyme does.",
        "coverage": {
            "types": len(types),
            "types_with_any_citation": len(with_any),
            "types_with_a_doi": len(with_doi),
            "types_with_none": len(types) - len(with_any),
            "records": len(flat),
            "records_by_source": by_source,
            "distinct_dois": len({r["doi"] for r in flat if r.get("doi")}),
        },
        "doi_resolution": {
            "checked": len([d for d in resolved if resolved[d] is not None]),
            "resolved": len([d for d in resolved if resolved[d] is True]),
            "unresolved": unresolved,
            "not_checked": unchecked,
            "method": "GET https://doi.org/api/handles/<doi>; responseCode 1 means the DOI is registered. The publisher's page is NOT opened: many publishers answer an automated request with an error, and reading that error as a dead DOI marked 79 of 125 valid DOIs as broken in the first version of this check",
        },
        "cache": {
            "path": args.cache,
            "hits": cache.hits,
            "misses_without_fetch": cache.misses,
            "fetch_errors": cache.errors,
        },
        "types": types,
    }
    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    tmp = args.out + ".partial"
    with open(tmp, "w") as fh:
        json.dump(report, fh, indent=1, ensure_ascii=False)
    os.replace(tmp, args.out)

    print("[literatur] %d tip, %d kayit, %d farkli DOI"
          % (len(types), len(flat), report["coverage"]["distinct_dois"]))
    print("  kaynaga gore:", by_source)
    print("  atfi olan tip: %d, DOI'si olan: %d, hicbiri olmayan: %d"
          % (len(with_any), len(with_doi), len(types) - len(with_any)))
    print("  DOI cozumu: %d bakildi, %d cozuldu, %d cozulmedi, %d bakilmadi"
          % (report["doi_resolution"]["checked"],
             report["doi_resolution"]["resolved"],
             len(unresolved), len(unchecked)))
    if cache.errors:
        print("  AG HATASI (onbellege yazilmadi):", len(cache.errors))
    if cache.misses and not args.fetch:
        print("  onbellekte olmayan %d istek atlandi (--fetch verilmedi)"
              % cache.misses)
    print("[yazildi]", args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
