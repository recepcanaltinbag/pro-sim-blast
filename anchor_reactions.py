"""Rhea kimligi olmayan tipler icin capa ara: EC numarasi ya da alternatif
UniProt kaydi uzerinden resmi reaksiyona baglan.

NEDEN. 53 tipin kimyasi yalnizca kurasyona dayaniyordu cunku UniProt kaydinda
katalitik aktivite yoktu ya da accession UniProt bicimli degildi. Ama bazi
tiplerin EC numarasi `source`/`notes` metninde duruyor ve EC'den Rhea'ya
dogrudan gidilebiliyor; bazilarinin da `accession_db` alani ayni diziyi
tasiyan bir UniProt kaydini adiyla aniyor.

CAPA UYDURULMAZ. Metinde EC yoksa ve alternatif UniProt yoksa tip `curated`
olarak KALIR -- uydurma bir EC, eksik EC'den kotudur.
"""
import json
import re
import sys
import time
import urllib.parse
import urllib.request

import yaml

UA = "ROAR-DB/1.0 (mailto:recocanextra@gmail.com)"
EC_RE = re.compile(r"\bEC[: ]\s*([1-6]\.[0-9]+\.[0-9]+\.[0-9]+)\b")
EC_BARE = re.compile(r"\b([1-6]\.[0-9]{1,2}\.[0-9]{1,2}\.[0-9]{1,3})\b")
UP_RE = re.compile(r"\b(?:[OPQ][0-9][A-Z0-9]{3}[0-9]"
                   r"|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})\b")


def get(url, tries=3):
    for i in range(tries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": UA})
            return urllib.request.urlopen(req, timeout=60).read()
        except Exception:
            if i == tries - 1:
                raise
            time.sleep(2 * (i + 1))


def rhea_by_ec(ec):
    url = ("https://www.rhea-db.org/rhea?query=" + urllib.parse.quote("ec:" + ec)
           + "&columns=rhea-id,equation&format=tsv")
    lines = [ln for ln in get(url).decode().splitlines() if ln.strip()]
    out = []
    for ln in lines[1:]:
        parts = ln.split("\t")
        if len(parts) >= 2 and parts[0].startswith("RHEA:"):
            out.append({"rhea": parts[0], "equation": parts[1]})
    return out


def uniprot_rhea(acc):
    d = json.loads(get(f"https://rest.uniprot.org/uniprotkb/{acc}.json"))
    out = []
    for c in d.get("comments") or []:
        if c.get("commentType") != "CATALYTIC ACTIVITY":
            continue
        r = c.get("reaction") or {}
        rid = next((x["id"] for x in (r.get("reactionCrossReferences") or [])
                    if x["database"] == "Rhea" and x["id"].startswith("RHEA:")),
                   None)
        if rid:
            out.append({"rhea": rid, "equation": r.get("name")})
    pd = d.get("proteinDescription") or {}
    ecs = [e["value"] for e in
           ((pd.get("recommendedName") or {}).get("ecNumbers") or [])]
    return out, ecs


def main():
    import glob
    rows = []
    for f in sorted(glob.glob("types/*.yaml")):
        m = yaml.safe_load(open(f))
        cl = m["cluster"]
        if m.get("chemistry_evidence") != "curated":
            continue
        row = {"cluster": cl, "anchor": None, "found": []}
        text = " ".join(str(m.get(k) or "") for k in
                        ("source", "notes", "reaction", "accession_db"))
        ecs = EC_RE.findall(text) or [e for e in EC_BARE.findall(text)
                                      if e.count(".") == 3
                                      and not e.startswith(("0.", "7.", "8.",
                                                            "9."))]
        alt = [u for u in UP_RE.findall(str(m.get("accession_db") or ""))
               if u != str(m.get("accession"))]
        try:
            if ecs:
                hits = rhea_by_ec(ecs[0])
                if hits:
                    row.update(anchor=f"EC {ecs[0]} (from the record text)",
                               found=hits)
            if not row["found"] and alt:
                hits, aecs = uniprot_rhea(alt[0])
                if hits:
                    row.update(anchor=f"UniProt {alt[0]} (named in "
                                      f"accession_db as holding the same "
                                      f"sequence)", found=hits)
                elif aecs:
                    hits = rhea_by_ec(aecs[0])
                    if hits:
                        row.update(anchor=f"UniProt {alt[0]} -> EC {aecs[0]}",
                                   found=hits)
        except Exception as exc:                        # noqa: BLE001
            row["error"] = str(exc)[:90]
        rows.append(row)
        mark = "CAPA" if row["found"] else "yok "
        print(f"[{mark}] {cl}" + (f" <- {row['anchor']}" if row["found"]
                                  else ""), file=sys.stderr)
        time.sleep(0.2)
    json.dump({"_about": "curated tipler icin resmi reaksiyon capasi aramasi. "
                         "Capa uydurulmaz; bulunamayan tip curated kalir.",
               "rows": rows},
              open("analysis_out/reaction_anchors.json", "w"), indent=1,
              ensure_ascii=False)
    got = [r for r in rows if r["found"]]
    print(f"\ncurated tip: {len(rows)} | CAPA BULUNAN: {len(got)} | "
          f"capasiz kalan: {len(rows) - len(got)}")


if __name__ == "__main__":
    main()
