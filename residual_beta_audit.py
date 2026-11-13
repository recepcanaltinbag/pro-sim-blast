"""Anotasyonun beta dedigi ama siniflanmayan 60 proteini INCELER.

NEDEN BU BETIK VAR. `check_component_scan_is_fresh` kapisinin ikinci
ayagi, anotasyonun "dioxygenase small/beta subunit" dedigi her komsunun
`beta` siniflanmasini arar ve 60 proteinlik bir ARTAKALAN butcesi tasiyor
(duzeltme oncesi 688). ROADMAP o 60'i soyle tarif ediyordu: "urun metni
beta diyor, hicbir model tutmuyor; ikisinden biri yanlis ve hangisi
oldugu HENUZ OLCULMEDI".

**ILK OLCUM O TARIFI DEGISTIRDI.** 60'in TAMAMINDA `neighbor_protein`
satiri YOK -- yani bu proteinler bilesen taramasinin aradigi kumede hic
bulunmadi. Dolayisiyla "hicbir model tutmuyor" cumlesi ne anotasyon ne
model hakkinda bir ifade: TARANMAMIS bir dizi hakkinda hicbir model bir
sey soyleyemez. Soru "hangisi yanlis" degil, "dizi neden yok".

Bu betik kaynaga bakar: GenBank kaydindaki CDS'in `/pseudo`,
`/translation` ve not alanlarini okur. Koordinat donusumu olcumun
parcasi -- DB 0-tabanli yari-acik, GenBank 1-tabanli kapsayici
(`start+1`), ve bu donusum yapilmadan kayit bulunamiyor.

HICBIR SEY DEGISTIRMEZ: siniflandirma, kapi butcesi ve tablolar aynidir.

Cikti: analysis_out/residual_beta_audit.json
"""

import argparse
import gzip
import json
import os
import re
import sqlite3
import sys
from collections import Counter

KEY = ("nb.nucleotide_id||':'||nb.start||'-'||nb.end||':'||nb.strand")


def open_record(gbk_dir, nucleotide_id):
    plain = os.path.join(gbk_dir, nucleotide_id + ".gbk")
    if os.path.exists(plain):
        return open(plain, "rt", encoding="utf-8", errors="replace")
    if os.path.exists(plain + ".gz"):
        return gzip.open(plain + ".gz", "rt", encoding="utf-8",
                         errors="replace")
    return None


def read_feature(gbk_dir, nucleotide_id, start, end):
    """CDS/gene blogunu bulur ve nitelikleri dondurur.

    `start+1` cunku DB BioPython'un 0-tabanli yari-acik araligini
    saklarken GenBank 1-tabanli kapsayici yazar. Bu donusum yapilmadan
    kayit BULUNAMIYOR -- ilk denemede tam olarak bu oldu.
    """
    handle = open_record(gbk_dir, nucleotide_id)
    if handle is None:
        return {"record": "missing"}
    want = f"{start + 1}..{end}"
    found, attrs, block = False, {"record": "read"}, []
    with handle:
        for line in handle:
            if line.startswith("     ") and not line.startswith("      "):
                if found and block:
                    break
                found = want in line
                if found:
                    block = [line]
                continue
            if found:
                block.append(line)
    text = "".join(block)
    attrs["pseudo"] = "/pseudo" in text
    attrs["has_translation"] = "/translation=" in text
    note = re.search(r'/note="([^"]+)"', text, re.S)
    if note:
        attrs["note"] = " ".join(note.group(1).split())[:160]
    attrs["frameshift_or_partial"] = bool(
        re.search(r"frameshift|internal stop|incomplete|partial",
                  text, re.I))
    return attrs


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--db", default="roar.sqlite")
    parser.add_argument("--gbk", default="gbk_files")
    parser.add_argument("--out",
                        default="analysis_out/residual_beta_audit.json")
    args = parser.parse_args()

    con = sqlite3.connect(args.db)
    rows = con.execute(f"""
        SELECT DISTINCT {KEY} AS pk, nb.nucleotide_id, nb.start, nb.end,
               nb.product, nb.locus_tag
          FROM neighbor nb
         WHERE (lower(nb.product) LIKE '%dioxygenase small subunit%'
             OR lower(nb.product) LIKE '%dioxygenase beta subunit%')
           AND {KEY} NOT IN (SELECT protein_key FROM neighbor_component
                              WHERE component = 'beta')
         ORDER BY pk""").fetchall()
    print(f"[artakalan] {len(rows)} protein", flush=True)

    tally, detail = Counter(), []
    for pk, nuc, start, end, product, locus in rows:
        has_row = con.execute(
            "SELECT 1 FROM neighbor_protein WHERE nucleotide_id=? "
            "AND start=? AND end=?", (nuc, start, end)).fetchone() is not None
        attrs = read_feature(args.gbk, nuc, start, end)
        if attrs.get("pseudo"):
            verdict = "pseudogene"
        elif not attrs.get("has_translation"):
            verdict = "no_translation_in_record"
        elif not has_row:
            verdict = "translation_exists_but_not_extracted"
        else:
            verdict = "scanned_and_no_model_hit"
        tally[verdict] += 1
        detail.append({"protein_key": pk, "locus_tag": locus,
                       "product": product, "verdict": verdict,
                       "pseudo": attrs.get("pseudo"),
                       "frameshift_or_partial":
                           attrs.get("frameshift_or_partial"),
                       "note": attrs.get("note")})

    report = {
        "question": ("the gate carries a residual budget of 60 proteins the "
                     "annotation calls a beta subunit but which are not "
                     "classified as beta. Why?"),
        "what_the_roadmap_said": (
            "that one of the two was wrong -- the annotation or the model "
            "set -- and that which one had not been measured"),
        "what_the_measurement_says": (
            "neither. None of the 60 has a row in neighbor_protein, so no "
            "model was ever given a sequence to match. A statement about an "
            "unscanned sequence is not a statement about the annotation"),
        "method": {
            "source": "the CDS feature in the GenBank record itself",
            "coordinates": ("the database stores BioPython's 0-based "
                            "half-open range while GenBank writes 1-based "
                            "inclusive, so the feature is looked up at "
                            "start+1..end; without that conversion the "
                            "record is not found at all, which is what "
                            "happened on the first attempt"),
            "changes_nothing": ("classification, gate budget and tables are "
                                "untouched; this measures and reports"),
        },
        "verdicts": dict(tally.most_common()),
        "proteins": detail,
    }

    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    with open(args.out, "w", encoding="utf-8") as handle:
        json.dump(report, handle, indent=2, ensure_ascii=False)
        handle.write("\n")

    print(f"[hukum] {dict(tally.most_common())}")
    print(f"[yazildi] {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
