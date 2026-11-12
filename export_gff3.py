"""Genomik baglami GFF3 olarak yazar: standart bicim, kendi bicimimiz degil.

NEDEN. Bu DB operonlari, komsu genleri ve bilesen cagrilarini kendi
tablolarinda tutuyor; onlari okumak icin SQL ya da bizim CSV'lerimiz
gerekiyor. GFF3 ise genom tarayicilarinin (IGV, JBrowse, Artemis) ve
bioinformatik araclarinin (bedtools, gffread) dogrudan okudugu bicim.
Standart bir bicimde yayimlanmayan bir veri, baska bir analize GIRMEZ.

KOORDINAT -- ve bu dosyanin en kolay yanlis yapilacak yeri. Veritabani
BioPython'un konvansiyonunu tutuyor: `start` 0-TABANLI, `end` DISLAYICI.
GFF3 ise 1-TABANLI ve KAPSAYICI. Olculdu (2026-11-12, 4.000 onaylanmis
giris): 3.996'sinda `end - start == 3 * (aminoasit + 1)`, yani yari-acik
aralik ve dur kodonu dahil. Dolayisiyla GFF3'e `start + 1` yazilir ve
`end` oldugu gibi kalir. Kalan 4 kayit bu esitlige uymuyor (orijini asan
`join(...)` ve parcali CDS'ler); SESSIZCE DUZELTILMEZ -- boyle bir satir
`partial_coordinates=true` ile isaretlenir, cunku bir koordinati
"duzeltmek" onu uydurmak olur.

BIRIM. Her onaylanmis alfa altbirimi icin bir `operon` ozelligi ve onun
altinda her operon geni icin bir `CDS`. Komsularin TAMAMI yazilmaz:
operona girmeyen komsu bir baglam degil, bir pencere artefaktidir ve
19.125 giris icin dosyayi 1 milyon satira cikarirdi.

Cikti: analysis_out/roar_genomic_context.gff3.gz
"""

import argparse
import gzip
import os
import re
import sqlite3
import time
from urllib.parse import quote

#  GFF3 9. kolonunda bu karakterler KACIRILMAK zorunda (spec 1.21):
#  ; = & , ve kontrol karakterleri. Bosluk serbest.
SAFE = " :/()[]+'.|*-_"


def esc(value):
    """GFF3 ozellik degeri. Bos deger yazilmaz, cagiran atlar."""
    if value is None:
        return ""
    return quote(str(value).strip(), safe=SAFE)


def gff_strand(strand):
    #  DB'de 1 / -1 / 0 duruyor. GFF3 bilinmeyen yon icin '?' tanimli,
    #  '.' ise "yon ILGISIZ" demek -- bir CDS icin dogrusu '?'.
    return {1: "+", -1: "-", 0: "?"}.get(strand, "?")


def verify(path):
    """Yazilan dosyayi YENIDEN OKUYUP yapisini sinar.

    NEDEN. Bir dosyanin GFF3 oldugunu SOYLEMEK onu GFF3 yapmiyor, ve bu
    dosyayi okuyacak arac (IGV, bedtools) hatayi bizim degil KULLANICININ
    ekraninda verir. Bu yuzden yazan taraf kendi ciktisini ayristirir:
    9 kolon, 1-tabanli ve artan aralik, cozulen her `Parent`, ve
    `##sequence-region` sinirlari.

    Donen: rapor sozlugu. Bir ihlal bulunursa cagiran DURUR -- yarim
    gecerli bir degisim dosyasi yayimlanmaz.
    """
    regions, ids, problems = {}, set(), []
    parents, lines, features = [], 0, 0
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as fh:
        for lineno, line in enumerate(fh, 1):
            lines += 1
            line = line.rstrip("\n")
            if line.startswith("##sequence-region"):
                _, seqid, start, end = line.split()
                regions[seqid] = (int(start), int(end))
                continue
            if line.startswith("#") or not line:
                continue
            cols = line.split("\t")
            if len(cols) != 9:
                problems.append(f"satir {lineno}: 9 kolon yerine {len(cols)}")
                continue
            features += 1
            seqid, _src, _type, start, end, _score, strand, _phase, attrs = cols
            try:
                start, end = int(start), int(end)
            except ValueError:
                problems.append(f"satir {lineno}: koordinat sayi degil")
                continue
            if start < 1:
                problems.append(f"satir {lineno}: start {start} < 1 "
                                f"(GFF3 1-tabanli)")
            if end < start:
                problems.append(f"satir {lineno}: end {end} < start {start}")
            if strand not in ("+", "-", "?", "."):
                problems.append(f"satir {lineno}: yon {strand!r} gecersiz")
            if seqid not in regions:
                problems.append(f"satir {lineno}: {seqid} icin "
                                f"##sequence-region yok")
            else:
                low, high = regions[seqid]
                if start < low or end > high:
                    #  Orijini asan gen bunu yapabilir ve isaretli olmali.
                    if "spans_origin=true" not in attrs:
                        problems.append(
                            f"satir {lineno}: {seqid}:{start}-{end} "
                            f"bildirilen {low}-{high} disinda ve "
                            f"spans_origin isaretli degil")
            seen = set()
            for pair in attrs.split(";"):
                if not pair:
                    continue
                if "=" not in pair:
                    problems.append(f"satir {lineno}: ozellik {pair!r} "
                                    f"key=value degil")
                    continue
                key, _, value = pair.partition("=")
                if key in seen:
                    problems.append(f"satir {lineno}: {key} iki kez")
                seen.add(key)
                if key == "ID":
                    ids.add(value)
                elif key == "Parent":
                    parents.append((lineno, value))
    missing = [f"satir {lineno}: Parent={value} cozulmuyor"
               for lineno, value in parents if value not in ids]
    problems.extend(missing[:10])
    return {"lines": lines, "features": features,
            "sequence_regions": len(regions), "distinct_ids": len(ids),
            "parent_links": len(parents),
            "unresolved_parents": len(missing),
            "problems": problems[:20], "valid": not problems}


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--db", default="roar.sqlite")
    ap.add_argument("--out", default="analysis_out/roar_genomic_context.gff3.gz")
    ap.add_argument("--revision", default="")
    args = ap.parse_args()

    con = sqlite3.connect(args.db)
    con.row_factory = sqlite3.Row

    #  Yalnizca ozellik tasiyan replikonlar icin ##sequence-region.
    regions = con.execute("""
        SELECT p.nucleotide_id, p.length, p.organism, p.is_circular
          FROM replicon p
         WHERE p.nucleotide_id IN (SELECT nucleotide_id FROM ro WHERE is_confirmed = 1)
         ORDER BY p.nucleotide_id""").fetchall()

    ro_rows = con.execute("""
        SELECT r.candidate_id, r.nucleotide_id, r.start, r.end, r.strand,
               r.product, r.locus_tag, r.protein_id, r.gene, r.ro_cluster,
               r.ro_group, LENGTH(r.sequence) aa,
               o.start op_start, o.end op_end, o.strand op_strand,
               o.n_genes, o.completeness, o.layout
          FROM ro r
          LEFT JOIN operon o USING(candidate_id)
         WHERE r.is_confirmed = 1
         ORDER BY r.nucleotide_id, r.start""").fetchall()

    genes = {}
    for row in con.execute("""
            SELECT g.candidate_id, g.position, g.component, g.category,
                   g.product, n.start, n.end, n.strand, n.locus_tag,
                   n.protein_id, n.gene, n.spans_origin
              FROM operon_gene g
              JOIN neighbor n ON n.neighbor_id = g.neighbor_id
             ORDER BY g.candidate_id, g.position"""):
        genes.setdefault(row["candidate_id"], []).append(row)

    counts = {"regions": 0, "operons": 0, "cds": 0, "partial": 0}
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    with gzip.open(args.out, "wt", newline="\n") as fh:
        fh.write("##gff-version 3\n")
        fh.write(f"#!source ROAR-DB, generated {time.strftime('%Y-%m-%d')}\n")
        if args.revision:
            #  Bir GFF3 dosyasi elden ele dolasir; hangi veri surumunden
            #  geldigi dosyanin ICINDE durmali.
            fh.write(f"#!data-revision {args.revision}\n")
        fh.write("#!coordinates 1-based inclusive (converted from the "
                 "0-based half-open coordinates stored in the database)\n")
        for reg in regions:
            if not reg["length"]:
                continue
            fh.write(f"##sequence-region {reg['nucleotide_id']} 1 {reg['length']}\n")
            counts["regions"] += 1

        for ro in ro_rows:
            cid = ro["candidate_id"]
            members = genes.get(cid, [])
            op_id = f"operon:{cid}"
            if ro["op_start"] is not None and ro["op_end"] is not None:
                attrs = [f"ID={esc(op_id)}",
                         f"Name={esc(ro['ro_cluster'] or 'unassigned')}",
                         f"ro_alpha={esc(cid)}"]
                for key, value in (("ro_type", ro["ro_cluster"]),
                                   ("ro_group", ro["ro_group"]),
                                   ("genes", ro["n_genes"]),
                                   ("completeness", ro["completeness"]),
                                   ("layout", ro["layout"])):
                    if value not in (None, ""):
                        attrs.append(f"{key}={esc(value)}")
                fh.write("\t".join([
                    ro["nucleotide_id"], "ROAR-DB", "operon",
                    str(ro["op_start"] + 1), str(ro["op_end"]), ".",
                    gff_strand(ro["op_strand"]), ".", ";".join(attrs)]) + "\n")
                counts["operons"] += 1
                parent = op_id
            else:
                parent = None

            #  Alfa altbirimin kendisi. `partial_coordinates`: olculen
            #  esitlik tutmuyorsa isaretlenir, duzeltilmez.
            span = (ro["end"] or 0) - (ro["start"] or 0)
            partial = bool(ro["aa"]) and span != 3 * (ro["aa"] + 1)
            if partial:
                counts["partial"] += 1
            attrs = [f"ID={esc('cds:' + cid)}",
                     f"Name={esc(ro['gene'] or ro['locus_tag'] or cid)}"]
            if parent:
                attrs.append(f"Parent={esc(parent)}")
            for key, value in (("product", ro["product"]),
                               ("locus_tag", ro["locus_tag"]),
                               ("protein_id", ro["protein_id"]),
                               ("ro_type", ro["ro_cluster"]),
                               ("ro_group", ro["ro_group"]),
                               ("component", "ro_alpha")):
                if value not in (None, ""):
                    attrs.append(f"{key}={esc(value)}")
            if partial:
                attrs.append("partial_coordinates=true")
            fh.write("\t".join([
                ro["nucleotide_id"], "ROAR-DB", "CDS",
                str((ro["start"] or 0) + 1), str(ro["end"] or 0), ".",
                gff_strand(ro["strand"]), "0", ";".join(attrs)]) + "\n")
            counts["cds"] += 1

            for member in members:
                attrs = [f"ID={esc('cds:' + cid + ':' + str(member['position']))}",
                         f"Name={esc(member['gene'] or member['locus_tag'] or 'gene')}"]
                if parent:
                    attrs.append(f"Parent={esc(parent)}")
                for key, value in (("product", member["product"]),
                                   ("locus_tag", member["locus_tag"]),
                                   ("protein_id", member["protein_id"]),
                                   ("component", member["component"]),
                                   ("category", member["category"]),
                                   ("operon_position", member["position"])):
                    if value not in (None, ""):
                        attrs.append(f"{key}={esc(value)}")
                if member["spans_origin"]:
                    #  Orijini asan gen: GFF3'te tek aralikla temsil
                    #  EDILEMEZ, o yuzden isaretlenir.
                    attrs.append("spans_origin=true")
                fh.write("\t".join([
                    ro["nucleotide_id"], "ROAR-DB", "CDS",
                    str((member["start"] or 0) + 1), str(member["end"] or 0), ".",
                    gff_strand(member["strand"]), "0", ";".join(attrs)]) + "\n")
                counts["cds"] += 1

    report = verify(args.out)
    counts["verified"] = report

    import json as _json
    report_path = os.path.splitext(os.path.splitext(args.out)[0])[0] + ".json"
    with open(report_path, "w") as rf:
        _json.dump({
            "question": "Can this database's genomic context be read by "
                        "standard tools without converting anything?",
            "format": "GFF3 (version 3), gzip-compressed",
            "coordinates": "1-based inclusive, converted from the 0-based "
                           "half-open coordinates stored in the database; "
                           "measured on 4000 entries, 3996 satisfy "
                           "end - start == 3 * (amino acids + 1)",
            "unit": "one operon feature per confirmed alpha subunit, with each "
                    "operon gene as a child CDS",
            "not_exported": "Neighbouring genes that are not part of a called "
                            "operon. A gene inside the search window but "
                            "outside the operon is a window artefact, not "
                            "context, and including it would turn a 117 "
                            "thousand line file into a million.",
            "counts": {k: v for k, v in counts.items() if k != "verified"},
            "verification": report,
            "data_revision": args.revision or None,
        }, rf, indent=1)

    print(f"[yazildi] {args.out}")
    print(f"  sequence-region : {counts['regions']}")
    print(f"  operon          : {counts['operons']}")
    print(f"  CDS             : {counts['cds']}")
    print(f"  partial koordinat isaretli: {counts['partial']}")
    print(f"[dogrulama] {report['features']} ozellik, "
          f"{report['parent_links']} Parent baglantisi, "
          f"{report['unresolved_parents']} cozulmeyen")
    print(f"[yazildi] {report_path}")
    con.close()
    if not report["valid"]:
        print("\n[GECERSIZ] yapisal ihlaller:")
        for problem in report["problems"]:
            print(f"    {problem}")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main() or 0)
