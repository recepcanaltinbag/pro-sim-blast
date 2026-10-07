"""Veritabaninin her surumunu kaydeder ve iki surum arasinda NE DEGISTIGINI soyler.

NEDEN. Pipeline yeniden kosunca her sayi degisir ve "neden degisti" sorusunun
cevabi hicbir yerde durmuyordu. Bir tipin uye sayisi 94'ten 310'a ciktiginda
bunun sebebi yeni indirilen genomlar mi, bir kuratorluk duzeltmesi mi, yoksa
referans setine eklenen bir dizi mi -- gecmise donup bakmanin yolu yoktu.

Bu dosya her surumde uc sey saklar:

  SAYIM   dogrulanmis giris, tip, yaprak, replikon, komsu, operon, kanit
          duzeyi dagilimi, VE tip basina uye sayisi. Asil kiyas bu sonuncusu;
          fark alinca hangi tipin kac uye kazandigi/kaybettigi dogrudan
          gorulur.
  GIRDI   o surumu ureten seylerin parmak izi: chemistry.csv, referans
          fasta, profil kutuphanesi, ve git commit'i. Boylece bir degisiklik
          gorulunce sebebi aranacak yer bellidir.
  NOT     surumu kaydeden kisinin tek cumlelik aciklamasi.

Bir surum DEGISMEZ: dosya bir kez yazilir ve bir daha dokunulmaz. Kayitlar
`analysis_out/revisions/` altinda, ozet satirlar `revisions.jsonl` icinde.

Kullanim:
    python3 db_revision.py --note "2.962 genom yeniden indirildikten ONCE"
    python3 db_revision.py --list
    python3 db_revision.py --diff rev-0003 rev-0004
    python3 db_revision.py --diff last     # son iki surum
"""

import argparse
import glob
import hashlib
import json
import os
import re
import subprocess
import sqlite3
import sys
import time
from collections import Counter

REVDIR = "analysis_out/revisions"
INDEX = "analysis_out/revisions.jsonl"

# Bu dosyalar o surumun NEDENIDIR; parmak izleri kayda girer.
INPUT_FILES = [
    "chemistry.csv",
    "cluster_ecology.csv",
    "ROs_71_Clean/refs71.fasta",
    "RieskeDB71.hmm",
    "ROs_71_Clean/ROmotif71.hmm",
]

COUNT_TABLES = ["ro", "replicon", "replicon_source", "neighbor", "neighbor_protein",
                "operon", "operon_gene", "gene_category", "leaf", "ro_leaf",
                "ro_evidence", "ro_etc", "ro_regulation", "ro_carboxylate",
                "ro_domain", "ro_subfamily"]


def file_digest(path):
    if not os.path.exists(path):
        return None
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()[:16]


def git_commit():
    try:
        out = subprocess.run(["git", "rev-parse", "--short", "HEAD"],
                             capture_output=True, text=True, timeout=20)
        dirty = subprocess.run(["git", "status", "--porcelain"],
                               capture_output=True, text=True, timeout=20)
        return {"commit": out.stdout.strip() or None,
                "uncommitted_files": len([l for l in dirty.stdout.splitlines() if l])}
    except (subprocess.SubprocessError, OSError):
        return {"commit": None, "uncommitted_files": None}


def snapshot(db):
    con = sqlite3.connect(f"file:{db}?mode=ro", uri=True)
    have = {r[0] for r in con.execute(
        "SELECT name FROM sqlite_master WHERE type='table'")}
    counts = {}
    for t in COUNT_TABLES:
        if t in have:
            counts[t] = con.execute(f"SELECT COUNT(*) FROM {t}").fetchone()[0]

    confirmed = con.execute(
        "SELECT COUNT(*) FROM ro WHERE is_confirmed=1").fetchone()[0]
    per_type = dict(con.execute(
        "SELECT ro_cluster, COUNT(*) FROM ro WHERE is_confirmed=1 "
        "GROUP BY 1").fetchall())
    per_group = dict(con.execute(
        "SELECT ro_group, COUNT(*) FROM ro WHERE is_confirmed=1 "
        "GROUP BY 1").fetchall())
    tiers = dict(con.execute(
        "SELECT COALESCE(e.tier,'none'), COUNT(*) FROM ro r "
        "LEFT JOIN ro_evidence e USING(candidate_id) "
        "WHERE r.is_confirmed=1 GROUP BY 1").fetchall()) if "ro_evidence" in have else {}
    # Genom kapsamasi: bos indirilmis ve dizisiz replikonlar -- bu projenin
    # en buyuk bilinen eksikligi, her surumde olculmeli.
    empty_replicons = con.execute(
        "SELECT COUNT(*) FROM replicon WHERE cds_count IS NULL OR cds_count=0"
    ).fetchone()[0]
    con.close()
    return {
        "counts": counts,
        "confirmed_entries": confirmed,
        "types_with_members": len(per_type),
        "members_per_type": per_type,
        "members_per_group": {str(k): v for k, v in per_group.items()},
        "evidence_tiers": tiers,
        "replicons_without_annotation": empty_replicons,
    }


def next_id():
    os.makedirs(REVDIR, exist_ok=True)
    used = [int(m.group(1)) for p in glob.glob(os.path.join(REVDIR, "rev-*.json"))
            for m in [re.search(r"rev-(\d+)", os.path.basename(p))] if m]
    return "rev-%04d" % ((max(used) + 1) if used else 1)


def record(db, note):
    rev = next_id()
    data = {
        "revision": rev,
        "recorded_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "note": note,
        "database": db,
        "inputs": {p: file_digest(p) for p in INPUT_FILES},
        "git": git_commit(),
        **snapshot(db),
    }
    path = os.path.join(REVDIR, rev + ".json")
    with open(path, "w") as fh:
        json.dump(data, fh, indent=1, sort_keys=True)
    summary = {k: data[k] for k in ("revision", "recorded_utc", "note",
                                    "confirmed_entries", "types_with_members",
                                    "replicons_without_annotation")}
    summary["git"] = data["git"]["commit"]
    with open(INDEX, "a") as fh:
        fh.write(json.dumps(summary) + "\n")
    print(f"[kaydedildi] {path}")
    print(f"  {rev}  {data['confirmed_entries']:,} dogrulanmis giris, "
          f"{data['types_with_members']} tip, "
          f"{data['replicons_without_annotation']:,} anotasyonsuz replikon")
    print(f"  not: {note}")
    return 0


def load(rev):
    path = os.path.join(REVDIR, rev + ".json")
    if not os.path.exists(path):
        raise SystemExit(f"surum yok: {path}")
    with open(path) as fh:
        return json.load(fh)


def diff(a_id, b_id):
    a, b = load(a_id), load(b_id)
    print(f"{a_id}  ->  {b_id}")
    print(f"  {a['recorded_utc']}  ->  {b['recorded_utc']}")
    print(f"  {a['note']}")
    print(f"  {b['note']}")
    print()

    # 1. NEDEN: hangi girdi degisti
    changed_in = [p for p in INPUT_FILES
                  if (a["inputs"].get(p) != b["inputs"].get(p))]
    print("Degisen girdiler (sebep burada aranir):")
    if changed_in:
        for p in changed_in:
            print(f"  {p}: {a['inputs'].get(p)} -> {b['inputs'].get(p)}")
    else:
        print("  hicbiri -- ayni kuratorlu veri ve ayni profil kutuphanesi")
    if a["git"].get("commit") != b["git"].get("commit"):
        print(f"  git: {a['git'].get('commit')} -> {b['git'].get('commit')}")
    print()

    # 2. NE: toplamlar
    def line(label, x, y, width=34):
        d = y - x
        arrow = f"{d:+,}" if d else "degismedi"
        pct = f"  ({100 * d / x:+.1f} %)" if x and d else ""
        print(f"  {label:<{width}} {x:>9,} -> {y:>9,}   {arrow}{pct}")

    print("Toplamlar:")
    line("dogrulanmis giris", a["confirmed_entries"], b["confirmed_entries"])
    line("uyesi olan tip", a["types_with_members"], b["types_with_members"])
    line("anotasyonsuz replikon", a["replicons_without_annotation"],
         b["replicons_without_annotation"])
    for t in COUNT_TABLES:
        if t in a["counts"] or t in b["counts"]:
            x, y = a["counts"].get(t, 0), b["counts"].get(t, 0)
            if x != y:
                line(t, x, y)
    print()

    tiers = sorted(set(a["evidence_tiers"]) | set(b["evidence_tiers"]))
    if tiers:
        print("Kanit duzeyi:")
        for t in tiers:
            x, y = a["evidence_tiers"].get(t, 0), b["evidence_tiers"].get(t, 0)
            if x != y:
                line(t, x, y)
        print()

    # 3. NEREDE: tip basina -- asil kiyas bu
    pa, pb = a["members_per_type"], b["members_per_type"]
    keys = set(pa) | set(pb)
    moved = sorted(((pb.get(k, 0) - pa.get(k, 0), k) for k in keys),
                   key=lambda kv: -abs(kv[0]))
    moved = [(d, k) for d, k in moved if d]
    gained_type = [k for k in keys if k not in pa]
    lost_type = [k for k in keys if k not in pb]
    print(f"Tip basina uye ({len(moved)} tip degisti):")
    for d, k in moved[:20]:
        tag = "  YENI" if k in gained_type else ("  BOSALDI" if k in lost_type else "")
        print(f"  {k:<22} {pa.get(k, 0):>7,} -> {pb.get(k, 0):>7,}   {d:+,}{tag}")
    if len(moved) > 20:
        print(f"  ... ve {len(moved) - 20} tip daha")
    if gained_type:
        print(f"\n  Uye KAZANAN bos tipler: {', '.join(sorted(gained_type))}")
    if lost_type:
        print(f"\n  Uyesini KAYBEDEN tipler: {', '.join(sorted(lost_type))}")
    return 0


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--db", default="roar.sqlite")
    ap.add_argument("--note", help="bu surumu kaydet, verilen aciklamayla")
    ap.add_argument("--list", action="store_true", help="surumleri listele")
    ap.add_argument("--diff", nargs="*", metavar="REV",
                    help="iki surumu karsilastir; 'last' son ikisini alir")
    args = ap.parse_args()

    if args.list:
        if not os.path.exists(INDEX):
            print("henuz surum yok")
            return 0
        for line in open(INDEX):
            r = json.loads(line)
            print(f"  {r['revision']}  {r['recorded_utc']}  "
                  f"{r['confirmed_entries']:>7,} giris  "
                  f"{r['types_with_members']:>3} tip  "
                  f"{r.get('git') or '-':>8}  {r['note']}")
        return 0

    if args.diff is not None:
        revs = args.diff
        if not revs or revs == ["last"]:
            found = sorted(glob.glob(os.path.join(REVDIR, "rev-*.json")))
            if len(found) < 2:
                raise SystemExit("karsilastirmak icin en az iki surum gerekli")
            revs = [os.path.basename(found[-2])[:-5],
                    os.path.basename(found[-1])[:-5]]
        if len(revs) != 2:
            raise SystemExit("--diff iki surum adi ister (ya da 'last')")
        return diff(revs[0], revs[1])

    if args.note:
        return record(args.db, args.note)
    ap.error("--note, --list ya da --diff gerekli")


if __name__ == "__main__":
    sys.exit(main())
