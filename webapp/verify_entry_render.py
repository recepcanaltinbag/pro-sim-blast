"""Komsuluk SVG'sinin JS PORTUNU Python aslina karsi dogrula.

NEDEN BU BETIK VAR
  Statik sitenin agirligini dusurmek icin komsuluk gen haritasi artik sunucuda
  HTML'e gomulmuyor, tarayicida ciziliyor (bkz build_entry_data.py). Bu, ayni
  mantigin IKI dilde yasamasi demek: app.neighborhood_svg (Python) ve
  static/neighsvg.js (JS). Iki kopya sessizce ayrisabilir -- ve bu depoda tam
  bu sinif hata birkac kez olculdu (motif kolonlarinin iki yerde tanimli
  olmasi, sozlugun yazici ile dogrulayicida ayrismasi).

  Burada ayrisma SESSIZ KALAMAZ: ayni girdiyle iki tarafin urettigi dizgi
  BIREBIR karsilastirilir. Ornekleme degil; --all ile butun dogrulanmis
  girisler taranir.

  Cikis kodu 1: en az bir giriste dizgiler farkli.

Kullanim:
    python3 verify_entry_render.py --limit 500      # hizli kontrol
    python3 verify_entry_render.py --all            # yayin oncesi
"""
import argparse
import json
import os
import subprocess
import sqlite3
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

NODE_RUNNER = r"""
const fs = require('fs');
const m = require(process.argv[2]);
const jobs = JSON.parse(fs.readFileSync(process.argv[3], 'utf8'));
const out = {};
for (const j of jobs.jobs) {
  out[j.id] = m.neighborhoodSvg(j.ro, j.genes, jobs.colors, jobs.window);
}
fs.writeFileSync(process.argv[4], JSON.stringify(out));
"""

GENE_FIELDS = ["start", "end", "strand", "component", "category", "gene",
               "product", "distance", "spans_origin"]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--db", default=os.path.join(os.path.dirname(HERE), "roar.sqlite"))
    ap.add_argument("--limit", type=int, default=500)
    ap.add_argument("--all", action="store_true", help="butun dogrulanmis girisler")
    ap.add_argument("--js", default=os.path.join(HERE, "static", "neighsvg.js"))
    ap.add_argument("--show", type=int, default=3, help="kac farki ayrintili yaz")
    a = ap.parse_args()

    import app as A                                   # noqa: PLC0415

    con = sqlite3.connect(f"file:{a.db}?mode=ro", uri=True)
    con.row_factory = sqlite3.Row
    q = "SELECT candidate_id FROM ro WHERE is_confirmed=1 ORDER BY candidate_id"
    if not a.all:
        q += f" LIMIT {a.limit}"
    ids = [r[0] for r in con.execute(q)]
    print(f"[giris] {len(ids)}")

    jobs, expected = [], {}
    for n, cid in enumerate(ids, 1):
        d = A.ro_detail(con, cid)
        if not d:
            continue
        ro, nbs, op_ids = d["ro"], d["neighbors"], d["operon_ids"]
        expected[cid] = A.neighborhood_svg(ro, nbs, op_ids)
        genes = []
        for nb in nbs:
            g = {f: (nb[f] if f in nb.keys() else None) for f in GENE_FIELDS}
            g["in_operon"] = 1 if nb["neighbor_id"] in op_ids else 0
            genes.append(g)
        jobs.append({"id": cid, "genes": genes,
                     "ro": {"start": ro["start"], "end": ro["end"], "strand": ro["strand"],
                            "gene": ro["gene"], "product": ro["product"]}})
        if n % 2000 == 0:
            print(f"   python {n}/{len(ids)}")
    con.close()

    with tempfile.TemporaryDirectory() as tmp:
        runner = os.path.join(tmp, "run.js")
        inp = os.path.join(tmp, "in.json")
        outp = os.path.join(tmp, "out.json")
        with open(runner, "w") as fh:
            fh.write(NODE_RUNNER)
        with open(inp, "w", encoding="utf-8") as fh:
            json.dump({"jobs": jobs, "colors": A.COMPONENT_COLORS,
                       "window": A.NEIGHBOR_WINDOW}, fh, ensure_ascii=False)
        print("[node] calisiyor...")
        r = subprocess.run(["node", runner, os.path.abspath(a.js), inp, outp],
                           capture_output=True, text=True)
        if r.returncode != 0:
            print("[HATA] node basarisiz:\n" + (r.stderr or "")[:2000])
            return 1
        with open(outp, encoding="utf-8") as fh:
            got = json.load(fh)

    bad = [cid for cid in expected if got.get(cid) != expected[cid]]
    print()
    if not bad:
        print(f"BIREBIR AYNI: {len(expected)} giriste Python ve JS ayni SVG'yi uretti.")
        return 0

    print(f"FARK VAR: {len(bad)}/{len(expected)} giriste dizgiler ayni degil.")
    for cid in bad[:a.show]:
        e, g = expected[cid], got.get(cid) or ""
        i = next((k for k in range(min(len(e), len(g))) if e[k] != g[k]), min(len(e), len(g)))
        print(f"\n  {cid}  (ilk fark {i}. karakterde; python {len(e)} B, js {len(g)} B)")
        print(f"    python: ...{e[max(0,i-45):i+45]}...")
        print(f"    js    : ...{g[max(0,i-45):i+45]}...")
    return 1


if __name__ == "__main__":
    sys.exit(main())
