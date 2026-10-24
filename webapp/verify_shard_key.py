"""Parca anahtarini (md5) JS ile Python arasinda dogrula.

NEDEN. build_entry_data.py girisleri md5(candidate_id)'nin ilk iki hex
hanesine gore parcaliyor; tarayici ayni anahtari uretemezse giris 404 olur.
Iki ayri dilde iki md5 var, ve sessizce ayrisabilirler. Burada GERCEK
kimliklerin hepsi uzerinde karsilastirilir.

Cikis kodu 1: en az bir kimlikte anahtarlar farkli.
"""
import argparse, hashlib, json, os, sqlite3, subprocess, sys, tempfile

HERE = os.path.dirname(os.path.abspath(__file__))

RUNNER = """
const fs = require('fs');
const md5 = require(process.argv[2]);
const ids = JSON.parse(fs.readFileSync(process.argv[3], 'utf8'));
fs.writeFileSync(process.argv[4], JSON.stringify(ids.map(md5)));
"""


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--db", default=os.path.join(os.path.dirname(HERE), "roar.sqlite"))
    ap.add_argument("--js", default=os.path.join(HERE, "static", "md5.js"))
    a = ap.parse_args()

    con = sqlite3.connect(f"file:{a.db}?mode=ro", uri=True)
    ids = [r[0] for r in con.execute(
        "SELECT candidate_id FROM ro WHERE is_confirmed=1 ORDER BY candidate_id")]
    con.close()
    print(f"[kimlik] {len(ids)}")

    with tempfile.TemporaryDirectory() as t:
        r_js, inp, outp = (os.path.join(t, n) for n in ("r.js", "in.json", "out.json"))
        open(r_js, "w").write(RUNNER)
        json.dump(ids, open(inp, "w", encoding="utf-8"), ensure_ascii=False)
        r = subprocess.run(["node", r_js, os.path.abspath(a.js), inp, outp],
                           capture_output=True, text=True)
        if r.returncode:
            print("[HATA] node:", (r.stderr or "")[:1500])
            return 1
        got = json.load(open(outp, encoding="utf-8"))

    bad = [(i, hashlib.md5(i.encode()).hexdigest(), g)
           for i, g in zip(ids, got) if hashlib.md5(i.encode()).hexdigest() != g]
    if bad:
        print(f"FARK VAR: {len(bad)}/{len(ids)}")
        for i, p, g in bad[:5]:
            print(f"   {i}\n     python={p}\n     js    ={g}")
        return 1
    print(f"BIREBIR AYNI: {len(ids)} kimligin hepsinde md5 ayni.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
