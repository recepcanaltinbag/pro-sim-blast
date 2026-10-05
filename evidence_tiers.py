"""
Kanit duzeyi (evidence tier) -- "benzerlik var" ile "tanimlanmis enzim" ayni sey degil.

HMM atamasi bir girisi en yakin referans TIPINE koyar; bu, girisin o referansin
substratini cevirdigi anlamina gelmez. Burada her dogrulanmis RO'nun en yakin
KURATORLU REFERANS PROTEINE (71 deneysel RO) diamond blastp kimligi olculur ve
kademelendirilir:

    characterized   >= 95% kimlik, >= 90% kaplama  -- referans enzimin kendisi / sus varyanti
    close_homolog   >= 60%                          -- ayni fonksiyon cok olasi
    family_member   40-60%                          -- ayni tip, substrat belirsiz
    distant         25-40%                          -- RO alpha, en yakin tip X, fonksiyon bilinmiyor
    novel           < 25% veya hit yok              -- bilinen tiplere uymuyor

Esikler genel enzim-fonksiyon aktarim kurallarina dayanir (>=60% kimlikte EC
transferi genelde guvenilir; 40%'in altinda zayif) ve --cut-* ile degistirilebilir.

ONEMLI SINIR (analysis_out/reference_pairs.csv): 71 kuratorlu referansin kendi
aralarinda, FARKLI substrat etiketli ciftler %99,8 kimlige kadar cikiyor
(EdoA1 etilbenzen / CumA1 kumen), naftalen vs nitrotoluen dioksijenazlari ~%90.
Yani RO'larda substrat secimi birkac aktif-bolge kalintisiyla belirlenir; HICBIR
global kimlik esigi "ayni substrat" garantisi vermez. Kademeler "ayni enzim TIPI"
guvenini olcer; substrat etiketi her kademede referansin substratidir, uyenin degil.

Cikti: ro_evidence tablosu + analysis_out/evidence_by_cluster.csv
"""

import argparse
import csv
import os
import sqlite3
import subprocess
import tempfile
from collections import Counter, defaultdict

TIERS = [("characterized", 95.0), ("close_homolog", 60.0), ("family_member", 40.0),
         ("distant", 25.0)]
TIER_ORDER = ["characterized", "close_homolog", "family_member", "distant", "novel"]
MIN_COVER_CHARACTERIZED = 90.0


def tier_of(identity, qcov):
    if identity is None:
        return "novel"
    for name, cut in TIERS:
        if identity >= cut and (name != "characterized" or qcov >= MIN_COVER_CHARACTERIZED):
            return name
    return "novel"


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--db", default="roar.sqlite")
    ap.add_argument("--refs", default="ROs_71_Clean/refs71.fasta")
    ap.add_argument("--out-dir", default="analysis_out")
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--cut-characterized", type=float, default=95.0)
    ap.add_argument("--cut-close", type=float, default=60.0)
    ap.add_argument("--cut-family", type=float, default=40.0)
    ap.add_argument("--cut-distant", type=float, default=25.0)
    args = ap.parse_args()
    TIERS[:] = [("characterized", args.cut_characterized), ("close_homolog", args.cut_close),
                ("family_member", args.cut_family), ("distant", args.cut_distant)]

    con = sqlite3.connect(args.db)
    con.executescript("""
        CREATE TABLE IF NOT EXISTS ro_evidence (
            candidate_id    TEXT PRIMARY KEY,
            nearest_ref     TEXT,     -- en yakin kuratorlu referans (model adi)
            ref_identity    REAL,     -- % kimlik (diamond)
            ref_qcov        REAL,     -- sorgu kaplamasi %
            same_as_hmm     INTEGER,  -- en yakin referans, HMM'in atadigi kume ile ayni mi
            tier            TEXT
        );""")

    with tempfile.TemporaryDirectory() as tmp:
        q = os.path.join(tmp, "q.fa")
        with open(q, "w") as fh:
            for cid, seq in con.execute(
                    "SELECT candidate_id, sequence FROM ro WHERE is_confirmed=1 AND sequence IS NOT NULL"):
                fh.write(f">{cid}\n{seq}\n")
        db = os.path.join(tmp, "refs")
        subprocess.run(["diamond", "makedb", "--in", args.refs, "-d", db, "--quiet"], check=True)
        out = os.path.join(tmp, "hits.tsv")
        subprocess.run(["diamond", "blastp", "-q", q, "-d", db, "-o", out, "--quiet",
                        "-p", str(args.threads), "-k", "1", "--max-hsps", "1",
                        "--ultra-sensitive", "-e", "1e-5",
                        "--outfmt", "6", "qseqid", "sseqid", "pident", "qcovhsp", "bitscore"],
                       check=True)
        best = {}
        with open(out) as fh:
            for line in fh:
                qid, sid, pid, qcov, bits = line.rstrip("\n").split("\t")
                best[qid] = (sid, float(pid), float(qcov))

    rows = []
    hmm_cluster = dict(con.execute("SELECT candidate_id, ro_cluster FROM ro WHERE is_confirmed=1"))
    for cid, cluster in hmm_cluster.items():
        sid, pid, qcov = best.get(cid, (None, None, 0.0))
        # referans basligi "1_101_OxoO_monooxygenase_pro" -> model adi ilk 3 parca
        ref_model = "_".join(sid.split("_")[:3]) if sid else None
        rows.append((cid, ref_model, pid, qcov,
                     int(ref_model == cluster) if ref_model else 0, tier_of(pid, qcov)))
    con.execute("DELETE FROM ro_evidence")
    con.executemany("INSERT INTO ro_evidence VALUES (?,?,?,?,?,?)", rows)
    con.commit()

    counts = Counter(r[5] for r in rows)
    print("[ro_evidence] kademe dagilimi:")
    for t in TIER_ORDER:
        print(f"   {t:15s} {counts[t]:6d}  ({100*counts[t]/len(rows):.1f}%)")
    disagree = sum(1 for r in rows if r[1] and not r[4])
    print(f"   HMM kumesi != en yakin referans: {disagree} ({100*disagree/len(rows):.1f}%)")

    os.makedirs(args.out_dir, exist_ok=True)
    by_cluster = defaultdict(Counter)
    for r in rows:
        by_cluster[hmm_cluster[r[0]]][r[5]] += 1
    path = os.path.join(args.out_dir, "evidence_by_cluster.csv")
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["cluster", "n"] + TIER_ORDER + ["frac_characterized_or_close"])
        for cl, c in sorted(by_cluster.items(), key=lambda kv: -sum(kv[1].values())):
            n = sum(c.values())
            w.writerow([cl, n] + [c[t] for t in TIER_ORDER]
                       + [round((c["characterized"] + c["close_homolog"]) / n, 3)])
    print(f"[yazildi] {path}")
    con.close()


if __name__ == "__main__":
    main()
