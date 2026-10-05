"""
Elektron tasima zinciri (ETC) bilesenleri -- RO siniflandirmasinin biyokimyasal temeli.

Batie (1991) ve Kweon (2008) siniflandirmalari RO'lari alpha'nin dizisine gore
degil, elektronu NAD(P)H'den alpha'ya tasiyan ortaklara gore ayirir:
    reduktaz tipi     FNR-tipi (FAD/NAD baglayan, PF00970/PF00175/PF08022)
                      GR-tipi  (glutatyon reduktaz katlanmasi, PF07992/PF02852)
    ferredoksin tipi  Rieske-tipi [2Fe-2S] (PF00355/PF13806, <=250 aa)
                      bitki-tipi [2Fe-2S]  (PF00111 Fer2, Rieske yok, <=250 aa)
    beta alt birimi   PF00866 (alpha3beta3 vs alpha3)

Bu script her dogrulanmis RO icin +-10 kb icindeki (ve operondaki) ortaklari
neighbor_domain'den okur ve tanimlayici bir ETC profili cikarir. Kesin "Tip I-V"
etiketi vermek yerine gozlenen bilesen kombinasyonunu yazar; eksik anotasyon ya da
kontig sinirlari yuzunden "gorulmedi" ile "yok" ayni sey degildir.

Cikti: ro_etc tablosu + analysis_out/etc_by_cluster.csv
"""

import argparse
import csv
import os
import sqlite3
from collections import Counter, defaultdict

FNR = {"NAD_binding_1", "FAD_binding_6", "FAD_binding_8"}
GR = {"Pyr_redox_2", "Pyr_redox_dim"}
RIESKE = {"Rieske", "Rieske_2"}
PLANT = {"Fer2"}
MAX_FD_LEN = 250


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--db", default="roar.sqlite")
    ap.add_argument("--out-dir", default="analysis_out")
    args = ap.parse_args()

    con = sqlite3.connect(args.db)
    con.executescript("""
        CREATE TABLE IF NOT EXISTS ro_etc (
            candidate_id     TEXT PRIMARY KEY,
            reductase_type   TEXT,   -- FNR | GR | FNR+GR | none
            ferredoxin_type  TEXT,   -- rieske | plant | rieske+plant | none
            has_beta         INTEGER,
            in_operon        INTEGER, -- bu bilesenlerin en az biri operonda mi
            etc_profile      TEXT     -- okunabilir kombinasyon
        );""")

    # protein_key -> domain seti, uzunluk
    domains = defaultdict(set)
    for key, hmm in con.execute("SELECT protein_key, hmm FROM neighbor_domain"):
        domains[key].add(hmm)
    length = dict(con.execute("SELECT protein_key, length FROM neighbor_protein"))
    in_operon = defaultdict(set)
    for cid, key in con.execute(
            "SELECT candidate_id, protein_key FROM operon_gene WHERE protein_key IS NOT NULL"):
        in_operon[cid].add(key)

    rows = []
    for cid, in con.execute("SELECT candidate_id FROM ro WHERE is_confirmed=1"):
        red, fd, beta, op = set(), set(), False, False
        for (key,) in con.execute(
                "SELECT nucleotide_id||':'||start||'-'||end||':'||strand FROM neighbor WHERE candidate_id=?",
                (cid,)):
            d = domains.get(key)
            if not d:
                continue
            hit = False
            if d & FNR:
                red.add("FNR"); hit = True
            if d & GR:
                red.add("GR"); hit = True
            if length.get(key, 10 ** 6) <= MAX_FD_LEN:
                if d & RIESKE:
                    fd.add("rieske"); hit = True
                elif d & PLANT:
                    fd.add("plant"); hit = True
            if "Ring_hydroxyl_B" in d:
                beta = True; hit = True
            if hit and key in in_operon.get(cid, ()):
                op = True
        red_t = "+".join(sorted(red)) or "none"
        fd_t = "+".join(sorted(fd)) or "none"
        profile = (("α3β3" if beta else "α3") + " · Fd:" + fd_t + " · Red:" + red_t)
        rows.append((cid, red_t, fd_t, int(beta), int(op), profile))

    con.execute("DELETE FROM ro_etc")
    con.executemany("INSERT INTO ro_etc VALUES (?,?,?,?,?,?)", rows)
    con.commit()

    cluster = dict(con.execute("SELECT candidate_id, ro_cluster FROM ro WHERE is_confirmed=1"))
    by_cluster = defaultdict(Counter)
    for r in rows:
        by_cluster[cluster[r[0]]][r[5]] += 1
    os.makedirs(args.out_dir, exist_ok=True)
    path = os.path.join(args.out_dir, "etc_by_cluster.csv")
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["cluster", "n", "profile", "count", "fraction"])
        for cl, c in sorted(by_cluster.items(), key=lambda kv: -sum(kv[1].values())):
            n = sum(c.values())
            for p, k in c.most_common():
                w.writerow([cl, n, p, k, round(k / n, 3)])
    tot = Counter(r[5] for r in rows)
    print("[ro_etc] en sik ETC profilleri:")
    for p, k in tot.most_common(10):
        print(f"   {k:6d}  {p}")
    print(f"[yazildi] {path}")
    con.close()


if __name__ == "__main__":
    main()
