"""Varyant kimliklerinin KARARLI oldugunu dogrula.

NEDEN BU BETIK VAR. Varyant kimligi artik uye dizilerinden turetiliyor
(recursive_homogenize.variant_id). Tasarimin isledigini soyleyecek tek sinama
sudur: ayni uye kumesi her kosuda AYNI kimligi almali, ve bir referans
eklendiginde yalnizca uyeligi GERCEKTEN degisen kumelerin kimligi degismeli.

Eski kimlik f"{cluster}#{leaf_index}" idi ve bu sinanamiyordu: tipin uyeligi
degisince butun kardes yapraklar yeniden numaralaniyordu ve kimse fark
etmiyordu. Bu depoda ayni sinif icin uc sinama var (md5 parca anahtari, SVG
portu, iki render yolunun metni); bu dorduncusu.

Iki kip:
    --record   su anki kimlik -> uye-ozeti eslemesini kaydet (temel)
    (varsayilan) kaydedilmis temele karsi dogrula

Cikis kodu 1: ayni uye kumesi farkli bir kimlik almis, ya da kimlik ayni
kalirken uye kumesi degismis (ikisi de kimlik uretiminin bozuldugu anlamina
gelir).
"""
import argparse
import hashlib
import json
import os
import sqlite3
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT_STATE = os.path.join(HERE, "analysis_out", "variant_ids.json")


def member_digest(con, variant):
    """Bir kumenin uye kumesinin ozeti -- kimlikten BAGIMSIZ hesaplanir."""
    rows = con.execute(
        "SELECT r.sequence FROM ro_leaf l JOIN ro r USING(candidate_id) "
        "WHERE l.variant_id=? AND r.sequence IS NOT NULL", (variant,)).fetchall()
    h = hashlib.sha256()
    for d in sorted(hashlib.sha256((r[0] or "").encode()).hexdigest() for r in rows):
        h.update(d.encode())
    return h.hexdigest()[:16], len(rows)


def snapshot(db):
    con = sqlite3.connect(f"file:{db}?mode=ro", uri=True)
    try:
        variants = [r[0] for r in con.execute(
            "SELECT DISTINCT variant_id FROM leaf WHERE variant_id IS NOT NULL")]
        out = {}
        for v in variants:
            digest, n = member_digest(con, v)
            out[v] = {"members": digest, "n": n}
        return out
    finally:
        con.close()


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--db", default=os.path.join(HERE, "roar.sqlite"))
    ap.add_argument("--state", default=DEFAULT_STATE)
    ap.add_argument("--record", action="store_true", help="temeli yaz")
    a = ap.parse_args()

    now = snapshot(a.db)
    print(f"[kume] {len(now)}")

    if a.record or not os.path.exists(a.state):
        os.makedirs(os.path.dirname(a.state), exist_ok=True)
        with open(a.state, "w") as fh:
            json.dump(now, fh, indent=1, sort_keys=True)
        print(f"[kaydedildi] {a.state}" if a.record
              else f"[temel kurulmadi, yazildi] {a.state}")
        return 0

    with open(a.state) as fh:
        before = json.load(fh)

    # 1. Ayni kimlik, FARKLI uye kumesi -> kimlik icerigi temsil etmiyor.
    drifted = [v for v in set(before) & set(now)
               if before[v]["members"] != now[v]["members"]]
    # 2. Ayni uye kumesi, FARKLI kimlik -> kimlik kararsiz.
    by_members_before = {}
    for v, d in before.items():
        by_members_before.setdefault(d["members"], set()).add(v)
    renamed = []
    for v, d in now.items():
        owners = by_members_before.get(d["members"])
        if owners and v not in owners:
            renamed.append((v, sorted(owners)[0]))

    gone = sorted(set(before) - set(now))
    fresh = sorted(set(now) - set(before))

    print(f"  korunan : {len(set(before) & set(now))}")
    print(f"  yeni    : {len(fresh)}")
    print(f"  kaybolan: {len(gone)}")
    print()
    bad = 0
    if drifted:
        bad = 1
        print(f"HATA: {len(drifted)} kimlikte uye kumesi degismis ama kimlik AYNI "
              "kalmis -- kimlik icerigi temsil etmiyor:")
        for v in drifted[:5]:
            print(f"   {v}: uye {before[v]['n']} -> {now[v]['n']}")
    if renamed:
        bad = 1
        print(f"HATA: {len(renamed)} kumede uye kumesi AYNI ama kimlik degismis "
              "-- kimlik kararsiz:")
        for new, old in renamed[:5]:
            print(f"   {old} -> {new}")
    if not bad:
        print("KARARLI: ayni uye kumesi ayni kimligi aliyor, ve kimlikler "
              "yalnizca uyelik gercekten degistiginde degisiyor.")
    return bad


if __name__ == "__main__":
    sys.exit(main())
