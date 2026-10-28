"""Kopru karboksilat kolonu GERCEKTEN kopruyu mu gosteriyor? Yapiya karsi olc.

NEDEN. Katalitik karboksilat kolonu (471) 27 kuratorlu referansin kristal
yapisina karsi sinandi ve cagirdigi kalinti 26/26'sinda demirden 10,9-14,7 A
uzakta cikti: kolon ligandi OKUMUYORDU. Kopru kolonu (286) ayni riski tasiyor
ve BUGUNE KADAR yapilara karsi HIC kontrol edilmedi -- yine de on sayfa onun
uzerine "Asp halka-hidroksilleyenlerde, Glu kuaterner amin dalinda, yani tek
bir kolon substrat kimyasini izliyor" diye bir iddia kuruyor.

KOPRUNUN ISI farkli: bir alt birimin Rieske [2Fe-2S] merkezini KOMSU alt
birimin mononukleer demirine baglamak. Dolayisiyla test de farkli: kalinti iki
merkeze de yakin olmali VE iki ZINCIRI birbirine baglamali. Naftalin
dioksijenazda bu Asp205'tir ve karboksilati bir yanda Rieske ligandi His104'e,
obur yanda komsu alt birimin Fe ligandi His208'e hidrojen bagi yapar.

Bu betik SINIFLANDIRMAYI DEGISTIRMEZ, olcer ve raporlar.
"""
import json
import math
import os
import sys
import time
import urllib.request

UA = "ROAR-DB/1.0 (mailto:recocanextra@gmail.com)"
OUT = "analysis_out/bridging_geometry_audit.json"
CARBOX_ATOMS = {"OD1", "OD2", "OE1", "OE2", "SG"}


def fetch(pdb, cache="/media/lin-bio/back2/roar_pdb"):
    os.makedirs(cache, exist_ok=True)
    path = os.path.join(cache, f"{pdb.lower()}.cif")
    if not os.path.exists(path) or os.path.getsize(path) < 2000:
        url = f"https://files.rcsb.org/download/{pdb.upper()}.cif"
        req = urllib.request.Request(url, headers={"User-Agent": UA})
        data = urllib.request.urlopen(req, timeout=120).read()
        with open(path, "wb") as fh:
            fh.write(data)
    return path


def atoms(path):
    """mmCIF atom_site satirlari -> (chain, comp, seq, atom, x, y, z)."""
    out = []
    cols, in_loop = {}, False
    with open(path) as fh:
        for line in fh:
            if line.startswith("_atom_site."):
                cols[line.strip().split(".")[1]] = len(cols)
                in_loop = True
                continue
            if in_loop and (line.startswith("ATOM") or
                            line.startswith("HETATM")):
                p = line.split()
                if len(p) < len(cols):
                    continue
                try:
                    out.append((
                        p[cols.get("auth_asym_id", cols.get("label_asym_id"))],
                        p[cols["label_comp_id"]],
                        p[cols.get("auth_seq_id", cols.get("label_seq_id"))],
                        p[cols["label_atom_id"]],
                        float(p[cols["Cartn_x"]]), float(p[cols["Cartn_y"]]),
                        float(p[cols["Cartn_z"]])))
                except (KeyError, ValueError, IndexError):
                    continue
            elif in_loop and line.startswith("#"):
                if out:
                    break
    return out


def dist(a, b):
    return math.sqrt(sum((a[i] - b[i]) ** 2 for i in (4, 5, 6)))


def main():
    targets = json.load(open(sys.argv[1] if len(sys.argv) > 1 else
                             "/tmp/bridge_targets.json"))
    rows = []
    for cluster, t in sorted(targets.items()):
        if not (t.get("pdb") and t.get("resnum")):
            continue
        try:
            at = atoms(fetch(t["pdb"]))
        except Exception as exc:                        # noqa: BLE001
            rows.append({"cluster": cluster, "error": str(exc)[:90]})
            print(f"[{cluster}] HATA {exc}", file=sys.stderr)
            continue
        # Rieske kumesi: FES/SF4 hetero gruplari; mononukleer demir: tek FE
        rieske = [a for a in at if a[1] in ("FES", "SF4")]
        mono = [a for a in at if a[1] == "FE" and a[3] == "FE"]
        mine = [a for a in at if a[2] == str(t["resnum"])
                and a[3] in CARBOX_ATOMS
                and (t.get("chain") is None or a[0] == t["chain"])]
        if not mine:
            mine = [a for a in at if a[2] == str(t["resnum"])
                    and a[3] in CARBOX_ATOMS]
        rec = {"cluster": cluster, "pdb": t["pdb"], "chain": t.get("chain"),
               "residue": f"{t['residue']}{t['resnum']}",
               "n_rieske_clusters": len(rieske) // 4 if rieske else 0,
               "n_mono_fe": len(mono)}
        if not mine:
            rec["verdict"] = "the called residue has no carboxylate atom in " \
                             "the structure"
            rows.append(rec)
            print(f"[{cluster}] kalinti yapida bulunamadi", file=sys.stderr)
            continue
        def closest(pool):
            best = None
            for m in mine:
                for p in pool:
                    d = dist(m, p)
                    if best is None or d < best[0]:
                        best = (d, p[0], f"{p[1]}{p[2]}")
            return best
        cr, cm = closest(rieske), closest(mono)
        rec["to_rieske"] = ({"distance": round(cr[0], 2), "chain": cr[1],
                             "partner": cr[2]} if cr else None)
        rec["to_mononuclear_fe"] = ({"distance": round(cm[0], 2),
                                     "chain": cm[1], "partner": cm[2]}
                                    if cm else None)
        # YON ONEMLI, ve ilk yazimda TERS kodlanmisti. 1NDO'da Asp205,
        # KOMSU alt birimin (C zinciri) Rieske kumesinden 6,70 A ve KENDI
        # zincirinin (A) demirinden 6,95 A uzakta. Yani kopru, demiri tasiyan
        # alt birimde oturur ve KOMSUNUN Rieske merkezine uzanir. Ters test
        # 12/12'de "hayir" demisti -- bir bulgu degil, benim hatam.
        rec["spans_two_chains"] = bool(
            cr and cm and t.get("chain") and
            (cm[1] == t["chain"]) and (cr[1] != t["chain"]))
        # Yapida tek alt birim varsa komsunun Rieske'si DOSYADA YOKTUR ve
        # Rieske uzakligi anlamsiz buyur (kendi merkezine olan mesafe).
        rec["partner_subunit_in_file"] = rec["n_rieske_clusters"] > 1
        rows.append(rec)
        print(f"[{cluster}] {rec['residue']} Rieske "
              f"{rec['to_rieske'] and rec['to_rieske']['distance']} A, Fe "
              f"{rec['to_mononuclear_fe'] and rec['to_mononuclear_fe']['distance']} A, "
              f"iki zincir={rec['spans_two_chains']}", file=sys.stderr)
        time.sleep(0.2)
    os.makedirs("analysis_out", exist_ok=True)
    json.dump({"_about": "Kopru karboksilat kolonunun (286) yapiya karsi "
                         "geometrik denetimi. Siniflandirma degismez.",
               "_generated": time.strftime("%Y-%m-%d %H:%M UTC",
                                           time.gmtime()),
               "rows": rows}, open(OUT, "w"), indent=1, ensure_ascii=False)
    good = [r for r in rows if r.get("to_rieske") and r.get("to_mononuclear_fe")]
    print(f"\nolculen: {len(good)} / {len(rows)}")
    if good:
        dr = sorted(r["to_rieske"]["distance"] for r in good)
        dm = sorted(r["to_mononuclear_fe"]["distance"] for r in good)
        print(f"Rieske'ye uzaklik  : {dr[0]}-{dr[-1]} A, ortanca "
              f"{dr[len(dr)//2]}")
        print(f"mononukleer Fe'ye  : {dm[0]}-{dm[-1]} A, ortanca "
              f"{dm[len(dm)//2]}")
        multi = [r for r in good if r.get("partner_subunit_in_file")]
        print(f"iki zinciri baglayan: "
              f"{sum(1 for r in good if r['spans_two_chains'])}/{len(good)}")
        if multi:
            dr2 = sorted(r["to_rieske"]["distance"] for r in multi)
            dm2 = sorted(r["to_mononuclear_fe"]["distance"] for r in multi)
            print(f"\nKOMSU ALT BIRIM DOSYADA OLAN {len(multi)} yapi:")
            print(f"  komsunun Rieske'sine : {dr2[0]}-{dr2[-1]} A, ortanca "
                  f"{dr2[len(dr2)//2]}")
            print(f"  kendi demirine       : {dm2[0]}-{dm2[-1]} A, ortanca "
                  f"{dm2[len(dm2)//2]}")
            print(f"  iki zinciri baglayan : "
                  f"{sum(1 for r in multi if r['spans_two_chains'])}/{len(multi)}")
    print("->", OUT)


if __name__ == "__main__":
    main()
