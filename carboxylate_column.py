# -*- coding: utf-8 -*-
"""Mononukleer demirin karboksilat KOLONU -- yapilardan olculur.

NEDEN BU DOSYA VAR
  Atlas sayfasinda duran bir paragraf ("the carboxylate column") uc dosyaya
  dayaniyordu: carboxylate_column_audit.json, ..._raw.json ve
  carboxylate_clade_shift.json. Ucunu de HICBIR adim uretmiyordu -- 2026-10-28
  tarihli tek seferlik bir oturumun ciktilariydi. Yani yayimlanan sayilar
  yeniden uretilemiyor, dogrulanamiyor ve referans seti degistiginde sessizce
  bayatliyordu. Kullanici (2026-11-15) "yapmak gerekir" dedi.

OLCUMU YENIDEN YAPMIYOR; ZATEN OLCULMUS OLANI TUREVLIYOR
  Bu betigin en onemli karari, yapilari yeniden ayristirMAmak. active_site.py
  her referans yapisi icin sunlari ZATEN olcuyor ve yayinliyor:
    - katalitik mononukleer demiri ve ona koordine olan kalintilari
      (`catalytic_iron.coordinating_contacts.carboxylate_oxygen`)
    - motif modelinin her metal-ligand kolonunun yapida hangi kalintiya
      dustugunu, o kalintinin dogru metali BAGLAYIP baglamadigini ve kolon
      kaymasini (`column_mapping_verification.checks`)
    - her havuz kalintisinin hizalama kolonunu (`pocket[r].residues`)
  Ayni olcumu ikinci bir kod yolundan yapmak, iki sayinin ayrismasi demekti;
  bu projede bir olcumun iki uygulamasi yasak. Bu yuzden girdi yapi dosyalari
  degil, active_site.json.

NE OLCULUYOR
  Karboksilat kolonu (ro_motif.CARBOXYLATE_COLUMN) yapida gercekten demiri
  koordine eden Asp/Glu'ya mi dusuyor? Dusmuyorsa kac kolon kayiyor? Kayan
  yapilarda gercek karboksilat hangi kolonda duruyor?

NE OLCULMUYOR
  Yapisi olmayan referanslar hakkinda hicbir sey. Yapisi olup active_site
  tarafindan cozulemeyenler de (indirilemeyen ya da demiri olmayan yapilar)
  AYRI sayilir ve sebebi yazilir: "olculemedi" ile "kolonda yok" ayni sey
  degil, ve bu projede tam bu karistirma daha once yayimlanmis bir sayiyi
  yanlis yapti.
"""
import argparse
import json
import os
import sys


def load(path):
    with open(path) as handle:
        return json.load(handle)


def carboxylate_check(record, column):
    """Karboksilat kolonunun bu yapidaki denetimi -- yoksa None."""
    for check in ((record.get("column_mapping_verification") or {})
                  .get("checks") or []):
        if check.get("column") == column:
            return check
    return None


def pocket_column_of(record, residue):
    """Bir kalintinin hizalama kolonu: en kucuk yaricapli havuzda arar."""
    pocket = record.get("pocket") or {}
    for key in sorted(pocket, key=lambda k: float(k)):
        for row in (pocket[key].get("residues") or []):
            if row.get("pdb_residue") == residue:
                return row.get("alignment_column")
    return None


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--active-site", default="analysis_out/active_site.json")
    ap.add_argument("--out", default="analysis_out/carboxylate_column.json")
    ap.add_argument("--root", default=".")
    args = ap.parse_args()

    sys.path.insert(0, args.root)
    import ro_motif                                      # noqa: PLC0415
    #  Kolon numarasi ro_motif'ten OKUNUR, buraya yazilmaz: kolon bir gun
    #  kayarsa iki yerde iki sayi olmasin.
    column = next(col for col, res, _label in ro_motif.CATALYTIC_SITES
                  if res == "DE")

    data = load(os.path.join(args.root, args.active_site))
    rows, unresolved = [], []

    #  `undetermined` LISTESI DE OKUNUR. active_site cozemedigi yapiyi ayri
    #  bir listeye koyuyor; yalnizca `structures`'a bakmak, o yapiyi
    #  olcumden sessizce dusururdu ve "25 yapidan 24'u" ile "24 yapi"
    #  arasindaki fark kaybolurdu. Sebebi de yazilir: HcaE/8K0A'da
    #  katalitik kriterleri saglayan hicbir demir yok, yani kolon hakkinda
    #  bir sey soylemiyor -- bu bir YOKLUK degil, olculemeyen bir durum.
    for record in (data.get("undetermined") or []):
        unresolved.append({"type": record.get("type"),
                           "pdb_id": record.get("pdb_id"),
                           "why": (record.get("problems") or
                                   ["active_site could not resolve this structure"])[0]})

    for record in (data.get("structures") or []):
        check = carboxylate_check(record, column)
        iron = record.get("catalytic_iron") or {}
        contacts = (iron.get("coordinating_contacts") or {}).get(
            "carboxylate_oxygen") or []
        if not iron:
            unresolved.append({"type": record.get("type"),
                               "pdb_id": record.get("pdb_id"),
                               "why": "no catalytic mononuclear iron resolved "
                                      "in this structure"})
            continue
        exact = (check or {}).get("exact_column") or {}
        windowed = None
        for shift in ((record.get("column_mapping_verification") or {})
                      .get("column_shifts") or []):
            if shift.get("model_column") == column:
                windowed = shift
        #  GERCEK karboksilat: demire koordine olan Asp/Glu. Birden fazla
        #  oksijen ayni kalintidan gelebilir (OD1 ve OD2), kalinti adi tekil.
        real = sorted({c[0].split("/")[-1].split(":")[0] for c in contacts})
        real_columns = {r: pocket_column_of(record, r) for r in real}
        rows.append({
            "type": record.get("type"),
            "pdb_id": record.get("pdb_id"),
            "model_column": column,
            "residue_at_model_column": exact.get("structure_residue"),
            "residue_is_asp_or_glu": exact.get("residue_matches_expectation"),
            "ligates_the_catalytic_iron": exact.get("ligates_the_expected_metal"),
            "offset_from_model_column": (check or {}).get(
                "offset_from_model_column"),
            "observed_at_column": (check or {}).get("observed_at_column"),
            "shift": windowed,
            "iron_coordinating_carboxylates": real,
            "columns_of_those_carboxylates": real_columns,
            "closest_carboxylate_oxygen_A": (min(c[1] for c in contacts)
                                             if contacts else None),
        })

    at_column = [r for r in rows if r.get("observed_at_column") == column]
    within_one = [r for r in rows
                  if isinstance(r.get("observed_at_column"), int)
                  and abs(r["observed_at_column"] - column) <= 1]
    ligating = [r for r in rows if r.get("ligates_the_catalytic_iron")]
    off_window = [r for r in rows if r not in within_one]
    other_columns = sorted({c for r in rows
                            for c in r["columns_of_those_carboxylates"].values()
                            if isinstance(c, int) and c != column})
    dists = [r["closest_carboxylate_oxygen_A"] for r in rows
             if isinstance(r.get("closest_carboxylate_oxygen_A"), (int, float))]
    no_carboxylate = [{"type": r["type"], "pdb_id": r["pdb_id"]}
                      for r in rows if not r["iron_coordinating_carboxylates"]]

    #  KIMLIK: demiri koordine eden karboksilat Asp mi Glu mu? Sayfa bunu
    #  "aspartate in N of N structures" diye yaziyor ve sayi OLCULMELI;
    #  kalinti adinin ilk uc harfi yapidan geliyor, tahmin degil.
    identity = {}
    for row in rows:
        for residue in row["iron_coordinating_carboxylates"]:
            identity[residue[:3]] = identity.get(residue[:3], 0) + 1

    out = {
        "module": "carboxylate_column.py",
        "what_was_measured":
            "Whether the motif model's carboxylate column lands, in each "
            "curated reference structure, on the Asp or Glu that actually "
            "coordinates the mononuclear iron. Derived from active_site.json, "
            "which is where the structures are read; nothing here parses a "
            "structure a second time.",
        "what_this_is_not":
            "It says nothing about references without a structure, and a "
            "structure that could not be resolved is counted separately with "
            "its reason. A structure that was never downloaded is not "
            "evidence that the column is wrong there.",
        "model_column": column,
        "structures_examined": len(rows),
        "structures_unresolved": unresolved,
        "at_the_model_column": len(at_column),
        "within_one_column": len(within_one),
        "outside_the_window": [{"type": r["type"], "pdb_id": r["pdb_id"],
                                "observed_at_column": r["observed_at_column"],
                                "columns_of_those_carboxylates":
                                    r["columns_of_those_carboxylates"]}
                               for r in off_window],
        "column_lands_on_an_iron_ligand": len(ligating),
        "other_columns_seen_for_an_iron_coordinating_carboxylate": other_columns,
        "closest_oxygen_distance_A": ({"n": len(dists), "min": min(dists),
                                       "max": max(dists)} if dists else None),
        "structures_with_no_carboxylate_on_the_iron": no_carboxylate,
        "identity_of_the_iron_coordinating_carboxylate": identity,
        "per_reference": rows,
    }
    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    with open(args.out, "w") as handle:
        json.dump(out, handle, indent=1)
    print("[yazildi] %s" % args.out)
    print("  kolon %d | %d yapi incelendi, %d'si kolonda, %d'si +-1 icinde, "
          "%d'si demiri baglayan bir D/E'ye dusuyor"
          % (column, len(rows), len(at_column), len(within_one), len(ligating)))
    if unresolved:
        print("  cozulemeyen yapi: %d (sebepleri dosyada)" % len(unresolved))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
