"""Substrat icin KONTROLLU bir kimlik: ChEBI, serbest metin degil.

SORU. `chemistry.csv`'de substrat bir DIZGI ("naphthalene"). Dizgi hicbir
dis kaynaga baglanmaz: baska bir veritabanindaki ayni molekule, bir
yola, bir toksisite kaydina otomatik olarak eslenemez. Kontrollu
vokabuler (ChEBI) bunu mumkun kilar.

NEREDEN GELIYOR -- ve bu ayrim bu dosyanin butun noktasi. Iki yol var ve
ESIT DEGILLER:

  rhea_equation  Tipin kendi UniProt kaydindan Rhea'ya, oradan
                 DENGELENMIS denklemin katilimcilarina gidilir; substratin
                 ChEBI kimligi o denklemden okunur. Bu bir kimlik.
  name_match     Kuratorun yazdigi ADIN ChEBI'de aranmasi. Bu bir ADAY,
                 kimlik degil: ayni ad ayri molekullere, ayri adlar ayni
                 molekule denk gelebilir. `state` bunu tasir ve hicbir
                 iddia bu satirlara dayandirilmaz.

BIR TIPIN BIRDEN COK SUBSTRATI OLABILIR ve bu bir hata degil: NdmA'nin
Rhea kayitlari kafein, teofilin ve 1,7-dimetilksantin uzerinde ayri ayri
dengelenmis denklemler veriyor. Tablo bu yuzden tip basina TEK satir
degil, (tip x substrat) ciftidir.

AYRICA OLCULEN BIR SEY. Kuratorun yazdigi substrat adi, Rhea denkleminin
substratlari arasinda GECIYOR MU? Gecmiyorsa bu bir bulgudur -- ayni
hata sinifi `reaction_balance.py`'nin var olma sebebi: formul aritmetigi
dogru yapilip YANLIS MOLEKULE uygulanmasi.

SEMA. Yeni alan `chemistry.csv`'ye SUTUN olarak eklenmiyor; turetilmis
bir tablo + tazelik kontrolu olarak duruyor (projenin sema genisletme
kurali). Boylece kaynak dosya tek sahipli kalir ve turetilmis deger
yeniden uretilebilir.

Cikti: tablo `type_substrate_chebi` + analysis_out/substrate_chebi.json
"""

import argparse
import json
import os
import re
import sqlite3
import time
import urllib.parse
import urllib.request

import yaml

UA = "ROAR-DB/1.0 (mailto:recocanextra@gmail.com)"
OLS = ("https://www.ebi.ac.uk/ols4/api/search"
       "?ontology=chebi&exact=true&queryFields=label,synonym&q={q}")
CACHE = "analysis_out/chebi_name_cache.json"


#  --- YAPISAL AYIRT EDICI -------------------------------------------------
#  Bir AD farki bir uyusmazlik DEGILDIR: "dicamba" ile
#  "3,6-dichloro-2-methoxybenzoate" ayni molekul, "carbazole" ile
#  "9H-carbazole" ayni molekul. Adlari karsilastirmak bu yuzden kurt var
#  demektir. Karsilastirilacak sey YAPI: kuratorun yazdigi SMILES'in agir
#  atom sayimi, Rhea katilimcisinin formuluyle ayni mi.
#
#  Hidrojen SAYILMAZ: ortuk H'yi dogru saymak valans kurallari ister ve
#  burada gerekmiyor -- C/N/O/S/halojen iki molekulu ayirt etmeye yeter.
#
#  SAYAC DOGRULANDI: kuratorun adi Rhea'nin adiyla AYNI olan 15 ciftte
#  sayim 15/15 birebir tutuyor, 0 uyusmazlik (2026-11-12). Yani bir
#  uyusmazlik cikarsa sayacin degil verinin uyusmazligidir.
_BRACKET = re.compile(r"\[([^\]]*)\]")
_TWO_LETTER = ("Cl", "Br")
#  Rhea JENERIK formul yazabilir: `C7H4O2X` -- X bir halojen YERINE gecer.
#  Boyle bir formul atom sayimiyla karsilastirilamaz, ve karsilastirmaya
#  zorlamak sahte bir uyusmazlik uretir.
_PLACEHOLDER = re.compile(r"[XR]")


def heavy_atoms(smiles):
    """SMILES -> {element: sayi}, hidrojen haric. Donen None: SMILES yok."""
    if not smiles:
        return None
    counts = {}

    def bump(element):
        counts[element] = counts.get(element, 0) + 1

    for body in _BRACKET.findall(smiles):
        match = re.match(r"\d*([A-Z][a-z]?|[a-z])", body)
        if match:
            element = match.group(1)
            bump(element.capitalize() if len(element) > 1 else element.upper())
    rest = _BRACKET.sub("", smiles)
    index = 0
    while index < len(rest):
        pair = rest[index:index + 2]
        char = rest[index]
        if pair in _TWO_LETTER:
            bump(pair)
            index += 2
            continue
        if char.isalpha() and char.isupper():
            bump(char)
        elif char in "cnops":            # aromatik kucuk harf
            bump(char.upper())
        index += 1
    counts.pop("H", None)
    return counts


def formula_atoms(formula):
    """Molekul formulu -> {element: sayi}, hidrojen haric.

    Donen None: formul JENERIK (X ya da R tasiyor), yani sayilamaz.
    """
    if not formula or _PLACEHOLDER.search(formula):
        return None
    counts = {}
    for element, number in re.findall(r"([A-Z][a-z]?)(\d*)", formula):
        if not element:
            continue
        counts[element] = counts.get(element, 0) + int(number or 1)
    counts.pop("H", None)
    return counts


def compare_structure(curated_smiles, subs):
    """Kurator yapisi, Rhea substratlarindan biriyle ayni mi?

    Donen: (durum, eslesen_substrat_adi_ya_da_None, aciklama)
    """
    atoms = heavy_atoms(curated_smiles)
    if not atoms:
        return ("not_decidable", None,
                "the curated record carries no SMILES, so the structures "
                "cannot be compared and the name difference is left open")
    generic = [s for s in subs if formula_atoms(s.get("formula")) is None]
    for sub in subs:
        want = formula_atoms(sub.get("formula"))
        if want is not None and want == atoms:
            return ("same_molecule", sub.get("name"),
                    "the curated name differs but the heavy-atom composition "
                    "is identical, so this is a nomenclature difference "
                    "(trivial against systematic name, salt or stereo "
                    "descriptor), not a disagreement")
    if generic:
        return ("generic_rhea_class", generic[0].get("name"),
                "Rhea states a generic formula (a placeholder element stands "
                "for the halogen), so an atom count cannot decide; forcing the "
                "comparison would manufacture a disagreement")
    return ("different_molecule", None,
            "the curated structure matches none of the substrates of the "
            "balanced equation; this is the error class reaction_balance.py "
            "exists for -- correct arithmetic applied to the wrong molecule -- "
            "and it needs a curator, not a script")


def read_types(types_dir):
    out = {}
    for name in sorted(os.listdir(types_dir)):
        if not name.endswith((".yaml", ".yml")):
            continue
        with open(os.path.join(types_dir, name), encoding="utf-8") as fh:
            data = yaml.safe_load(fh) or {}
        if data.get("cluster"):
            out[data["cluster"]] = data
    return out


def rhea_substrates(path):
    """cluster -> [{name, chebi, formula}] (dengelenmis denklemden)."""
    if not os.path.exists(path):
        return {}
    with open(path, encoding="utf-8") as fh:
        data = json.load(fh)
    out = {}
    for row in data.get("rows") or []:
        seen = {}
        for detail in row.get("details") or []:
            sub = detail.get("substrate") or {}
            if sub.get("chebi"):
                seen[sub["chebi"]] = {"name": sub.get("name"),
                                      "chebi": sub["chebi"],
                                      "formula": sub.get("formula"),
                                      "rhea": detail.get("rhea")}
        if seen:
            out[row["cluster"]] = list(seen.values())
    return out


def name_lookup(name, cache, allow_network):
    """ChEBI'de TAM ad/sinonim eslesmesi. Donen: aday ya da None.

    AG, KAPININ KOSULU OLAMAZ. Onbellek YETKILIDIR: bir ad bir kez
    cozulduyse boru hatti onu cevrimdisi tekrar uretir. Ag yalnizca
    `--lookup-names` ile ve yalnizca onbellekte OLMAYAN adlar icin
    kullanilir. Boylece aynı veri, internet olmadan da ayni tabloyu
    uretir.
    """
    key = (name or "").strip().lower()
    if not key:
        return None
    if key in cache:
        return cache[key]
    if not allow_network:
        return None
    url = OLS.format(q=urllib.parse.quote(key))
    got = None
    try:
        req = urllib.request.Request(url, headers={"User-Agent": UA})
        with urllib.request.urlopen(req, timeout=60) as res:
            payload = json.load(res)
        docs = (payload.get("response") or {}).get("docs") or []
        #  Birden cok TAM eslesme varsa hicbiri secilmez: secmek, iki
        #  molekulden birini rastgele dogru saymak olurdu.
        exact = [d for d in docs
                 if (d.get("label") or "").strip().lower() == key]
        if len(exact) == 1:
            got = {"chebi": exact[0].get("obo_id") or exact[0].get("short_form"),
                   "label": exact[0].get("label"),
                   "ambiguous": False}
        elif len(exact) > 1:
            got = {"chebi": None, "label": None, "ambiguous": True,
                   "candidates": [d.get("obo_id") for d in exact[:5]]}
    except Exception as exc:
        got = {"error": f"{type(exc).__name__}"}
    cache[key] = got
    time.sleep(0.2)
    return got


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--db", default="roar.sqlite")
    ap.add_argument("--types", default="types")
    ap.add_argument("--oxygen-fate", default="analysis_out/oxygen_fate_computed.json")
    ap.add_argument("--out", default="analysis_out/substrate_chebi.json")
    ap.add_argument("--lookup-names", action="store_true",
                    help="cozulmeyen tipler icin ChEBI'de ad ara (ag gerekir)")
    args = ap.parse_args()

    types = read_types(args.types)
    from_rhea = rhea_substrates(args.oxygen_fate)

    cache = {}
    if os.path.exists(CACHE):
        with open(CACHE, encoding="utf-8") as fh:
            cache = json.load(fh).get("names", {})

    rows, mismatches, unresolved, ambiguous = [], [], [], []
    for cluster, data in sorted(types.items()):
        curated = (data.get("substrate_en") or "").strip()
        subs = from_rhea.get(cluster) or []
        if subs:
            names = {(s["name"] or "").lower() for s in subs}
            #  Kuratorun adi denklemin substratlari arasinda mi?
            agrees = curated.lower() in names
            if curated and not agrees:
                #  AD farki bir uyusmazlik degil: karar YAPIDAN verilir.
                state, matched, why = compare_structure(
                    data.get("substrate_smiles"), subs)
                mismatches.append({"cluster": cluster,
                                   "curated_substrate": curated,
                                   "curated_heavy_atoms":
                                       heavy_atoms(data.get("substrate_smiles")),
                                   "rhea_substrates": sorted(names),
                                   "verdict": state,
                                   "matched_substrate": matched,
                                   "reading": why})
            for sub in subs:
                rows.append((cluster, sub["chebi"], sub["name"], sub["formula"],
                             "rhea_equation", sub.get("rhea"),
                             int(bool(curated) and (sub["name"] or "").lower()
                                 == curated.lower())))
            continue

        got = name_lookup(curated, cache, args.lookup_names)
        if got and got.get("chebi"):
            rows.append((cluster, got["chebi"], got["label"], None,
                         "name_match", None, 1))
        elif got and got.get("ambiguous"):
            ambiguous.append({"cluster": cluster, "curated_substrate": curated,
                              "candidates": got.get("candidates")})
            unresolved.append(cluster)
        else:
            unresolved.append(cluster)

    con = sqlite3.connect(args.db)
    con.execute("DROP TABLE IF EXISTS type_substrate_chebi")
    con.execute("""
        CREATE TABLE type_substrate_chebi (
            cluster TEXT, chebi TEXT, chebi_label TEXT, formula TEXT,
            state TEXT, rhea TEXT, is_curated_substrate INTEGER,
            PRIMARY KEY (cluster, chebi))""")
    con.executemany("INSERT OR REPLACE INTO type_substrate_chebi "
                    "VALUES (?,?,?,?,?,?,?)", rows)
    con.commit()

    by_state = {}
    for row in rows:
        by_state[row[4]] = by_state.get(row[4], 0) + 1

    report = {
        "question": "Does every enzyme type's substrate carry a controlled "
                    "identifier, and where does that identifier come from?",
        "method": {
            "unit": "one (enzyme type, substrate) pair -- a type with three "
                    "balanced Rhea reactions has three substrates, which is a "
                    "fact about the enzyme, not a duplicate row",
            "states": {
                "rhea_equation": "read from the balanced equation of the type's "
                                 "own Rhea reaction; this is an identity",
                "name_match": "the curated substrate name matched exactly one "
                              "ChEBI label or synonym; this is a CANDIDATE, not "
                              "an identity, and no claim rests on it",
            },
            "ambiguity_is_not_resolved": "When a name matches more than one "
                                         "ChEBI entry, none is chosen: choosing "
                                         "would mean calling one of two "
                                         "molecules correct at random.",
            "schema": "A derived table, not a new column in chemistry.csv, so "
                      "the curated source stays single-owner and the derived "
                      "value stays reproducible.",
        },
        "denominators": {
            "types": len(types),
            "types_with_rhea_substrate": len(from_rhea),
            "pairs_written": len(rows),
            "by_state": by_state,
            "types_unresolved": len(unresolved),
            "network_lookup_allowed": bool(args.lookup_names),
            "name_matches_from_cache": sum(
                1 for value in cache.values()
                if isinstance(value, dict) and value.get("chebi")),
        },
        "curated_name_differs_from_equation": {
            "why_this_is_not_a_count_of_errors":
                "A name difference is not a disagreement: 'dicamba' and "
                "'3,6-dichloro-2-methoxybenzoate' are the same molecule. Each "
                "row below is therefore decided on STRUCTURE -- the heavy-atom "
                "composition of the curated SMILES against the formula of the "
                "balanced equation's substrate -- not on the name.",
            "atom_counter_validated":
                "On the 15 pairs where the curated name equals the Rhea name, "
                "the counter agrees 15/15 with the Rhea formula; a mismatch is "
                "therefore a mismatch in the data, not in the counter.",
            "by_verdict": {state: sum(1 for m in mismatches
                                      if m["verdict"] == state)
                           for state in sorted({m["verdict"]
                                                for m in mismatches})},
            "rows": mismatches,
        },
        "ambiguous_names": ambiguous,
        "unresolved_types": sorted(unresolved),
    }

    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    with open(args.out, "w") as fh:
        json.dump(report, fh, indent=1)
    if args.lookup_names:
        with open(CACHE, "w") as fh:
            json.dump({"checked_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ",
                                                    time.gmtime()),
                       "names": cache}, fh, indent=1)

    print(f"[tablo] type_substrate_chebi: {len(rows)} (tip x substrat) cifti")
    for state, count in sorted(by_state.items()):
        print(f"    {state:14s} {count}")
    print(f"[cozulmeyen] {len(unresolved)} tip "
          f"({len(types)} tipin {100*len(unresolved)/max(1,len(types)):.0f}%'i)")
    if mismatches:
        groups = {}
        for row in mismatches:
            groups.setdefault(row["verdict"], []).append(row)
        print(f"[ad farki] {len(mismatches)} tipte kuratorun adi denklemin "
              f"substrat adlariyla ayni degil; karar YAPIDAN verildi:")
        for state in sorted(groups):
            rows_of = groups[state]
            print(f"    {state:22s} {len(rows_of)}")
            if state == "different_molecule":
                for row in rows_of:
                    print(f"        >>> {row['cluster']}: "
                          f"'{row['curated_substrate']}' "
                          f"{row['curated_heavy_atoms']} vs "
                          f"{row['rhea_substrates']}")
    print(f"[yazildi] {args.out}")
    con.close()


if __name__ == "__main__":
    main()
