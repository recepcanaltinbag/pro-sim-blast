"""Oksijen kaderini RESMI DENGELENMIS DENKLEMDEN HESAPLA -- okuyarak degil.

NEDEN BU MODUL VAR. Oksijen kaderi 2026-10-28'e kadar kuratorun (ve benim)
OKUMASINA dayaniyordu ve kolin monooksijenazda IKI KEZ yanlis cikti. Hatanin
sinifi her seferinde ayni: formul aritmetigi dogru yapildi ama YANLIS
MOLEKULUN formuluyle, cunku kaydimizdaki urun adi enzimin biraktigi molekul
degildi. Insan okumasi bu hatayi uretmeye devam eder; atom sayimi uretmez.

TEMEL. Uc resmi kaynak zincirlenir:
  1. UniProt  -> tipin EC numarasi ve Rhea reaksiyon kimligi
  2. Rhea     -> DENGELENMIS denklem ve katilimcilarin ChEBI kimlikleri
                 (ChEBI listesi denklem SIRASINDA gelir, yani yanlara bolunur)
  3. ChEBI    -> her katilimcinin molekul formulu (OLS4 uzerinden)
Sonra oksijen atomlari sayilir. Kurator adlandirmasi bu zincire HIC girmez.

SINIR. Zincir yalnizca Rhea kimligi olan tipler icin kurulabilir. Kimligi
olmayan tipin kimyasi hala yalnizca kurasyona dayanir ve bu GORUNUR olmak
zorunda -- `chemistry_evidence` alani bunu kaydeder.
"""
import json
import os
import re
import time
import urllib.parse
import urllib.request

UA = "ROAR-DB/1.0 (mailto:recocanextra@gmail.com)"
CACHE = "analysis_out/chebi_formula_cache.json"

# Reaksiyonun her iki yaninda da gorev yapan, ISKELET SAYILMAYACAK
# katilimcilar. Oksijen muhasebesi bunlari ayri tutar: O2 girdi, H2O ve
# yukseltgenmis tasiyicilar cikti tarafinda beklenir.
COFACTOR_CHEBI = {
    "CHEBI:15379",  # O2
    "CHEBI:15377",  # H2O
    "CHEBI:15378",  # H(+)
    "CHEBI:16474", "CHEBI:18009",  # NADPH / NADP(+)
    "CHEBI:57540", "CHEBI:57945",  # NAD(+) / NADH
    "CHEBI:33737", "CHEBI:33738",  # oxidized / reduced [2Fe-2S] ferredoxin
    "CHEBI:29033", "CHEBI:29034",  # Fe(2+) / Fe(3+)
    "CHEBI:57783", "CHEBI:58349",  # NADPH / NADP(+) (alt forms)
}


def _get(url, tries=3):
    for i in range(tries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": UA})
            return urllib.request.urlopen(req, timeout=60).read()
        except Exception:
            if i == tries - 1:
                raise
            time.sleep(2 * (i + 1))


def _load_cache():
    if os.path.exists(CACHE):
        try:
            return json.load(open(CACHE))
        except Exception:                               # noqa: BLE001
            return {}
    return {}


_CACHE = _load_cache()


def save_cache():
    os.makedirs(os.path.dirname(CACHE), exist_ok=True)
    json.dump(_CACHE, open(CACHE, "w"), indent=1, sort_keys=True)


def chebi_formula(chebi):
    """ChEBI kimliginin molekul formulu. OLS4'ten okunur ve ONBELLEKLENIR."""
    if chebi in _CACHE:
        return _CACHE[chebi]
    url = ("https://www.ebi.ac.uk/ols4/api/ontologies/chebi/terms?obo_id="
           + urllib.parse.quote(chebi))
    val = None
    try:
        d = json.loads(_get(url))
        terms = (d.get("_embedded") or {}).get("terms") or []
        if terms:
            ann = terms[0].get("annotation") or {}
            for key in ("generalized_empirical_formula", "formula"):
                if ann.get(key):
                    val = ann[key][0]
                    break
            _CACHE[chebi] = {"formula": val, "label": terms[0].get("label")}
    except Exception as exc:                            # noqa: BLE001
        _CACHE[chebi] = {"formula": None, "error": str(exc)[:80]}
    return _CACHE.get(chebi)


ATOM = re.compile(r"([A-Z][a-z]?)(\d*)")


def atoms(formula):
    """Formulden atom sayimi. Yuk isaretleri ve parantezler yok sayilir."""
    if not formula:
        return {}
    out = {}
    for el, n in ATOM.findall(re.sub(r"[^A-Za-z0-9]", "", formula)):
        if not el:
            continue
        out[el] = out.get(el, 0) + (int(n) if n else 1)
    return out


def rhea_reaction(rhea_id):
    """Dengelenmis denklem ve katilimcilar, yanlara bolunmus halde."""
    url = ("https://www.rhea-db.org/rhea?query=" + urllib.parse.quote(rhea_id)
           + "&columns=rhea-id,equation,chebi-id&format=tsv")
    txt = _get(url).decode()
    lines = [ln for ln in txt.splitlines() if ln.strip()]
    if len(lines) < 2:
        return None
    _rid, equation, chebis = lines[1].split("\t")[:3]
    ids = [c.strip() for c in chebis.split(";") if c.strip()]
    if "=" not in equation:
        return None
    left_txt, right_txt = equation.split("=", 1)
    split = lambda s: [p.strip() for p in s.split(" + ") if p.strip()]
    lnames, rnames = split(left_txt), split(right_txt)
    # ChEBI listesi denklem SIRASINDA gelir; sayilar uyusmuyorsa bolme
    # GUVENILMEZ ve None donulur -- yanlis bolme sessiz bir hata olurdu.
    if len(ids) != len(lnames) + len(rnames):
        return {"equation": equation, "chebi": ids, "left": None,
                "right": None,
                "note": "ChEBI count does not match the participant count, "
                        "so the sides cannot be assigned safely"}
    pair = lambda names, offset: [
        {"name": names[i], "chebi": ids[offset + i]} for i in range(len(names))]
    return {"equation": equation, "left": pair(lnames, 0),
            "right": pair(rnames, len(lnames))}


def _skeleton(parts):
    """Iskelet katilimcilar: kofaktor olmayanlar, karbon sayisina gore buyuk
    once."""
    out = []
    for p in parts:
        if p["chebi"] in COFACTOR_CHEBI:
            continue
        f = (chebi_formula(p["chebi"]) or {}).get("formula")
        a = atoms(f)
        out.append({**p, "formula": f, "C": a.get("C", 0), "O": a.get("O", 0)})
    out.sort(key=lambda x: -x["C"])
    return out


def derive_fate(rhea_id):
    """Resmi denklemden oksijen kaderini HESAPLA.

    product       : ana urun (en buyuk karbon iskeleti) oksijen KAZANMIS
    leaving_group : ana urun kazanmamis ama baska bir urun oksijen tasiyor
    oxidant_only  : urun tarafinda oksijen kazanci YOK
    ambiguous     : iskelet parcalanmis ve substrat oksijen tasiyor -- sayim
                    atom izleyemez, deger HESAPLANMAZ
    """
    r = rhea_reaction(rhea_id)
    if not r or not r.get("left"):
        return {"rhea": rhea_id, "fate": None,
                "why": (r or {}).get("note", "reaction could not be read")}
    L, R = _skeleton(r["left"]), _skeleton(r["right"])
    if not L or not R or any(x["formula"] is None for x in L + R):
        return {"rhea": rhea_id, "fate": None, "equation": r["equation"],
                "why": "a participant has no formula in ChEBI"}
    sub, main = L[0], R[0]
    gain_main = main["O"] - sub["O"]
    others = sum(x["O"] for x in R[1:])
    total_gain = sum(x["O"] for x in R) - sum(x["O"] for x in L)

    # SAYIM ATOM IZLEMEZ -- ve bunu bilmek zorunda.
    #
    # Iskelet PARCALANDIYSA (ana urunun karbon sayisi substrattan kucukse) ve
    # substrat oksijen tasiyorsa, substratin KENDI oksijenleri bir parcayla
    # cikmis olabilir; o zaman ana urundeki duz sayim "yeni oksijen yok"
    # demeye yetmez. Antranilat 1,2-dioksijenaz tam bu durumda:
    # antranilat C7H6NO2 (O=2) -> katekol C6H6O2 (O=2) + CO2. Sayim "ana
    # urunde kazanc yok" diyor, oysa katekolun IKI OH'i da O2'den gelir ve
    # antranilatin iki oksijeni CO2 olarak cikar. Dolayisiyla sayim burada
    # EMIN BIR YANLIS uretiyordu.
    #
    # Bu durumda deger hesaplanmaz: "ambiguous" donulur ve kuratorun
    # gerekcesiyle kaydettigi deger ayakta kalir. Atom duzeyinde cozum Rhea'nin
    # atom eslemeli (RXN) verisini ya da 18O deneyini gerektirir.
    if main["C"] < sub["C"] and sub["O"] > 0:
        return {"rhea": rhea_id, "equation": r["equation"],
                "fate": "ambiguous",
                "why": (f"the skeleton fragments ({sub['name']} C{sub['C']} -> "
                        f"{main['name']} C{main['C']}) and the substrate "
                        f"already carries {sub['O']} oxygen, so a flat count "
                        f"in the main product cannot tell a retained oxygen "
                        f"from a new one; atom counting cannot resolve this "
                        f"reaction"),
                "substrate": sub, "main_product": main,
                "other_products": R[1:]}

    if gain_main > 0:
        fate = "product"
        why = (f"the main product {main['name']} ({main['formula']}) carries "
               f"{gain_main} more oxygen than the substrate {sub['name']} "
               f"({sub['formula']})")
    elif total_gain > 0 and others > 0:
        fate = "leaving_group"
        why = (f"the main product {main['name']} ({main['formula']}) gains no "
               f"oxygen, but the product side gains {total_gain} overall, "
               f"carried by the released fragment(s)")
    else:
        fate = "oxidant_only"
        why = (f"no product gains oxygen: {sub['name']} ({sub['formula']}) -> "
               f"{main['name']} ({main['formula']})")
    return {"rhea": rhea_id, "equation": r["equation"], "fate": fate,
            "why": why, "substrate": sub, "main_product": main,
            "other_products": R[1:]}


if __name__ == "__main__":
    import sys
    for rid in (sys.argv[1:] or ["RHEA:17769"]):
        d = derive_fate(rid)
        print(f"{rid}: {d.get('fate')}")
        print("   ", d.get("equation", "")[:110])
        print("   ", d.get("why"))
    save_cache()
