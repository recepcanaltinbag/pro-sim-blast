"""Her referansin reaksiyonunu RESMI kayda karsi dogrula (UniProt/Rhea/EC).

NEDEN BU BETIK VAR. Kolin monooksijenazin oksijen kaderi 2026-10-28'de IKI KEZ
yanlis kaydedildi ve ikisi de ayni hatadan geldi: denge, enzimin BIRAKTIGI
urunden degil, onun dehidratlanmis denge formundan hesaplandi. Bizim
`product_en` alani "betaine aldehyde" diyordu; UniProt/Rhea ise resmi
reaksiyonu "choline + O2 + 2 reduced [2Fe-2S]-[ferredoxin] + 2 H(+) = betaine
aldehyde HYDRATE + H2O + ..." diye veriyor. Urun hidrat olunca bir oksijen
URUNE giriyor ve "oksijen kazanci yok" sonucu curuyor.

Dolayisiyla kural: BIR OKSIJEN KADERI, ONCE URUNUN KIMLIGI RESMI KAYITTAN
DOGRULANMADAN belirlenmez. Bu betik o dogrulamayi 76 tipin tamami icin yapar.
Kendi basina HICBIR SEYI DUZELTMEZ -- uyusmazligi raporlar, karari kurator
verir.

  python3 verify_reactions.py              # rapor
  python3 verify_reactions.py --out X.json # baska dosyaya
"""
import argparse
import glob
import json
import os
import re
import sys
import time
import urllib.request

import yaml

UA = "ROAR-DB/1.0 (mailto:recocanextra@gmail.com)"
OUT = "analysis_out/reaction_verification.json"
UNIPROT_RE = re.compile(r"^(?:[OPQ][0-9][A-Z0-9]{3}[0-9]"
                        r"|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})$")


def get(url, tries=3):
    for i in range(tries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": UA})
            return urllib.request.urlopen(req, timeout=60).read()
        except Exception:
            if i == tries - 1:
                raise
            time.sleep(2 * (i + 1))


def uniprot_reactions(acc):
    """EC numaralari ve KATALITIK AKTIVITE kayitlari (Rhea denklemi dahil)."""
    d = json.loads(get(f"https://rest.uniprot.org/uniprotkb/{acc}.json"))
    pd = d.get("proteinDescription") or {}
    ecs = [e["value"] for e in
           ((pd.get("recommendedName") or {}).get("ecNumbers") or [])]
    for sub in (pd.get("submissionNames") or []):
        ecs += [e["value"] for e in (sub.get("ecNumbers") or [])]
    acts = []
    for c in d.get("comments") or []:
        if c.get("commentType") != "CATALYTIC ACTIVITY":
            continue
        r = c.get("reaction") or {}
        xr = r.get("reactionCrossReferences") or []
        acts.append({
            "equation": r.get("name"),
            "ec": r.get("ecNumber"),
            "rhea": next((x["id"] for x in xr if x["database"] == "Rhea"
                          and x["id"].startswith("RHEA:")), None),
            "chebi": [x["id"] for x in xr if x["database"] == "ChEBI"],
        })
    return {"ec_numbers": sorted(set(ecs)), "activities": acts,
            "reviewed": d.get("entryType")}


def sides(equation):
    """Denklemi iki yana ayir ve katilimcilari listele."""
    if not equation or "=" not in equation:
        return None, None
    left, right = equation.split("=", 1)
    split = lambda s: [p.strip() for p in s.split(" + ") if p.strip()]
    return split(left), split(right)


# Reaksiyonun her iki yaninda da bulunan, URUN SAYILMAYACAK katilimcilar.
COFACTORS = re.compile(
    r"^(?:\d+\s+)?(?:O2|H2O|H\(\+\)|NADH|NAD\(\+\)|NADPH|NADP\(\+\)|"
    r"(?:2 )?(?:reduced|oxidized) \[.*|FAD|FADH2|AH2|A)$", re.I)


def product_candidates(right):
    return [p for p in right if not COFACTORS.match(p)]


def norm(s):
    s = (s or "").lower()
    s = re.sub(r"\(.*?\)", " ", s)
    return {w for w in re.split(r"[^a-z0-9]+", s) if len(w) > 3}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=OUT)
    ap.add_argument("--only", default="")
    a = ap.parse_args()
    want = {c.strip() for c in a.only.split(",") if c.strip()}

    rows = []
    for f in sorted(glob.glob("types/*.yaml")):
        m = yaml.safe_load(open(f))
        cl = m["cluster"]
        if want and cl not in want:
            continue
        acc = str(m.get("accession") or "").strip()
        row = {"cluster": cl, "accession": acc,
               "recorded_substrate": m.get("substrate_en"),
               "recorded_product": m.get("product_en"),
               "recorded_reaction_class": m.get("reaction_class"),
               "recorded_oxygen_fate": m.get("oxygen_fate")}
        if not acc or not UNIPROT_RE.match(acc):
            row["status"] = ("no UniProt accession -- the official reaction "
                             "cannot be looked up this way")
            rows.append(row)
            print(f"[{cl}] atlandi (UniProt accession yok)", file=sys.stderr)
            continue
        try:
            u = uniprot_reactions(acc)
        except Exception as exc:                        # noqa: BLE001
            row["status"] = f"lookup failed: {exc}"
            rows.append(row)
            continue
        row["ec_numbers"] = u["ec_numbers"]
        row["activities"] = u["activities"]
        if not u["activities"]:
            row["status"] = ("the UniProt entry records no catalytic activity, "
                             "so there is no official equation to check "
                             "against")
        else:
            # KARSILASTIRMA IKI YONLU OLMAK ZORUNDA. Ilk yazimda yalnizca
            # ORTAK kelime sayiliyordu ve o kontrol, bu betigi yazdiran
            # hatayi YAKALAMIYORDU: bizim "betaine aldehyde" ile resmi
            # "betaine aldehyde hydrate" iki kelime paylasiyor, yani
            # "uyusuyor" cikiyordu -- oysa eksik olan tek kelime ("hydrate")
            # oksijen kaderini tersine ceviren kelimeydi. Artik RESMIDE OLUP
            # BIZDE OLMAYAN kelimeler de bildiriliyor ve bos olmayan her
            # fark INCELEME gerektirir.
            eqs = []
            ours = norm(m.get("product_en"))
            for act in u["activities"]:
                _left, right = sides(act["equation"])
                prods = product_candidates(right or [])
                scored = []
                for pr in prods:
                    o = norm(pr)
                    scored.append({"official_product": pr,
                                   "shared": sorted(ours & o),
                                   "missing_from_ours": sorted(o - ours)})
                best = max(scored, key=lambda x: (len(x["shared"]),
                                                  -len(x["missing_from_ours"])),
                           default=None)
                eqs.append({"equation": act["equation"], "rhea": act["rhea"],
                            "official_products": prods,
                            "closest": best})
            row["checked"] = eqs
            matched = [e for e in eqs if e["closest"] and e["closest"]["shared"]]
            row["product_agrees"] = bool(matched)
            # Resmi urunde olup bizde olmayan kelimeler -- bos olmayan her
            # kume, oksijen kaderine guvenmeden once okunmak zorunda.
            row["qualifiers_missing_from_our_product"] = sorted(
                {w for e in matched for w in e["closest"]["missing_from_ours"]})
            if not matched:
                row["status"] = ("PRODUCT DOES NOT MATCH any official "
                                 "equation -- review before trusting the "
                                 "oxygen fate")
            elif row["qualifiers_missing_from_our_product"]:
                row["status"] = ("product matches but the official name "
                                 "carries words ours does not: "
                                 + ", ".join(
                                     row["qualifiers_missing_from_our_product"])
                                 + " -- REVIEW")
            else:
                row["status"] = "product name agrees with an official equation"
        rows.append(row)
        print(f"[{cl}] {row['status'][:60]}", file=sys.stderr)
        time.sleep(0.2)

    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    json.dump({"_generated": time.strftime("%Y-%m-%d %H:%M UTC",
                                           time.gmtime()),
               "_about": "Her referansin kayitli urunu, UniProt/Rhea'nin "
                         "resmi reaksiyon denklemine karsi. Duzeltme "
                         "YAPILMAZ.",
               "_why": "Kolin monooksijenazin oksijen kaderi iki kez yanlis "
                       "kaydedildi cunku denge, enzimin biraktigi urun yerine "
                       "dehidratlanmis denge formundan hesaplandi.",
               "rows": rows}, open(a.out, "w"), indent=1, ensure_ascii=False)
    clean = [r["cluster"] for r in rows if r.get("product_agrees")
             and not r.get("qualifiers_missing_from_our_product")]
    review = [(r["cluster"], r["qualifiers_missing_from_our_product"])
              for r in rows if r.get("product_agrees")
              and r.get("qualifiers_missing_from_our_product")]
    bad = [r["cluster"] for r in rows if r.get("product_agrees") is False]
    skip = [r["cluster"] for r in rows if "product_agrees" not in r]
    print(f"\ntemiz: {len(clean)} | INCELE: {len(review)} | UYUSMAYAN: "
          f"{len(bad)} | kontrol edilemeyen: {len(skip)}")
    for cl, w in review:
        print(f"  INCELE {cl}: resmi adda var, bizde yok -> {', '.join(w)}")
    if bad:
        print("UYUSMAYANLAR:", ", ".join(bad))
    print("->", a.out)


if __name__ == "__main__":
    main()
