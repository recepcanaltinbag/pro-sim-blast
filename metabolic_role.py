"""Her tip katabolik mi biyosentetik mi -- tek yerden coz.

NEDEN. Kullanici biyodegradasyon yapmayan uyeleri suzmek istedi. Ilk akla
gelen olcut "oksijen urune giriyor mu" idi ve o YANLIS EKSEN:
  * SxtT, GxtA, PrnD, CAO oksijeni urune sokar ama biyodegradasyon YAPMAZ --
    oksijen kaderine gore suzulse kalirlardi;
  * vanillat, guaiakol, dikamba, kafein, 2,4-D demetilasyonlarinda oksijen
    kopan formaldehite girer ama bunlar klasik biyodegradasyon -- oksijen
    kaderine gore suzulse 13 tip silinirdi.
Iki eksen DIKTIR. Bu modul ikincisini verir.

Rol, enzimin RIESKE OLUP OLMADIGI hakkinda hicbir sey SOYLEMEZ. CmoS
biyosentetiktir ve ayni zamanda tam bir Rieske oksijenazdir (motif kapisini
Rieske 4/4 ile geciyor); etiket yalnizca yolun yonunu kaydeder.
"""
import argparse
import collections
import csv
import os
import yaml

VOCAB_FILE = "metabolic_roles.yaml"
ECOLOGY_FILE = "cluster_ecology.csv"
OUT_FILE = "metabolic_role.csv"
OUT_FIELDS = ["cluster", "substrate_class", "metabolic_role", "resolved_from",
              "basis"]
ROLES = ("catabolic", "biosynthetic", "unknown")


def load_vocab(root="."):
    return yaml.safe_load(open(os.path.join(root, VOCAB_FILE))) or {}


def load_substrate_class(root="."):
    path = os.path.join(root, ECOLOGY_FILE)
    if not os.path.exists(path):
        return {}
    with open(path) as fh:
        return {r["cluster"]: (r.get("substrate_class") or "").strip()
                for r in csv.DictReader(fh)}


def resolve(manifest, vocab, substrate_class):
    """(rol, kaynak) dondur.

    Beyan eksikse SESSIZCE katabolik sayilmaz: 'unknown' doner ve dogrulama
    bagirir. Varsayilan yalnizca rolun TEK ANLAMLI oldugu substrat siniflari
    icin gecerlidir.
    """
    declared = str(manifest.get("metabolic_role") or "").strip()
    defaults = vocab.get("defaults_by_substrate_class") or {}
    rule = defaults.get(substrate_class)
    if rule == "per_type":
        return (declared, "type") if declared in ROLES else (
            "unknown", "missing_declaration")
    if declared:
        if rule and declared != rule:
            ok = bool(str(manifest.get("metabolic_role_basis") or "").strip())
            return declared, ("type_override_allowed" if ok
                              else "type_overrides_default")
        return declared, "type"
    if rule in ROLES:
        return rule, "substrate_class"
    return "unknown", "no_rule"


def table(root="."):
    import glob
    vocab = load_vocab(root)
    sclass = load_substrate_class(root)
    rows = []
    for f in sorted(glob.glob(os.path.join(root, "types", "*.yaml"))):
        m = yaml.safe_load(open(f))
        cl = m.get("cluster")
        sc = sclass.get(cl, "")
        role, src = resolve(m, vocab, sc)
        basis = str(m.get("metabolic_role_basis") or "").strip()
        if not basis and src == "substrate_class":
            basis = (f"assigned from the substrate class '{sc}', for which "
                     "the direction is unambiguous: the pathway degrades the "
                     "compound.")
        rows.append({"cluster": cl, "substrate_class": sc,
                     "metabolic_role": role, "resolved_from": src,
                     "basis": " ".join(basis.split())})
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=".")
    ap.add_argument("--write", action="store_true")
    a = ap.parse_args()
    rows = table(a.root)
    print("rol:", dict(collections.Counter(r["metabolic_role"] for r in rows)))
    print("kaynak:", dict(collections.Counter(r["resolved_from"]
                                              for r in rows)))
    for r in rows:
        if r["resolved_from"] in ("missing_declaration", "no_rule",
                                  "type_overrides_default"):
            print(f"  [dikkat] {r['cluster']}: {r['resolved_from']}")
    if a.write:
        path = os.path.join(a.root, OUT_FILE)
        with open(path, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=OUT_FIELDS)
            w.writeheader()
            w.writerows(rows)
        print(f"[yazildi] {path} ({len(rows)} satir)")


if __name__ == "__main__":
    main()
