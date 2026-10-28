"""O2'den gelen oksijenin kaderini tek yerden coz.

NEDEN. "Rieske oksijenaz" bir AILE adidir ve ailenin icinde urune hic oksijen
sokmayan uyeler var: desaturasyon ve C-C bagi kurma, O2'yi yalnizca yukseltgen
olarak kullanir. "Dioksijenaz" diye filtreleyen biri onlari da aliyordu. Bu
modul sinif sozlugunu (reaction_classes.yaml) ve tiplerin kendi beyanlarini
birlestirip her tip icin TEK bir deger verir.

Siniflandirma reaksiyon sinifindan TURETILEMEZ: dioxygenation_with_release
homojen degildir -- nitroarenlerde oksijen ana urune, fenoksiasetatlarda kopan
glioksilata gider. O yuzden sozluk o siniflar icin `per_type` der ve tipin
kendi `oxygen_fate` alanini beyan etmesini ZORUNLU kilar.
"""
import os
import yaml

VOCAB_FILE = "reaction_classes.yaml"
# transient: oksijen gercekten sokulur ama urunde KALMAZ. CMO bunun ornegi:
# kolin hidroksillenip gem-diol olusur, ara urun su kaybeder ve betain
# aldehit cikar -- kolin C5H14NO+ -> betain aldehit C5H12NO+, yani -2H ve
# sifir O. "hydroxylation" sinifinin varsayilani product'tir ve CMO icin
# YANLISTI; deger kullanicinin "CmoS icin Rieske diyorlar ama" sorusu
# uzerine formulle kontrol edilince bulundu.
FATES = ("product", "leaving_group", "transient", "oxidant_only",
         "unknown")

# Filtrede ne gosterilecegi: oksijen AKTARAN enzimler bu ikisidir.
# Oksijeni bir urunde BIRAKAN enzimler. transient buraya girmez: oksijen
# sokulur ama su olarak geri cikar, yani hicbir uronde bulunmaz.
TRANSFERRING = ("product", "leaving_group")


def load_vocab(root="."):
    path = os.path.join(root, VOCAB_FILE)
    data = yaml.safe_load(open(path)) or {}
    return {c["id"]: c for c in (data.get("reaction_classes") or [])}


def resolve(manifest, vocab):
    """(fate, kaynak) dondur. Kaynak: 'type' beyandan, 'class' sozlukten.

    Bilinmeyen sinif ya da eksik beyan SESSIZCE 'product' sayilmaz -- tahmin
    yazmak, bos birakmaktan kotudur. 'unknown' doner ve dogrulama bagirir.
    """
    rc = str(manifest.get("reaction_class") or "").strip()
    declared = str(manifest.get("oxygen_fate") or "").strip()
    cls = vocab.get(rc)
    if cls is None:
        return "unknown", "undefined_class"
    if cls["oxygen_fate"] == "per_type":
        if declared in FATES:
            return declared, "type"
        return "unknown", "missing_declaration"
    if declared:
        # Sinif homojen oldugu halde tip ayri beyan etmis. Bu SESSIZCE
        # gecerli sayilmaz: sozluk o sinif icin `overridable: true` demis VE
        # tip gerekcesini yazmis olmak zorunda. Aksi halde celiski olarak
        # bildirilir ve dogrulama bagirir.
        if declared != cls["oxygen_fate"]:
            ok = bool(cls.get("overridable")) and bool(
                str(manifest.get("oxygen_fate_basis") or "").strip())
            return declared, ("type_override_allowed" if ok
                              else "type_overrides_class")
        return declared, "type_redundant"
    return cls["oxygen_fate"], "class"


def transfers_oxygen(fate):
    return fate in TRANSFERRING


OUT_FILE = "oxygen_fate.csv"
OUT_FIELDS = ["cluster", "reaction_class", "oxygen_fate", "resolved_from",
              "basis"]


def table(root="."):
    """Her tip icin cozulmus satir. chemistry.csv'ye SUTUN EKLENMEDI: o dosyayi
    yirmiden fazla betik okuyor ve sema degisikligi gereksiz risk. Bunun yerine
    turetilmis bir tablo uretilir, bagimlilik grafigine girer ve bayatlarsa
    dogrulama kapisi yakalar."""
    import glob
    vocab = load_vocab(root)
    rows = []
    for f in sorted(glob.glob(os.path.join(root, "types", "*.yaml"))):
        m = yaml.safe_load(open(f))
        fate, src = resolve(m, vocab)
        basis = str(m.get("oxygen_fate_basis") or "").strip()
        if not basis:
            cls = vocab.get(str(m.get("reaction_class") or "").strip())
            basis = (cls or {}).get("evidence", "")
        rows.append({"cluster": m.get("cluster"),
                     "reaction_class": m.get("reaction_class"),
                     "oxygen_fate": fate, "resolved_from": src,
                     "basis": " ".join(str(basis).split())})
    return rows


def main():
    import argparse
    import collections
    import csv
    import glob
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=".")
    ap.add_argument("--write", action="store_true",
                    help=f"{OUT_FILE} dosyasini yaz")
    a = ap.parse_args()
    vocab = load_vocab(a.root)
    counts = collections.Counter()
    srcs = collections.Counter()
    for f in sorted(glob.glob(os.path.join(a.root, "types", "*.yaml"))):
        m = yaml.safe_load(open(f))
        fate, src = resolve(m, vocab)
        counts[fate] += 1
        srcs[src] += 1
        if src in ("missing_declaration", "undefined_class",
                   "type_overrides_class"):
            print(f"  [dikkat] {m.get('cluster')}: {src}")
    print("oksijen kaderi:", dict(counts))
    print("cozum kaynagi :", dict(srcs))
    n = sum(v for k, v in counts.items() if transfers_oxygen(k))
    print(f"oksijen AKTARAN: {n} | aktarmayan: {counts['oxidant_only']} | "
          f"bilinmeyen: {counts['unknown']}")
    if a.write:
        rows = table(a.root)
        path = os.path.join(a.root, OUT_FILE)
        with open(path, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=OUT_FIELDS)
            w.writeheader()
            w.writerows(rows)
        print(f"[yazildi] {path} ({len(rows)} satir)")


if __name__ == "__main__":
    main()
