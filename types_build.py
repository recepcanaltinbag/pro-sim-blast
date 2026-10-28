"""`types/` dizinini TEK KAYNAK yapar: chemistry.csv ve referans fasta ondan URETILIR.

NEDEN. Bir enzim tipi su anda iki ayri yerde tanimli -- kimyasi
`chemistry.csv`'de, dizisi `ROs_71_Clean/refs71.fasta`'da -- ve ikisini
birbirine baglayan tek sey bir adlandirma gelenegi. Hicbir sey ikisinin
tutarli oldugunu dogrulamiyor, ve kimin ne zaman hangi makaleden ekledigi
hicbir yerde yazmiyor. Bu, tek basina calisan biri icin katlanilir; baskalari
yeni Rieske'ler ekleyecekse degil.

Bundan sonra bir tip TEK dosyadir: `types/<kume>.yaml`. Dizi, kimya, yapi,
kaynak, ekleyen ve tarih orada. `chemistry.csv` ve referans fasta bu dizinden
URETILIR ve elle duzenlenmez.

GECIS GUVENLIGI. Kaynak degistirmek sessizce veri kaybettirebilir, bu yuzden
arac once GIDIS-DONUS testi yapar: manifestolardan uretilen chemistry.csv
mevcut dosyayla hucre hucre karsilastirilir ve TEK bir fark bile varsa yazma
reddedilir. Yirmi betik bu dosyayi okuyor; uretilen surum birebir ayni
degilse gecis yapilmaz.

BUTUN SETI DOGRULAR, tek tipi degil. Tek tip dogrulamasi `add_type.py`'de;
burada yalnizca kumenin TAMAMINA bakilarak gorulebilecek hatalar aranir:

    ayni kume kimligi iki dosyada
    ayni DIZI iki tipte           <- bu veritabanindaki 4 bos tipin sebebi
    kume adindaki grup ile group alani uyusmuyor
    zorunlu alan bos

Kullanim:
    python3 types_build.py --check          # dogrula + gidis-donus testi
    python3 types_build.py --write          # chemistry.csv ve fasta uret
"""

import argparse
import csv
import difflib
import hashlib
import os
import re
import shutil
import sys
from collections import Counter, defaultdict

TYPES = "types"
CHEMISTRY = "chemistry.csv"
REFS = "ROs_71_Clean/refs71.fasta"

CLUSTER_RE = re.compile(r"^(\d+)_(\d+)_([A-Za-z0-9_]+)$")
GROUPS_FILE = "groups.yaml"


CITATIONS_FILE = "citations.yaml"


def known_citations(path=CITATIONS_FILE):
    """citations.yaml icindeki atif anahtarlari. Dosya yoksa bos kume ve
    atif kontrolu yapilmaz -- davranis eskisi gibi kalir."""
    if not os.path.exists(path):
        return set()
    try:
        import yaml
        data = yaml.safe_load(open(path)) or {}
        return {str(c.get("id")) for c in (data.get("citations") or [])
                if c.get("id")}
    except Exception:
        return set()


def known_groups(path=GROUPS_FILE):
    """groups.yaml'daki grup kimlikleri. Dosya yoksa bos kume -- dogrulama
    o zaman grup kontrolu yapmaz ve davranis eskisi gibi kalir."""
    if not os.path.exists(path):
        return set()
    ids = set()
    try:
        import yaml
        data = yaml.safe_load(open(path)) or {}
        for g in data.get("groups") or []:
            if g.get("id") is not None:
                ids.add(str(g["id"]))
        return ids
    except ImportError:
        pass
    for line in open(path):
        m = re.match(r"\s*-\s*id:\s*(\d+)", line)
        if m:
            ids.add(m.group(1))
    return ids
CHEMISTRY_FIELDS = ["cluster", "substrate_en", "substrate_smiles", "product_en",
                    "reaction", "reaction_class", "family", "pdb", "source",
                    "source_kind", "curation_confidence", "notes"]
REQUIRED = ["cluster", "group", "substrate_en", "product_en", "reaction",
            "reaction_class", "source", "source_kind"]
DUPLICATE_IDENTITY = 99.0


def load_manifest(path):
    """Duz `anahtar: deger`. PyYAML varsa o, yoksa yerlesik ayristirici.

    Bagimlilik SART KOSULMAZ: bir katkici icin kurulum engeli olmamali.
    """
    text = open(path).read()
    try:
        import yaml
        data = yaml.safe_load(text)
        if isinstance(data, dict):
            # LISTELER KORUNUR. Eskiden her deger str() ile duzlestiriliyordu
            # ve `citations:` listesi "['a', 'b']" metnine donuyordu; o metin
            # harf harf gezilince her karakter bir atif anahtari sayiliyordu.
            # Liste disindaki her sey metin kalir -- chemistry.csv gidis-donusu
            # buna dayaniyor.
            return {k: ("" if v is None else v if isinstance(v, list)
                        else str(v)) for k, v in data.items()}
    except ImportError:
        pass
    data, key = {}, None
    for line in text.splitlines():
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        if line[0] not in " \t" and ":" in line:
            key, _, value = line.partition(":")
            key = key.strip()
            data[key] = value.strip().strip('"').strip("'")
        elif key and line.lstrip().startswith("- "):
            # PyYAML yoksa da liste alanlari liste olarak okunur.
            item = line.lstrip()[2:].strip().strip('"').strip("'")
            if not isinstance(data.get(key), list):
                data[key] = [] if not str(data.get(key) or "").strip() else [data[key]]
            data[key].append(item)
        elif key and isinstance(data.get(key), str):
            data[key] = (data[key] + line.strip()).strip()
    return data


def read_all(types_dir):
    out = {}
    for path in sorted(os.listdir(types_dir)):
        if not path.endswith((".yaml", ".yml", ".json")):
            continue
        m = load_manifest(os.path.join(types_dir, path))
        m["_file"] = path
        out[path] = m
    return out


def validate_set(manifests, groups=None):
    """Yalnizca BUTUN sete bakilarak gorulebilecek hatalar."""
    errors, warnings = [], []
    groups = known_groups() if groups is None else groups
    citations = known_citations()

    by_cluster = defaultdict(list)
    for fname, m in manifests.items():
        by_cluster[(m.get("cluster") or "").strip()].append(fname)
    for cluster, files in sorted(by_cluster.items()):
        if len(files) > 1:
            errors.append(f"cluster '{cluster}' {len(files)} dosyada: "
                          f"{', '.join(files)}")
        if not cluster:
            errors.append(f"cluster adi bos: {', '.join(files)}")
            continue
        mm = CLUSTER_RE.match(cluster)
        if not mm:
            errors.append(f"cluster adi kalibi tutmuyor: {cluster}")
            continue
        m = manifests[files[0]]
        if m.get("group") and str(m["group"]).strip() != mm.group(1):
            errors.append(f"{cluster}: group alani {m['group']} ama kume adi "
                          f"{mm.group(1)} diyor")
        # Grup groups.yaml'da TANIMLI olmali. Yoksa kume kimliginin ilk hanesi
        # sessizce yeni bir grup icat eder ve hicbir yerde ne oldugu yazmaz --
        # gruplarin bugune kadarki hali tam olarak buydu.
        if groups and mm.group(1) not in groups:
            errors.append(f"{cluster}: grup {mm.group(1)} {GROUPS_FILE} "
                          f"icinde tanimli degil. Once orada tanimla "
                          f"(etiket ve gerekcesiyle), sonra tip ekle")
        for field in REQUIRED:
            if not str(m.get(field) or "").strip():
                errors.append(f"{cluster}: zorunlu alan bos -- {field}")
        # COZULEBILIR TANIMLAYICI. Desen once yalnizca DOI/PMID/PMC/PDB
        # ariyordu ve UniProt ile GenBank accession'larini TANIMIYORDU; oysa
        # ikisi de cozulebilir. Sonucu: 23 tipe UniProt accession'i eklendikten
        # sonra bile uyari devam etti (2026-10-28). Artik ikisi de sayiliyor.
        #   UniProt : P12345 / Q6UP13 / A0A059U1S8
        #   GenBank : AAQ91914.2 / AEC21006.1 / WVH02003.1
        # 2026-10-28: kontrol artik SADECE `source` metnine bakmiyor.
        # Atiflar citations.yaml'a tasindi ve tipler oraya `citations:`
        # anahtar listesiyle baglanir; dizinin geldigi kayit ise `accession:`
        # alaninda durur. Yalnizca prozayi tarayan eski kontrol, 67 tipe
        # yapisal atif eklendikten SONRA bile uyari vermeye devam ediyordu --
        # yani tanimlayici vardi, kontrol onu okumuyordu.
        cits = [str(c) for c in (m.get("citations") or []) if str(c).strip()]
        for c in cits:
            if citations and c not in citations:
                errors.append(f"{cluster}: citations icinde "
                              f"{CITATIONS_FILE} tarafindan taninmayan "
                              f"anahtar -- {c}")
        has_structured = bool(str(m.get("accession") or "").strip()) or bool(
            [c for c in cits if not citations or c in citations])
        has_prose = bool(re.search(
            r"10\.\d{4,}/|PMID|PMC|PDB|\b\d{4};"
            r"|\b[OPQ][0-9][A-Z0-9]{3}[0-9]\b"
            r"|\b[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2}\b"
            r"|\b[A-Z]{3}[0-9]{5,7}\.[0-9]\b",
            str(m.get("source") or "")))
        if not has_structured and not has_prose:
            warnings.append(f"{cluster}: cozulebilir bir tanimlayici yok -- "
                            "ne citations/accession alani, ne source "
                            "icinde DOI/PMID/PDB/UniProt/GenBank")

    # AYNI DIZI: bu veritabanindaki dort bos tipin dogrudan sebebi. Ayni
    # diziye sahip iki referans yarismayi beraberlikle bitirir ve biri keyfi
    # olarak butun uyeleri alir.
    by_seq = defaultdict(list)
    for fname, m in manifests.items():
        seq = (m.get("sequence") or "").strip().upper()
        if seq:
            by_seq[seq].append((m.get("cluster"), fname))
    for seq, holders in by_seq.items():
        if len(holders) > 1:
            names = ", ".join(c for c, _ in holders)
            errors.append(f"AYNI DIZI {len(holders)} tipte: {names}. Iki "
                          "referans ayni diziye sahipse yarisma beraberlikle "
                          "biter ve biri keyfi olarak butun uyeleri alir -- bu "
                          "veritabanindaki bos tiplerin sebebi budur")
    # alt dizi olma durumu (NdmC/NdmB vakasi: biri otekinin etiketi kirpilmis hali)
    seqs = [(m.get("cluster"), (m.get("sequence") or "").strip().upper())
            for m in manifests.values()]
    seqs = [(c, s) for c, s in seqs if s]
    for i, (ca, sa) in enumerate(seqs):
        for cb, sb in seqs[i + 1:]:
            if sa == sb:
                continue
            if sa in sb or sb in sa:
                short, long_ = (ca, cb) if len(sa) < len(sb) else (cb, ca)
                errors.append(f"{short} dizisi {long_} dizisinin TAM ALT "
                              f"DIZISI -- biri otekinin kirpilmis hali olabilir "
                              f"(bu veritabaninda NdmC/NdmB tam olarak boyleydi)")
    return errors, warnings


def to_chemistry_rows(manifests):
    rows = []
    for m in manifests.values():
        rows.append({f: str(m.get(f, "") or "") for f in CHEMISTRY_FIELDS})
    rows.sort(key=lambda r: r["cluster"])
    return rows


def write_csv(rows, path, fields):
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)


def roundtrip(manifests, chemistry_path):
    """Uretilen chemistry.csv mevcut dosyayla birebir ayni mi?"""
    if not os.path.exists(chemistry_path):
        return [], []
    with open(chemistry_path) as fh:
        reader = csv.DictReader(fh)
        fields = list(reader.fieldnames)
        current = {r["cluster"]: r for r in reader}
    generated = {r["cluster"]: r for r in to_chemistry_rows(manifests)}
    diffs = []
    for cluster in sorted(set(current) | set(generated)):
        a, b = current.get(cluster), generated.get(cluster)
        if a is None:
            diffs.append((cluster, "(whole row)", "missing from chemistry.csv", "would be added"))
            continue
        if b is None:
            diffs.append((cluster, "(whole row)", "present in chemistry.csv", "missing from types/"))
            continue
        for f in fields:
            x, y = (a.get(f) or "").strip(), (b.get(f) or "").strip()
            if x != y:
                diffs.append((cluster, f, x, y))
    return diffs, fields


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--types", default=TYPES)
    ap.add_argument("--chemistry", default=CHEMISTRY)
    ap.add_argument("--refs", default=REFS)
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--check", action="store_true",
                   help="dogrula ve gidis-donus testi yap; yazma")
    g.add_argument("--write", action="store_true",
                   help="chemistry.csv ve referans fastayi uret")
    ap.add_argument("--allow-drift", action="store_true",
                    help="gidis-donus farki olsa bile yaz. Yalnizca farkin "
                         "KASITLI oldugunu dogruladiysan")
    args = ap.parse_args()

    manifests = read_all(args.types)
    print(f"{len(manifests)} manifesto okundu: {args.types}/")

    errors, warnings = validate_set(manifests)
    for w in warnings:
        print(f"  [uyari] {w}")
    if errors:
        print()
        print(f"{len(errors)} HATA -- hicbir sey yazilmadi:")
        for e in errors:
            print(f"  - {e}")
        return 1
    print("  set dogrulamasi gecti: kimlik cakismasi yok, ayni dizi yok, "
          "alt dizi yok, zorunlu alanlar dolu")

    diffs, fields = roundtrip(manifests, args.chemistry)
    print()
    if diffs:
        print(f"GIDIS-DONUS FARKI: {len(diffs)} hucre")
        for cluster, field, cur, gen in diffs[:25]:
            print(f"  {cluster}  {field}")
            print(f"      chemistry.csv: {cur[:80]!r}")
            print(f"      types/:        {gen[:80]!r}")
        if len(diffs) > 25:
            print(f"  ... ve {len(diffs) - 25} fark daha")
    else:
        print("Gidis-donus TEMIZ: manifestolardan uretilen chemistry.csv "
              "mevcut dosyayla birebir ayni.")

    if args.check:
        print()
        if diffs:
            print("Fark varken kaynak degistirilemez. Once manifestolari "
                  "duzelt, ya da farkin kasitli oldugunu dogrulayip "
                  "--write --allow-drift kullan.")
            return 1
        print("types/ dizini kaynak olarak kullanilabilir.")
        return 0

    if diffs and not args.allow_drift:
        print()
        print("YAZILMADI: gidis-donus farki var. Yirmi betik chemistry.csv'yi "
              "okuyor; uretilen surum birebir ayni olmadan gecis yapilmaz.")
        return 1

    for path in (args.chemistry, args.refs):
        if os.path.exists(path):
            shutil.copy2(path, path + ".bak")
    write_csv(to_chemistry_rows(manifests), args.chemistry,
              fields or CHEMISTRY_FIELDS)
    with open(args.refs, "w") as fh:
        for m in sorted(manifests.values(), key=lambda x: x.get("cluster") or ""):
            seq = (m.get("sequence") or "").strip().upper()
            if not seq:
                continue
            gene = m.get("gene") or (m.get("cluster") or "").split("_")[-1]
            fh.write(f">{m['cluster']}_{gene}_pro\n")
            for i in range(0, len(seq), 60):
                fh.write(seq[i:i + 60] + "\n")
    n_seq = sum(1 for m in manifests.values() if (m.get("sequence") or "").strip())
    print()
    print(f"[yazildi] {args.chemistry} ({len(manifests)} satir)")
    print(f"[yazildi] {args.refs} ({n_seq} dizi)")
    print(f"[yedek]   her ikisinin .bak dosyasi")
    print()
    print("Referans seti degistiyse sirayla: python3 build_reference_models.py, "
          "sonra python3 rebuild.py --run --include-heavy")
    return 0


if __name__ == "__main__":
    sys.exit(main())
