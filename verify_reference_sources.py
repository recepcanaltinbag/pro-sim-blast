"""Kaynak alanindaki "verified byte-identical" IDDIASINI sina.

NEDEN BU BETIK VAR. 76 kuratorlu referansin 35'inin `source` alaninda
"stored sequence verified byte-identical to ..." yaziyor. Ama bu bir CUMLE --
bir insanin yazdigi iddia. Dizi sonradan degisirse ya da yanlis accession
yazildiysa cumle oldugu yerde durur ve kimse fark etmez. Bu depoda tam bu
sinif hata birkac kez olculdu: sabit referans sayilari veri degisince
yanlisa dustu, motif kolonlari iki yerde tanimliydi.

Burada iddia kaynagina karsi YENIDEN cozulur: accession cekilir, dizi
karsilastirilir. Iddia yanlissa kontrol kirilir.

Kullanim:
    python3 verify_reference_sources.py            # orneklem
    python3 verify_reference_sources.py --all      # hepsi (yavas, ag erisimi)
"""
import argparse, glob, json, os, random, re, sys, time
import urllib.parse, urllib.request
import yaml

UA = "ROAR-DB-provenance (mailto:recocanextra@gmail.com)"
CLAIM = "verified byte-identical"


def uniprot(acc):
    u = f"https://rest.uniprot.org/uniprotkb/{acc}.fasta"
    t = urllib.request.urlopen(urllib.request.Request(u, headers={"User-Agent": UA}),
                               timeout=90).read().decode()
    return "".join(t.split("\n")[1:]).strip()


def genbank(acc):
    q = {"db": "protein", "id": acc, "rettype": "fasta", "retmode": "text",
         "tool": "ROAR-DB-provenance", "email": "recocanextra@gmail.com"}
    t = urllib.request.urlopen(
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?"
        + urllib.parse.urlencode(q), timeout=120).read().decode()
    return "".join(l.strip() for l in t.split("\n")[1:] if l.strip())


def pdb(code):
    """Yatirilan polimer varliklarinin dizileri."""
    base = "https://data.rcsb.org/rest/v1/core"
    e = json.load(urllib.request.urlopen(f"{base}/entry/{code}", timeout=90))
    out = []
    for pid in e.get("rcsb_entry_container_identifiers", {}).get("polymer_entity_ids", []):
        p = json.load(urllib.request.urlopen(f"{base}/polymer_entity/{code}/{pid}", timeout=90))
        s = p.get("entity_poly", {}).get("pdbx_seq_one_letter_code_can", "").replace("\n", "")
        if s:
            out.append(s)
    return out


def resolve(src):
    """Kaynak metninden (veritabani, accession) cikar."""
    m = re.search(r"\bPDB\s+([0-9A-Za-z]{4})\b", src)
    if m:
        return "PDB", m.group(1).upper()
    m = re.search(r"\bGenBank\s+([A-Z]{3}[0-9]{5,7}\.[0-9])\b", src)
    if m:
        return "GenBank", m.group(1)
    m = re.search(r"\bUniProt\s+([A-Z0-9]{6,10})\b", src)
    if m:
        return "UniProt", m.group(1)
    return None, None


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--types-dir", default="types")
    ap.add_argument("--sample", type=int, default=8)
    ap.add_argument("--all", action="store_true")
    a = ap.parse_args()

    claims = []
    for f in sorted(glob.glob(os.path.join(a.types_dir, "*.yaml"))):
        d = yaml.safe_load(open(f)) or {}
        src = str(d.get("source") or "")
        if CLAIM in src:
            db, acc = resolve(src)
            claims.append((d["cluster"], db, acc, d["sequence"]))
    print(f"[iddia] {len(claims)} referansin kaynaginda '{CLAIM}' yaziyor")

    unresolvable = [(c, s) for c, db, acc, s in claims if not acc]
    if unresolvable:
        print(f"[HATA] {len(unresolvable)} iddiada accession AYRISTIRILAMADI: "
              + ", ".join(c for c, _ in unresolvable))

    todo = claims if a.all else random.Random(0).sample(claims, min(a.sample, len(claims)))
    todo = [t for t in todo if t[2]]
    print(f"[sinanan] {len(todo)}")
    print()
    bad = 1 if unresolvable else 0
    for cluster, db, acc, ours in todo:
        try:
            if db == "UniProt":
                got = [uniprot(acc)]
            elif db == "GenBank":
                got = [genbank(acc)]
            else:
                got = pdb(acc)
        except Exception as exc:                      # noqa: BLE001
            print(f"[COZULEMEDI] {cluster:16s} {db} {acc}  ({type(exc).__name__})")
            bad = 1
            time.sleep(0.4)
            continue
        if any(s == ours for s in got):
            print(f"[OK   ] {cluster:16s} {db} {acc}  {len(ours)} aa")
        else:
            lens = ", ".join(str(len(s)) for s in got) or "-"
            print(f"[YANLIS] {cluster:16s} {db} {acc}  iddia birebir, ama DEGIL "
                  f"(bizde {len(ours)} aa, kaynakta {lens})")
            bad = 1
        time.sleep(0.4)

    print()
    if bad:
        print("IDDIA DOGRULANAMADI -- source alanindaki cumle kaynakla ortusmuyor.")
        return 1
    print("Butun sinanan iddialar dogru: kaynaktaki dizi bizdekiyle birebir ayni.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
