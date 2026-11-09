"""
Pfam modellerini AD ile degil ERISIM NUMARASI ile indirip kutuphane kurar.

NEDEN AYRI KUTUPHANE: `component_hmms/components.hmm` operon bilesenlerini
(beta, ferredoksin, reduktaz) tanimlar ve `assign_component` o kutuphanedeki
HER modeli bir bilesen etiketine baglar. Hareketlilik ya da duzenleyici
modelleri oraya karistirilirsa bilesen taramasinin anlami degisir: ayni
domtbl iki farkli soruya cevap vermis olur. Bu yuzden her soru kendi
kutuphanesinde durur ve kendi tablosunu yazar.

AD DEGIL ERISIM: model adlari Pfam surumleri arasinda degisebilir
(ornegin HTH_1 <-> LysR). Indirme PFxxxxx ile yapilir, inen dosyanin
NAME/ACC satirlari okunup manifestoya yazilir; beklenen ad tutmazsa
hata verilir, sessizce yanlis model kullanilmaz.

Cikti: family_hmms/<kutuphane>.hmm + family_hmms/<kutuphane>.manifest.json
"""

import argparse
import gzip
import hashlib
import json
import os
import sys
import urllib.request

PFAM_URL = "https://www.ebi.ac.uk/interpro/wwwapi/entry/pfam/{acc}/?annotation=hmm"
OUT_DIR = "family_hmms"

# (erisim, beklenen ad, ne ise yaradigi)
LIBRARIES = {
    # HAREKETLI ELEMAN: kendini kesip tasiyan transpozaz cekirdekleri.
    "mobility": [
        ("PF01609", "DDE_Tnp_1", "IS4/IS5 benzeri transpozaz DDE cekirdegi"),
        ("PF13843", "DDE_Tnp_1_7", "DDE_Tnp_1 ile birlikte gorulen transpozaz alani"),
        ("PF00665", "rve", "integraz cekirdegi (IS3, IS481, retroviral)"),
        ("PF01527", "HTH_Tnp_1", "IS3/IS150 transpozazin HTH ucu"),
        ("PF00872", "Transposase_mut", "MULE/Mutator tipi transpozaz"),
        ("PF03050", "DDE_Tnp_IS66", "IS66 ailesi transpozaz"),
        ("PF01526", "DDE_Tnp_Tn3", "Tn3 ailesi transpozaz"),
        ("PF02371", "Transposase_20", "IS110/IS492 ailesi transpozaz"),
        ("PF01797", "Y1_Tnp", "IS200/IS605 ailesi transpozaz"),
        ("PF01359", "Transposase_1", "Tc1/mariner (IS630) transpozaz"),
        ("PF03811", "Zn_ribbon_InsA", "IS1 transpozaz InsA"),
        ("PF05717", "TnpB_IS66", "IS66 yardimci proteini TnpB"),
    ],
    # REKOMBINAZ: DNA'yi kesip baglar ama kendisi hareketli eleman DEGILDIR.
    # Ayri kutuphanede durmasinin tek sebebi budur: urun adi regex'i ikisini
    # ayni kovaya koyuyor (XerC, RuvX, "recombinase-like HTH"), ayirmak icin
    # once ayri olculmeleri gerekir.
    "recombinase": [
        ("PF00589", "Phage_integrase", "tirozin rekombinaz / bolgeye ozgu integraz"),
        ("PF00239", "Resolvase", "serin rekombinaz cekirdegi"),
        ("PF02796", "HTH_7", "rekombinaz benzeri HTH ucu"),
        ("PF07508", "Recombinase", "rekombinaz ailesi orta alani"),
        ("PF13408", "Zn_ribbon_recom", "rekombinaz Zn-serit alani"),
        ("PF02075", "RuvC", "Holliday kavsagi cozucu RuvC/RuvX"),
        ("PF03652", "RuvX", "Holliday kavsagi cozucu RuvX/YqgF"),
    ],
    "regulator": [
        ("PF00440", "TetR_N", "TetR ailesi DNA baglanma alani"),
        ("PF00392", "GntR", "GntR ailesi HTH"),
        ("PF00126", "HTH_1", "LysR ailesi HTH ucu"),
        ("PF03466", "LysR_substrate", "LysR ailesi efektor baglanma alani"),
        ("PF01047", "MarR", "MarR ailesi HTH"),
        ("PF02311", "AraC_binding", "AraC ailesi efektor baglanma alani"),
        ("PF12833", "HTH_18", "AraC ailesi HTH ucu"),
        ("PF09339", "HTH_IclR", "IclR ailesi HTH"),
        ("PF01614", "IclR_C", "IclR ailesi efektor baglanma alani"),
    ],
}


def fetch(acc):
    url = PFAM_URL.format(acc=acc)
    with urllib.request.urlopen(url, timeout=120) as response:
        raw = response.read()
    if raw[:2] == b"\x1f\x8b":
        raw = gzip.decompress(raw)
    return raw.decode()


def header(text, key):
    for line in text.splitlines():
        if line.startswith(key):
            return line.split(None, 1)[1].strip()
    return None


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--library", choices=sorted(LIBRARIES) + ["all"],
                    required=True)
    ap.add_argument("--out-dir", default=OUT_DIR)
    ap.add_argument("--force", action="store_true",
                    help="tek tek inmis modeller varsa bile yeniden indir")
    args = ap.parse_args()

    if args.library == "all":
        for lib in sorted(LIBRARIES):
            build(lib, args.out_dir, args.force)
        return
    build(args.library, args.out_dir, args.force)


def build(library, out_dir, force):
    models = LIBRARIES[library]
    os.makedirs(out_dir, exist_ok=True)
    cache = os.path.join(out_dir, "pfam_cache")
    os.makedirs(cache, exist_ok=True)

    manifest, chunks, problems = [], [], []
    for acc, expected, purpose in models:
        path = os.path.join(cache, f"{acc}.hmm")
        if force or not os.path.exists(path) or os.path.getsize(path) == 0:
            print(f"[iniyor] {acc} ({expected})")
            try:
                text = fetch(acc)
            except Exception as exc:                      # noqa: BLE001
                problems.append(f"{acc}: inmedi ({exc})")
                continue
            with open(path, "w") as fh:
                fh.write(text if text.endswith("\n") else text + "\n")
        text = open(path).read()
        name = header(text, "NAME")
        full_acc = header(text, "ACC") or ""
        if name != expected:
            problems.append(f"{acc}: beklenen ad {expected}, inen ad {name}")
            continue
        if not full_acc.startswith(acc):
            problems.append(f"{acc}: ACC satiri {full_acc}")
            continue
        manifest.append({"accession": full_acc, "name": name, "purpose": purpose,
                         "length": int(header(text, "LENG") or 0),
                         "description": header(text, "DESC")})
        chunks.append(text if text.endswith("\n") else text + "\n")

    if problems:
        for p in problems:
            sys.stderr.write(f"[HATA] {p}\n")
        sys.exit(1)

    hmm_path = os.path.join(out_dir, f"{library}.hmm")
    body = "".join(chunks)
    digest = hashlib.sha256(body.encode()).hexdigest()[:16]
    # Icerik aynıysa dosyaya dokunma: mtime tetikler, ozet karar verir.
    if os.path.exists(hmm_path) and \
            hashlib.sha256(open(hmm_path).read().encode()).hexdigest()[:16] == digest:
        print(f"[korundu] {hmm_path} degismedi ({len(manifest)} model)")
    else:
        with open(hmm_path, "w") as fh:
            fh.write(body)
        print(f"[yazildi] {hmm_path} ({len(manifest)} model)")

    man_path = os.path.join(out_dir, f"{library}.manifest.json")
    text = json.dumps({"library": library, "digest": digest,
                       "models": manifest}, indent=1)
    # Manifesto da ayni kurala tabi: icerik aynıysa DOKUNMA. Her kosuda
    # yeniden yazilirsa mtime'i oynar ve onu girdi sayan her adim sonsuza
    # kadar bayat gorunur -- "mtime tetikler, ozet karar verir"in tersi.
    if os.path.exists(man_path) and open(man_path).read() == text:
        print(f"[korundu] {man_path} degismedi")
    else:
        with open(man_path, "w") as fh:
            fh.write(text)
        print(f"[yazildi] {man_path}")
    for m in manifest:
        print(f"   {m['accession']:<12} {m['name']:<16} {m['length']:>4} aa  {m['purpose']}")


if __name__ == "__main__":
    main()
