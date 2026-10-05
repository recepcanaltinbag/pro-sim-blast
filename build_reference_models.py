"""
Referans modellerini yeniden kur -- yeni bir literatur enzimi eklendiginde.

NEDEN BU SCRIPT VAR
  Referans seti sabit degil: literaturde karakterize edilmis ve bu sette olmayan
  RO'lar var, kullanici bunlari ekleyebilmeli. Ama eklemek tek dosya degistirmek
  degildi; su adimlarin hepsi elle yapiliyordu:
     1. refs.fasta'ya diziyi ekle
     2. hizalamayi yeniden kur (clustalo)
     3. ROmotif.hmm'i yeniden kur (hmmbuild)
     4. VE EN KRITIK OLANI: ro_motif.py icindeki 8 motif KOLON NUMARASINI
        elle yeniden cikar, cunku hizalama degisince kolonlar kayar
  Dorduncu adim hem zahmetli hem hata kaynagiydi: kolonlar kodda SABIT yaziliydi
  ve referans seti degisirse sessizce yanlis pozisyonlara bakilirdi.

BU SCRIPT dort adimi birlestirir ve kolonlari VERIDEN cikarir:
  - refs.fasta -> clustalo -> hmmbuild -> ROmotif<N>.hmm
  - referanslari kendi modeline hmmalign'layip her motif icin EN KORUNMUS
    kolonu arar (Rieske C-x-H...C-x-x-H ve katalitik 2-His-1-karboksilat)
  - sonucu ROs_<N>_Clean/motif_columns.json dosyasina yazar
  - ro_motif.py bu dosya varsa onu okur, yoksa gomulu varsayilanlari kullanir
    (geriye donuk uyumluluk: dosya yoksa davranis birebir ayni kalir)

KOLON ARAMA MANTIGI
  Motifler bagil yapiyla tanimlidir, mutlak pozisyonla degil:
     Rieske:    C-x-H ... (14-18 kalinti) ... C-x-x-H
     katalitik: H ... (5 kalinti) ... H ... (uzak) ... D/E
  Script once her kolonda beklenen kalintinin korunma oranini hesaplar, sonra
  bu bagil kisitlari saglayan kolon dortlusu/uclusunu arar ve toplam korunmayi
  en yuksek yapan kombinasyonu secer. Boylece "en korunmus C" gibi tek basina
  yaniltici bir olcute dayanilmaz.

Kullanim:
    python3 build_reference_models.py --refs ROs_71_Clean/refs71.fasta \\
        --out-dir ROs_71_Clean --name ROmotif71 [--skip-align]
    python3 build_reference_models.py --derive-only --out-dir ROs_71_Clean \\
        --name ROmotif71     # modeli yeniden kurmadan sadece kolonlari cikar
"""

import argparse
import json
import os
import subprocess
import sys
from collections import Counter

from ro_motif import read_stockholm_matchcols

# Bagil kisitlar: referans hizalamasinda olculen araliklar, genis tutuldu.
RIESKE_CXH_GAP = (1, 1)        # C ile H arasi kolon sayisi (C-x-H)
RIESKE_CXXH_GAP = (2, 2)       # ikinci motif C-x-x-H
RIESKE_BETWEEN = (10, 24)      # birinci H ile ikinci C arasi
CATALYTIC_HH = (3, 8)          # iki katalitik His arasi
CATALYTIC_TAIL = (100, 220)    # ikinci His ile karboksilat arasi
MIN_CONSERVATION = 0.70        # bir kolonun aday olabilmesi icin alt sinir


def run(cmd, **kw):
    print("[calisiyor]", " ".join(cmd))
    subprocess.run(cmd, check=True, **kw)


def column_conservation(aligned, residues):
    """Her kolon icin verilen kalinti kumesinin orani."""
    names = list(aligned)
    length = len(aligned[names[0]])
    out = []
    for i in range(length):
        hits = sum(1 for n in names if aligned[n][i].upper() in residues)
        out.append(hits / len(names))
    return out


def pick_rieske(cys, his):
    """C-x-H ... C-x-x-H kisitlarini saglayan en korunmus dortluyu bul."""
    best, best_score = None, -1.0
    cys_candidates = [i for i, v in enumerate(cys) if v >= MIN_CONSERVATION]
    his_candidates = [i for i, v in enumerate(his) if v >= MIN_CONSERVATION]
    his_set = set(his_candidates)
    for c1 in cys_candidates:
        for gap1 in range(RIESKE_CXH_GAP[0], RIESKE_CXH_GAP[1] + 1):
            h1 = c1 + gap1 + 1
            if h1 not in his_set:
                continue
            for c2 in cys_candidates:
                sep = c2 - h1
                if not (RIESKE_BETWEEN[0] <= sep <= RIESKE_BETWEEN[1]):
                    continue
                for gap2 in range(RIESKE_CXXH_GAP[0], RIESKE_CXXH_GAP[1] + 1):
                    h2 = c2 + gap2 + 1
                    if h2 not in his_set:
                        continue
                    score = cys[c1] + his[h1] + cys[c2] + his[h2]
                    if score > best_score:
                        best_score, best = score, (c1, h1, c2, h2)
    return best, best_score


def pick_catalytic(his, acid, rieske_end):
    """H ... H ... D/E kisitlarini saglayan en korunmus ucluyu bul."""
    best, best_score = None, -1.0
    his_candidates = [i for i, v in enumerate(his)
                      if v >= MIN_CONSERVATION and i > rieske_end]
    acid_candidates = [i for i, v in enumerate(acid) if v >= MIN_CONSERVATION]
    for h1 in his_candidates:
        for h2 in his_candidates:
            sep = h2 - h1
            if not (CATALYTIC_HH[0] <= sep <= CATALYTIC_HH[1]):
                continue
            for d in acid_candidates:
                tail = d - h2
                if not (CATALYTIC_TAIL[0] <= tail <= CATALYTIC_TAIL[1]):
                    continue
                score = his[h1] + his[h2] + acid[d]
                if score > best_score:
                    best_score, best = score, (h1, h2, d)
    return best, best_score


def derive_columns(sto_path):
    """hmmalign ciktisindan motif kolonlarini cikar. Donen: dict."""
    aligned = read_stockholm_matchcols(sto_path)
    if not aligned:
        raise SystemExit(f"hizalama okunamadi: {sto_path}")
    cys = column_conservation(aligned, "C")
    his = column_conservation(aligned, "H")
    acid = column_conservation(aligned, "DE")

    rieske, r_score = pick_rieske(cys, his)
    if rieske is None:
        raise SystemExit("Rieske motifi bulunamadi: kisitlari gevsetmek gerekebilir")
    catalytic, c_score = pick_catalytic(his, acid, rieske[3])
    if catalytic is None:
        raise SystemExit("katalitik triad bulunamadi")

    # Kopru D/E: ilk katalitik His'ten hemen once, en korunmus asidik kolon
    window = range(max(0, catalytic[0] - 6), catalytic[0])
    bridging = max(window, key=lambda i: acid[i]) if window else catalytic[0] - 3

    one = lambda i: i + 1          # noqa: E731  (0-tabanli -> 1-tabanli)
    out = {
        "rieske": [[one(rieske[0]), "C", "Rieske Cys-1"], [one(rieske[1]), "H", "Rieske His-1"],
                   [one(rieske[2]), "C", "Rieske Cys-2"], [one(rieske[3]), "H", "Rieske His-2"]],
        "catalytic": [[one(catalytic[0]), "H", "mononukleer Fe His-1"],
                      [one(catalytic[1]), "H", "mononukleer Fe His-2"],
                      [one(catalytic[2]), "DE", "mononukleer Fe karboksilat"]],
        "bridging": [one(bridging), "DE", "alt-birimler arasi elektron transfer Asp"],
        "model_length": len(next(iter(aligned.values()))),
        "n_references": len(aligned),
        "conservation": {
            str(one(rieske[0])): round(cys[rieske[0]], 4),
            str(one(rieske[1])): round(his[rieske[1]], 4),
            str(one(rieske[2])): round(cys[rieske[2]], 4),
            str(one(rieske[3])): round(his[rieske[3]], 4),
            str(one(catalytic[0])): round(his[catalytic[0]], 4),
            str(one(catalytic[1])): round(his[catalytic[1]], 4),
            str(one(catalytic[2])): round(acid[catalytic[2]], 4),
            str(one(bridging)): round(acid[bridging], 4),
        },
    }
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--refs", default="ROs_71_Clean/refs71.fasta")
    ap.add_argument("--out-dir", default="ROs_71_Clean")
    ap.add_argument("--name", default="ROmotif71")
    ap.add_argument("--derive-only", action="store_true",
                    help="modeli yeniden kurma, var olan hizalamadan kolonlari cikar")
    ap.add_argument("--threads", type=int, default=8)
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    aln = os.path.join(args.out_dir, f"{args.name}_aln.sto")
    hmm = os.path.join(args.out_dir, f"{args.name}.hmm")
    hmmaln = os.path.join(args.out_dir, f"{args.name}_hmmaln.sto")

    if not args.derive_only:
        n = sum(1 for line in open(args.refs) if line.startswith(">"))
        print(f"[bilgi] {n} referans dizisi: {args.refs}")
        run(["clustalo", "-i", args.refs, "-o", aln, "--outfmt=st",
             "--threads", str(args.threads), "--force"])
        run(["hmmbuild", "--amino", "-n", args.name, hmm, aln])
        run(["hmmalign", "--trim", "--amino", "-o", hmmaln, hmm, args.refs])
    else:
        # Mevcut kurulumda hizalama baska adla duruyor olabilir
        for candidate in (hmmaln, os.path.join(args.out_dir, "refs71_hmmaln.sto")):
            if os.path.exists(candidate):
                hmmaln = candidate
                break
        else:
            raise SystemExit("hmmalign ciktisi bulunamadi; --derive-only kullanma")

    columns = derive_columns(hmmaln)
    out_path = os.path.join(args.out_dir, "motif_columns.json")
    payload = {args.name: columns}
    if os.path.exists(out_path):
        with open(out_path) as fh:
            payload = {**json.load(fh), args.name: columns}
    with open(out_path, "w") as fh:
        json.dump(payload, fh, indent=1)

    print(f"\n[cikarilan kolonlar] model {args.name}, uzunluk {columns['model_length']}, "
          f"{columns['n_references']} referans")
    for group in ("rieske", "catalytic"):
        for col, residues, label in columns[group]:
            print(f"   kolon {col:>4}  {residues:3s}  {label:34s} "
                  f"korunma {columns['conservation'][str(col)]:.3f}")
    col, residues, label = columns["bridging"]
    print(f"   kolon {col:>4}  {residues:3s}  {label:34s} "
          f"korunma {columns['conservation'][str(col)]:.3f}")
    print(f"[yazildi] {out_path}")
    print("\nro_motif.py bu dosyayi otomatik okur. Kolonlar beklenenden farkli "
          "ciktiysa referans setini ve hizalamayi gozden gecir; sayilar sessizce "
          "kullanilmaz, once burada basilir.")


if __name__ == "__main__":
    sys.exit(main())
