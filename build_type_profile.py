"""Yeni bir tip icin ATAMA PROFILI kurar ve kutuphaneye ekler.

NEDEN BU SCRIPT VAR -- ve neden olmamasi sessiz bir hataydi.

Yeni bir tip eklemenin zinciri su ana kadar su noktada kopuyordu:

    types/<id>.yaml           tek kaynak            vardi
    types_build.py            dogrulama + uretim    vardi
    build_reference_models.py motif modeli          vardi
    RieskeDB.hmm              ATAMA PROFILI         YOKTU

Atama, adayin kapsama barajini gectigi profiller arasinda en yuksek bit
skorlusunu secer. Kutuphanede profili olmayan bir tip hicbir yarismaya
GIREMEZ. Yani bugun biri yeni bir tip eklese, manifestoya girerdi,
chemistry.csv'de gorunurdu, motif modeli guncellenirdi -- ve tip sessizce
BOS kalirdi. Bu, bu hafta emekliye ayrilan dort bos tiple ayni basarisizlik
bicimi, yalnizca ters yonden: orada iki referans ayni diziyi paylasiyordu,
burada referansin hic profili olmuyor.

PROFIL NEDEN TEK DIZIDEN KURULMAZ. Kutuphanedeki mevcut profiller AILE
duzeyinde: her biri 2 ile 900 arasi diziden kurulmus. Tek bir diziden kurulan
profil daha dar ve daha dusuk skor uretir, yani yeni tip mevcut tiplere karsi
yarismayi sistematik olarak kaybeder -- gercekten daha iyi eslesse bile. Bu
yuzden burada da aile profili kurulur: referansin yakin homologlari zaten
elimizdeki aday havuzundan toplanir, hizalanir, ve profil ondan kurulur.

Profil buyuklugu ile kazanilan uye sayisi arasinda guclu bir iliski OLCULDU ve
bulunamadi (Spearman rho 0,23, p=0,06), yani buyuk profil tek basina avantaj
degil. Ama tek dizilik bir profil ile 900 dizilik bir profil arasindaki fark o
olcumun kapsadigi araligin disinda; amac esitlik degil, KARSILASTIRILABILIRLIK.

Kullanim:
    python3 build_type_profile.py --manifest types/6_601_CMO.yaml            # rapor
    python3 build_type_profile.py --manifest types/6_601_CMO.yaml --apply    # ekle
"""

import argparse
import os
import re
import shutil
import subprocess
import sys
import tempfile

LIBRARY = "RieskeDB71.hmm"
CANDIDATES = "ro_alpha.fasta"

# Homolog toplama esikleri.
#   MIN_IDENTITY  bunun altindaki eslesme ayni tipin uyesi sayilmaz; profil
#                 genisledikce komsu tiplerin alanina tasar.
#   MAX_HOMOLOGS  mevcut profillerin en buyugu 900 dizi; ustune cikmak
#                 kutuphaneyi dengesizlestirir.
#   MIN_HOMOLOGS  bu kadar bile bulunamazsa aile profili kurulamaz ve script
#                 tek dizilik profile DUSMEZ -- sessizce zayif bir profil
#                 uretmektense durur ve sebebini soyler.
MIN_IDENTITY = 55.0
MIN_COVERAGE = 0.70
MAX_HOMOLOGS = 400
MIN_HOMOLOGS = 5

CLUSTER_RE = re.compile(r"^(\d+)_(\d+)_([A-Za-z0-9_]+)$")


def load_manifest(path):
    text = open(path).read()
    try:
        import yaml
        data = yaml.safe_load(text)
        if isinstance(data, dict):
            return {k: ("" if v is None else str(v)) for k, v in data.items()}
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
        elif key:
            data[key] = (data[key] + line.strip()).strip()
    return data


def library_names(path):
    return [l.split()[1] for l in open(path) if l.startswith("NAME")]


def collect_homologs(seq, cluster, candidates, workdir):
    """Aday havuzundan referansin yakin homologlarini topla.

    diamond kullanilir cunku 73 bin diziye karsi hizli ve kimlik/kapsama
    dogrudan okunur. Esiklerin altinda kalanlar alinmaz: profil komsu tiplerin
    alanina tasarsa atamayi bozar, ki duzeltmeye calistigimiz sey tam budur.
    """
    q = os.path.join(workdir, "ref.fasta")
    with open(q, "w") as fh:
        fh.write(f">{cluster}\n{seq}\n")
    db = os.path.join(workdir, "cand")
    out = os.path.join(workdir, "hits.tsv")
    try:
        subprocess.run(["diamond", "makedb", "--in", candidates, "-d", db,
                        "--quiet"], check=True, capture_output=True)
        subprocess.run(["diamond", "blastp", "-q", q, "-d", db, "-o", out,
                        "--quiet", "--max-target-seqs", str(MAX_HOMOLOGS * 3),
                        "--outfmt", "6", "qseqid", "sseqid", "pident",
                        "length", "qlen", "slen"],
                       check=True, capture_output=True)
    except (subprocess.CalledProcessError, FileNotFoundError) as exc:
        return None, f"diamond calistirilamadi: {exc}"

    keep, seen = [], set()
    import csv as _csv
    for row in _csv.reader(open(out), delimiter="\t"):
        _, sid, pid, ln, ql, sl = (row[0], row[1], float(row[2]), int(row[3]),
                                   int(row[4]), int(row[5]))
        if pid < MIN_IDENTITY or ln / min(ql, sl) < MIN_COVERAGE:
            continue
        if sid in seen:
            continue
        seen.add(sid)
        keep.append((pid, sid))
    keep.sort(reverse=True)
    return keep[:MAX_HOMOLOGS], None


def write_homolog_fasta(names, candidates, seq, cluster, path):
    """Secilen homologlarin dizilerini aday fastasindan cikar."""
    wanted = {n for _, n in names}
    found, cur, keep = {}, None, False
    for line in open(candidates):
        if line.startswith(">"):
            cur = line[1:].split()[0]
            keep = cur in wanted
            if keep:
                found[cur] = []
        elif keep:
            found[cur].append(line.strip())
    with open(path, "w") as fh:
        fh.write(f">{cluster}_reference\n{seq}\n")
        for name, parts in found.items():
            fh.write(f">{name}\n{''.join(parts)}\n")
    return len(found)


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--library", default=LIBRARY)
    ap.add_argument("--candidates", default=CANDIDATES,
                    help="homologlarin toplanacagi havuz")
    ap.add_argument("--apply", action="store_true",
                    help="profili kutuphaneye EKLE. Varsayilan rapordur")
    ap.add_argument("--cpu", type=int, default=8)
    args = ap.parse_args()

    m = load_manifest(args.manifest)
    cluster = (m.get("cluster") or "").strip()
    seq = (m.get("sequence") or "").strip().upper().replace(" ", "")
    print(f"=== {cluster} ===")
    if not CLUSTER_RE.match(cluster):
        raise SystemExit(f"kume adi kalibi tutmuyor: {cluster}")
    if not seq:
        raise SystemExit("manifestoda dizi yok")
    print(f"dizi: {len(seq)} aa")

    names = library_names(args.library)
    if cluster in names:
        print(f"\n{cluster} ZATEN kutuphanede. Once cikarilmali "
              f"(hmmfetch ile yeniden kur), yoksa iki profil ayni adi tasir.")
        return 1
    print(f"kutuphane: {len(names)} profil")

    work = tempfile.mkdtemp(prefix="typeprof_")
    try:
        hits, err = collect_homologs(seq, cluster, args.candidates, work)
        if err:
            print(f"  {err}")
            return 1
        print()
        print(f"homolog taramasi (kimlik >={MIN_IDENTITY:.0f} %, "
              f"kapsama >={MIN_COVERAGE:.2f}): {len(hits)} dizi")
        if hits:
            print(f"  en yakin : {hits[0][0]:.1f} %")
            print(f"  en uzak  : {hits[-1][0]:.1f} %")
        if len(hits) < MIN_HOMOLOGS:
            print()
            print(f"  DURDU: {MIN_HOMOLOGS} homolog bulunamadi.")
            print("  Tek dizilik bir profile DUSULMEZ: boyle bir profil mevcut")
            print("  aile profillerine karsi sistematik olarak dusuk skorlar ve")
            print("  tip, gercekten daha iyi eslesse bile bos kalir. Bu referans")
            print("  bu havuzda temsil edilmiyor demektir -- havuzu genislet ya")
            print("  da tipin gercekten yeni bir dal oldugunu kabul et.")
            return 1

        hf = os.path.join(work, "family.fasta")
        n = write_homolog_fasta(hits, args.candidates, seq, cluster, hf)
        print(f"  profil {n + 1} diziden kurulacak (referans dahil)")

        aln = os.path.join(work, "family.sto")
        hmm = os.path.join(work, "family.hmm")
        try:
            subprocess.run(["clustalo", "-i", hf, "-o", aln, "--outfmt=st",
                            "--threads", str(args.cpu), "--force"],
                           check=True, capture_output=True)
            subprocess.run(["hmmbuild", "--amino", "-n", cluster, hmm, aln],
                           check=True, capture_output=True)
        except (subprocess.CalledProcessError, FileNotFoundError) as exc:
            print(f"  profil kurulamadi: {exc}")
            return 1
        leng = next((l.split()[1] for l in open(hmm) if l.startswith("LENG")), "?")
        nseq = next((l.split()[1] for l in open(hmm) if l.startswith("NSEQ")), "?")
        print(f"  profil hazir: LENG={leng}, NSEQ={nseq}")

        print()
        if not args.apply:
            print("Yalnizca rapor. Eklemek icin --apply.")
            print("Ekledikten SONRA atama yeniden kosmali, yoksa yeni tip")
            print("veritabaninda gorunmez: python3 rebuild.py --run --include-heavy")
            return 0

        shutil.copy2(args.library, args.library + ".bak")
        with open(args.library, "a") as dst, open(hmm) as src:
            dst.write(src.read())
        subprocess.run(["hmmpress", "-f", args.library], check=True,
                       capture_output=True)
        after = len(library_names(args.library))
        print(f"[eklendi] {args.library}: {len(names)} -> {after} profil")
        print(f"[yedek]   {args.library}.bak")
        print()
        print("Simdi: python3 rebuild.py --run --include-heavy")
        return 0
    finally:
        shutil.rmtree(work, ignore_errors=True)


if __name__ == "__main__":
    sys.exit(main())
