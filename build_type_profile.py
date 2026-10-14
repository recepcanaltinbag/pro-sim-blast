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

# ESIK SABIT DEGIL, REFERANSA GORE. Nedeni olculdu.
#
# Ilk surum mutlak bir %55 esigi kullaniyordu ve bu, kullanicinin isaret ettigi
# hatayi uretiyordu: referansin havuzdaki en yakin komsusu %60'taysa, %55 esigi
# referanstan %40 UZAK dizileri toplar ve kurulan profil o enzimi degil,
# referans ile etrafta ne varsa onun arasindaki bulanik bir bulutu temsil eder.
# Profil referanstan kayar, komsu tiplerin alanina tasar, ve referansin kendisi
# kendi profiline en iyi eslesme olmayabilir.
#
# 67 referansin havuzdaki EN YAKIN komsusu olculdu:
#     >=%95   37 referans   aile profili saglam
#     %80-95  15 referans   saglam
#     %65-80   8 referans   KAYMA BOLGESI
#     <%65     7 referans   gercekten izole (CndA %42,6, qxyA %45,8, Dpe %50,6)
#
# Yani sorun nadir degil: referanslarin dortte biri mutlak esikle yanlis
# profil uretirdi. Pencere artik EN YAKIN KOMSUYA gore acilir.
IDENTITY_WINDOW = 20.0    # en yakin komsudan bu kadar asagiya kadar al
ABSOLUTE_FLOOR = 50.0     # ne olursa olsun bunun altina inme
MIN_COVERAGE = 0.70
MAX_HOMOLOGS = 400
MIN_HOMOLOGS = 5

# Rejim sinirlari: en yakin komsunun kimligi hangi davranisi gerektirir.
TIGHT_ENOUGH = 85.0       # ustu: aile profili sorunsuz
DRIFT_ZONE = 65.0         # arasi: aile profili kurulur ama kayma sinanir
                          # alti: izole; tek dizilik profil DURUST olandir

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

    hits, seen = [], set()
    import csv as _csv
    for row in _csv.reader(open(out), delimiter="\t"):
        _, sid, pid, ln, ql, sl = (row[0], row[1], float(row[2]), int(row[3]),
                                   int(row[4]), int(row[5]))
        if ln / min(ql, sl) < MIN_COVERAGE or sid in seen:
            continue
        seen.add(sid)
        hits.append((pid, sid))
    hits.sort(reverse=True)
    if not hits:
        return [], None
    # Pencere en yakin komsuya gore: mutlak esik, referansin nerede durduguna
    # bakmadan uzak dizileri iceri alir ve profili kaydirir.
    top = hits[0][0]
    cutoff = max(ABSOLUTE_FLOOR, top - IDENTITY_WINDOW)
    keep = [h for h in hits if h[0] >= cutoff]
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


def drift_check(seq, cluster, hmm, family_fasta, workdir):
    """Kurulan profil hala REFERANSI en iyi puanliyor mu?

    Aile genisledikce profil referanstan uzaklasabilir. O noktada tip, adini
    tasidigi enzimi artik tanimlamaz -- kullanicinin isaret ettigi hata tam
    olarak budur. Profil kendi aile dosyasina karsi kosturulur ve referansin
    skoru en yuksek skorla karsilastirilir.
    """
    out = os.path.join(workdir, "drift.domtbl")
    try:
        subprocess.run(["hmmsearch", "--domtblout", out, "--noali",
                        "-E", "100", hmm, family_fasta],
                       check=True, capture_output=True)
    except (subprocess.CalledProcessError, FileNotFoundError):
        return None
    ref_score, best_score, best_name = None, -1e9, None
    with open(out) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            p = line.split()
            if len(p) < 14:
                continue
            target, score = p[0], float(p[13])
            if target.endswith("_reference"):
                ref_score = score
            if score > best_score:
                best_score, best_name = score, target
    if ref_score is None:
        return None
    return ref_score, best_score, (None if best_name.endswith("_reference")
                                   else best_name)


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
        top = hits[0][0] if hits else 0.0
        print()
        print(f"havuzdaki en yakin komsu: {top:.1f} %")
        if top >= TIGHT_ENOUGH:
            regime = "tight"
            print("  -> SIKI aile. Aile profili saglam.")
        elif top >= DRIFT_ZONE:
            regime = "drift"
            print("  -> KAYMA BOLGESI. Aile profili kurulur ama referanstan")
            print("     kayabilir; asagida kayma sinamasi var.")
        else:
            regime = "isolated"
            print("  -> IZOLE. Bu referansin havuzda yakin akrabasi YOK.")
            print("     Aile profili kurmak, enzimi degil etraftaki uzak")
            print("     dizileri modellemek olur. Tek dizilik profil kurulur:")
            print("     dar ama DURUST. Bu tipin az uye kazanmasi beklenir ve")
            print("     bu bir kusur degil, dogru cevaptir.")

        hf = os.path.join(work, "family.fasta")
        if regime == "isolated" or len(hits) < MIN_HOMOLOGS:
            if regime != "isolated":
                print()
                print(f"  {MIN_HOMOLOGS} homolog bulunamadi ({len(hits)} var);")
                print("  tek dizilik profile dusuluyor.")
            with open(hf, "w") as fh:
                fh.write(f">{cluster}_reference\n{seq}\n")
            n = 0
            single = True
        else:
            n = write_homolog_fasta(hits, args.candidates, seq, cluster, hf)
            single = False
            print(f"  pencere: {hits[-1][0]:.1f} % - {top:.1f} %  "
                  f"({len(hits)} dizi)")
            print(f"  profil {n + 1} diziden kurulacak (referans dahil)")

        aln = os.path.join(work, "family.sto")
        hmm = os.path.join(work, "family.hmm")
        try:
            if single:
                # Tek dizi hizalanamaz; clustalo bos olmayan bir hizalama
                # bekler ve hata verir. hmmbuild tek diziyi dogrudan kabul
                # eder ve onu bir dizilik hizalama sayar.
                subprocess.run(["hmmbuild", "--amino", "-n", cluster, hmm, hf],
                               check=True, capture_output=True)
            else:
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

        # KAYMA SINAMASI. Profil referansi hala en iyi temsil ediyor mu?
        # Aile genisledikce profil referanstan uzaklasabilir; o zaman
        # kuratorlu enzim kendi profiline en iyi eslesme OLMAZ ve tip, adini
        # tasidigi enzimi artik tanimlamaz.
        if not single:
            drift = drift_check(seq, cluster, hmm, hf, work)
            print()
            if drift is None:
                print("  [uyari] kayma sinamasi kosulamadi")
            else:
                ref_s, best_s, best_name = drift
                print(f"  kayma sinamasi: referans {ref_s:.1f} bit, "
                      f"profilin en iyi puanladigi dizi {best_s:.1f} bit")
                if best_name and ref_s < best_s * 0.90:
                    print("  UYARI: profil referansi en iyi dizi olarak")
                    print("  GORMUYOR -- aile referanstan kaymis demektir.")
                    print("  Pencereyi daraltmak (IDENTITY_WINDOW) ya da bu")
                    print("  tipi izole kabul etmek gerekir.")
                else:
                    print("  referans profilin merkezinde: kayma yok.")

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
