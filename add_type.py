"""Yeni bir RO tipi ekle -- once NE OLACAGINI goster, sonra yaz.

NEDEN BU ARAC VAR

Yeni Rieske oksijenazlar kesfedilmeye devam edecek ve bu sete eklenebilmeli.
Bugun eklemek tek bir dosya degistirmek degil: dizi `refs<N>.fasta`'ya,
kimya `chemistry.csv`'ye gidiyor, hizalama ve profil kutuphanesi yeniden
kuruluyor, ve ikisinin BIRBIRIYLE tutarli oldugunu hicbir sey dogrulamiyor.
Kimin, ne zaman, hangi makaleden ekledigi de hicbir yerde yazmiyor.

Ama asil sorun bu degil. ASIL SORUN SU: atama YARISMALI. Bir protein,
kapsama barajini gectigi profiller arasinda EN YUKSEK BIT SKORLU olani alir
(`ro_filter.classify_target`). Yani yeni bir referans eklemek toplama bir sey
EKLEMEZ -- mevcut tiplerden uye CALAR, ve bazen yeni tip kendisi bos kalir.

Bu varsayim degil, olculmus gecmis: bu veritabaninda 10 tipin hic uyesi yok ve
hicbiri kapsama barajina takilmiyor; hepsi yarismayi kaybediyor. Dordu ise bir
referansin BIREBIR KOPYASI (CarAa=CARDO, NahAc=NDO 3_315, NdmC=NdmB,
OxoO=OMO). Ayni diziye sahip iki referans yarismayi beraberlikle bitirir ve
biri keyfi olarak her seyi alir.

Bu yuzden aracin varsayilan davranisi YAZMAMAKTIR. Once sunu soyler:

    bu tip kac aday kazanir
    hangi tiplerden kac uye calar
    kendisi bos kalacak mi
    mevcut bir referansla cok mu benziyor

ve ancak --apply verilirse dosyalara dokunur.

UCUZ YOL. Pahali adim 17 bin genomdan baglam cikarmaktir; atama ise
`ro_alpha.fasta` icindeki 73.272 hazir aday proteine karsi kosar. Bu arac
yalnizca ikincisini yapar: yeni referans icin tek profil kurar, aday setine
karsi `hmmsearch` kosturur ve yarismayi mevcut atamalarla karsilastirir.
Dakikalar surer, gunler degil.

BILDIRIM DOSYASI (types/<cluster>.yaml ya da .json) -- tipin TEK kaynagi:

    cluster:          5_509_AbcD          # grup_numara_gen
    gene:             AbcD
    group:            5
    sequence:         MTSS...             # ya da sequence_file: yol
    substrate_en:     benzalkonium
    substrate_smiles: "CCCC..."
    product_en:       ...
    reaction:         ...
    reaction_class:   hydroxylation
    family:           quaternary_amines
    pdb:              7ABC                # yoksa bos
    source:           "doi:10.1234/xyz"
    source_kind:      paper               # paper | structure | paper+structure
    curation_confidence: high             # high | medium | low
    added_by:         "ad soyad"
    added_on:         2026-10-07
    notes:            ...

Kullanim:
    python3 add_type.py types/5_509_AbcD.yaml              # yalnizca rapor
    python3 add_type.py types/5_509_AbcD.yaml --apply      # yaz
    python3 add_type.py --export-manifests types/          # mevcut 71'i disa aktar
"""

import argparse
import csv
import json
import os
import re
import shutil
import subprocess
import sqlite3
import sys
import tempfile
from collections import Counter, defaultdict

# Tekrar esigi. 24. maddede olculen gercek vakalar %99,2 ile %100 arasinda;
# esik oraya gore konuldu. Altinda kalan ama yine de yuksek olan ciftler
# (ornegin EdoA1/cumA1 %99,8) uyari uretir, reddedilmez -- deneysel olarak
# farkli substratlari olabilir.
DUPLICATE_REJECT = 99.0
DUPLICATE_WARN = 90.0

CLUSTER_RE = re.compile(r"^(\d+)_(\d+)_([A-Za-z0-9_]+)$")

CHEMISTRY_FIELDS = ["cluster", "substrate_en", "substrate_smiles", "product_en",
                    "reaction", "reaction_class", "family", "pdb", "source",
                    "source_kind", "curation_confidence", "notes"]

REQUIRED = ["cluster", "group", "substrate_en", "product_en", "reaction",
            "reaction_class", "source", "source_kind"]


# --------------------------------------------------------------- bildirim
def load_manifest(path):
    """YAML ya da JSON. PyYAML yoksa kucuk bir ayristirici devreye girer.

    Bagimlilik eklememek icin: bildirim dosyasi duz `anahtar: deger`
    satirlarindan olusuyor, ic ice yapi yok. Bu kadarini okumak icin PyYAML
    sart degil ve sart kosmak araci kurulum sorununa cevirirdi.
    """
    text = open(path).read()
    if path.endswith(".json"):
        return json.loads(text)
    try:
        import yaml
        return yaml.safe_load(text)
    except ImportError:
        data, key = {}, None
        for line in text.splitlines():
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            if line[0] not in " \t" and ":" in line:
                key, _, value = line.partition(":")
                key = key.strip()
                data[key] = value.strip().strip('"').strip("'")
            elif key:                      # devam satiri (katlanmis dizi)
                data[key] = (data[key] + line.strip()).strip()
        return data


def validate(manifest, chemistry, refs):
    """Bildirimi dogrula. Hatalar listesi doner; bos liste = gecerli."""
    errors, warnings = [], []
    cluster = (manifest.get("cluster") or "").strip()

    if not cluster:
        errors.append("cluster alani bos")
    elif not CLUSTER_RE.match(cluster):
        errors.append(f"cluster adi '{cluster}' kalibi tutmuyor: grup_numara_gen, "
                      "ornegin 5_509_AbcD")
    elif cluster in chemistry:
        errors.append(f"'{cluster}' zaten chemistry.csv'de var")
    elif cluster in refs:
        errors.append(f"'{cluster}' zaten referans setinde var")

    m = CLUSTER_RE.match(cluster) if cluster else None
    if m and manifest.get("group") and str(manifest["group"]) != m.group(1):
        errors.append(f"group alani {manifest['group']} ama cluster adi "
                      f"{m.group(1)} diyor -- ikisi ayni olmali")

    for field in REQUIRED:
        if not str(manifest.get(field) or "").strip():
            errors.append(f"zorunlu alan bos: {field}")

    seq = (manifest.get("sequence") or "").strip().replace(" ", "").upper()
    if not seq and manifest.get("sequence_file"):
        path = manifest["sequence_file"]
        if os.path.exists(path):
            seq = "".join(l.strip() for l in open(path)
                          if not l.startswith(">")).upper()
        else:
            errors.append(f"sequence_file bulunamadi: {path}")
    if not seq:
        errors.append("dizi yok (sequence ya da sequence_file gerekli)")
    else:
        bad = set(seq) - set("ACDEFGHIKLMNPQRSTVWY")
        if bad:
            errors.append(f"dizide amino asit olmayan karakterler: "
                          f"{''.join(sorted(bad))}")
        if len(seq) < 250:
            warnings.append(f"dizi kisa ({len(seq)} aa). RO alfa alt birimleri "
                            "tipik olarak 330-470 aa; bu bir parca olabilir")
        if len(seq) > 700:
            warnings.append(f"dizi uzun ({len(seq)} aa) -- kaynasmis bir "
                            "protein mi?")

    kind = (manifest.get("source_kind") or "").strip()
    if kind and kind not in ("paper", "structure", "paper+structure", "general"):
        warnings.append(f"source_kind '{kind}' alisilmis degerlerden biri degil")
    src = (manifest.get("source") or "")
    if not re.search(r"10\.\d{4,}/|PMID|PMC|PDB|\b\d{4};", src):
        warnings.append("source alaninda cozulebilir bir tanimlayici yok "
                        "(DOI, PMID ya da PDB). 71 tipin 50'si bu yuzden "
                        "takip edilemez durumda -- yenisini oyle ekleme")

    return errors, warnings, seq


# ---------------------------------------------------------------- tekrar
def duplicate_check(seq, refs, workdir):
    """Yeni dizi mevcut referanslara ne kadar benziyor?

    Bu tek kontrol, bugunku bos tiplerin dordunu bastan onlerdi.
    """
    qpath = os.path.join(workdir, "new.fasta")
    with open(qpath, "w") as fh:
        fh.write(">new\n" + seq + "\n")
    dbpath = os.path.join(workdir, "refs")
    refpath = os.path.join(workdir, "refs.fasta")
    with open(refpath, "w") as fh:
        for name, s in refs.items():
            fh.write(f">{name}\n{s}\n")
    try:
        subprocess.run(["diamond", "makedb", "--in", refpath, "-d", dbpath,
                        "--quiet"], check=True, capture_output=True)
        out = os.path.join(workdir, "hits.tsv")
        subprocess.run(["diamond", "blastp", "-q", qpath, "-d", dbpath,
                        "-o", out, "--quiet", "--ultra-sensitive",
                        "--max-target-seqs", "20", "--outfmt", "6",
                        "qseqid", "sseqid", "pident", "length", "qlen", "slen"],
                       check=True, capture_output=True)
    except (subprocess.CalledProcessError, FileNotFoundError) as exc:
        return None, f"diamond calistirilamadi ({exc}); tekrar kontrolu ATLANDI"

    best = {}
    for row in csv.reader(open(out), delimiter="\t"):
        _, sid, pid, ln, ql, sl = row[0], row[1], float(row[2]), int(row[3]), \
            int(row[4]), int(row[5])
        cov = ln / min(ql, sl)
        if cov < 0.8:
            continue
        if pid > best.get(sid, (0, 0))[0]:
            best[sid] = (pid, cov)
    hits = sorted(((p, c, n) for n, (p, c) in best.items()), reverse=True)
    return hits, None


# ------------------------------------------------------- yarisma onizleme
def current_winners(domtbl):
    """Her adayin BUGUNKU kazanani ve skoru.

    Ilk surum aday basligindaki `cluster=` etiketine bakiyordu ve yalnizca
    "bu adaylar su an kimde" diyebiliyordu -- yani bir UST SINIR. Oysa mevcut
    skorlar zaten diskte: `cand_dom.out` her aday x profil satirini tasiyor.
    Oradan okununca onizleme ust sinir olmaktan cikip KESIN cevaba doner:
    yeni profil bir adayi ancak mevcut kazanandan daha yuksek skorlarsa alir.
    """
    best = {}
    if not os.path.exists(domtbl):
        return None
    with open(domtbl) as fh:
        for line in fh:
            if line[0] == "#":
                continue
            p = line.split()
            target, query, score = p[0], p[3], float(p[13])
            if score > best.get(target, (None, -1e9))[1]:
                best[target] = (query, score)
    return best


def competition_preview(seq, cluster, candidates_fasta, workdir, cpu=4):
    """Yeni profil aday setine karsi kosar."""
    hmm = os.path.join(workdir, "new.hmm")
    seqfile = os.path.join(workdir, "new_ref.fasta")
    with open(seqfile, "w") as fh:
        fh.write(f">{cluster}\n{seq}\n")
    try:
        subprocess.run(["hmmbuild", "--amino", "-n", cluster, hmm, seqfile],
                       check=True, capture_output=True)
    except (subprocess.CalledProcessError, FileNotFoundError) as exc:
        return None, f"hmmbuild calistirilamadi ({exc})"

    tbl = os.path.join(workdir, "new.domtbl")
    try:
        subprocess.run(["hmmsearch", "--domtblout", tbl, "--cpu", str(cpu),
                        "-E", "1e-5", hmm, candidates_fasta],
                       check=True, capture_output=True)
    except (subprocess.CalledProcessError, FileNotFoundError) as exc:
        return None, f"hmmsearch calistirilamadi ({exc})"

    n_candidates = sum(1 for l in open(candidates_fasta) if l.startswith(">"))

    # yeni profilin her aday icin en iyi skoru
    new_score, new_cov = {}, {}
    with open(tbl) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            p = line.split()
            target, tlen, score = p[0], int(p[2]), float(p[13])
            hmm_from, hmm_to, qlen = int(p[15]), int(p[16]), int(p[5])
            cov = (hmm_to - hmm_from + 1) / qlen
            if score > new_score.get(target, -1e9):
                new_score[target] = score
                new_cov[target] = cov
    return {"n_candidates": n_candidates, "new_score": new_score,
            "new_cov": new_cov}, None


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("manifest", nargs="?", help="types/<cluster>.yaml")
    ap.add_argument("--db", default="roar.sqlite")
    ap.add_argument("--chemistry", default="chemistry.csv")
    ap.add_argument("--refs", default="ROs_71_Clean/refs71.fasta")
    ap.add_argument("--candidates", default="genomic_context/ro_candidates.fasta",
                    help="atamanin kostugu aday seti -- genomlardan cikarilan "
                         "proteinler, UniProt taramasi degil")
    ap.add_argument("--current-domtbl", default="genomic_context/cand_dom.out",
                    help="mevcut yarismanin sonucu; her adayin bugunku en "
                         "yuksek skoru buradan okunur")
    ap.add_argument("--min-coverage", type=float, default=0.45,
                    help="ro_filter ile ayni baraj")
    ap.add_argument("--cpu", type=int, default=4)
    ap.add_argument("--apply", action="store_true",
                    help="dosyalara YAZ. Varsayilan yalnizca rapordur")
    ap.add_argument("--export-manifests", metavar="DIR",
                    help="mevcut 71 tipi bildirim dosyasi olarak disa aktar")
    args = ap.parse_args()

    chemistry = {}
    with open(args.chemistry) as fh:
        for row in csv.DictReader(fh):
            chemistry[row["cluster"]] = row

    refs = {}
    name = None
    for line in open(args.refs):
        if line.startswith(">"):
            name = line[1:].strip().split()[0]
            refs[name] = ""
        elif name:
            refs[name] += line.strip()
    # referans basligi "1_101_OxoO_monooxygenase_pro" -> kume kimligi
    ref_by_cluster = {}
    for header, s in refs.items():
        mm = re.match(r"^(\d+_\d+_[A-Za-z0-9]+)", header)
        ref_by_cluster[mm.group(1) if mm else header] = s

    if args.export_manifests:
        return export_manifests(args.export_manifests, chemistry, ref_by_cluster)

    if not args.manifest:
        ap.error("bildirim dosyasi gerekli (ya da --export-manifests)")

    manifest = load_manifest(args.manifest)
    errors, warnings, seq = validate(manifest, chemistry, ref_by_cluster)
    cluster = (manifest.get("cluster") or "?").strip()

    print(f"=== {cluster} ===")
    print(f"bildirim: {args.manifest}")
    if seq:
        print(f"dizi:     {len(seq)} aa")
    print()

    if errors:
        print("GECERSIZ:")
        for e in errors:
            print(f"  - {e}")
        print("\nHicbir sey yazilmadi.")
        return 1
    for w in warnings:
        print(f"  [uyari] {w}")
    if warnings:
        print()

    workdir = tempfile.mkdtemp(prefix="addtype_")
    try:
        # 1. tekrar
        hits, err = duplicate_check(seq, ref_by_cluster, workdir)
        blocked = False
        if err:
            print(f"  [uyari] {err}")
        elif hits:
            top = hits[:5]
            print("En benzer mevcut referanslar:")
            for pid, cov, nm in top:
                flag = ("  <-- KOPYA" if pid >= DUPLICATE_REJECT else
                        "  <-- cok benzer" if pid >= DUPLICATE_WARN else "")
                print(f"  {pid:6.1f} %  kapsama {cov:.2f}  {nm}{flag}")
            if top[0][0] >= DUPLICATE_REJECT:
                blocked = True
                print()
                print(f"  REDDEDILDI: {top[0][2]} ile %{top[0][0]:.1f} ayni.")
                print("  Ayni diziye sahip iki referans yarismayi beraberlikle")
                print("  bitirir ve biri keyfi olarak butun uyeleri alir --")
                print("  bu veritabanindaki 10 bos tipin dordu tam olarak boyle")
                print("  olustu. Es anlamli ad olarak eklemek istiyorsan")
                print("  chemistry.csv'deki mevcut satira not dus.")
        else:
            print("Mevcut referanslarla kayda deger benzerlik yok.")
        print()

        # 2. yarisma
        if not os.path.exists(args.candidates):
            print(f"  [atlandi] aday seti yok: {args.candidates}")
        else:
            res, err = competition_preview(seq, cluster, args.candidates,
                                           workdir, args.cpu)
            if err:
                print(f"  [uyari] {err}")
            else:
                report_competition(res, args.min_coverage, cluster,
                                   current_winners(args.current_domtbl))
    finally:
        shutil.rmtree(workdir, ignore_errors=True)

    print()
    if blocked:
        print("Hicbir sey yazilmadi (tekrar).")
        return 1
    if not args.apply:
        print("Yalnizca rapor. Yazmak icin --apply ekle.")
        print("Yazilacak olanlar:")
        print(f"  {args.chemistry}  <- yeni satir {cluster}")
        print(f"  {args.refs}       <- yeni dizi")
        print("  sonra: python3 build_reference_models.py  (hizalama, profil,")
        print("         motif kolonlari) ve ardindan atama adimi yeniden")
        return 0

    apply_changes(args, manifest, seq, cluster)
    return 0


def report_competition(res, min_cov, cluster, winners):
    """Yeni tip kimden kac uye ALIR, ve kendisi ayakta kalir mi."""
    new_score, new_cov = res["new_score"], res["new_cov"]
    print(f"Aday seti: {res['n_candidates']:,} protein "
          f"(genomlardan cikarilanlar)")
    eligible = [t for t, c in new_cov.items() if c >= min_cov]
    print(f"Yeni profilin kapsama barajini ({min_cov}) gectigi: {len(eligible):,}")
    if not eligible:
        print("  Bu tip hicbir adayi tutmuyor -- eklenirse BOS kalir.")
        return
    if winners is None:
        print("  [uyari] mevcut skor tablosu yok; kesin yarisma hesaplanamadi")
        return

    taken = Counter()
    margins = defaultdict(list)
    unclaimed = 0
    for t in eligible:
        cur = winners.get(t)
        if cur is None:
            unclaimed += 1
            taken["(bugun hicbir tipte degil)"] += 1
            continue
        owner, owner_score = cur
        if new_score[t] > owner_score:
            taken[owner] += 1
            margins[owner].append(new_score[t] - owner_score)

    total = sum(taken.values())
    print()
    if not total:
        print("  Bu tip TEK BIR adayi bile kazanmiyor: kapsama barajini gectigi")
        print("  her aday, mevcut sahibinde daha yuksek skorluyor. Eklenirse")
        print("  BOS KALIR -- bu veritabanindaki 10 bos tip boyle olustu.")
        return

    print(f"Yeni tip {total:,} aday KAZANIR, su tiplerden:")
    for owner, n in taken.most_common(15):
        if owner.startswith("("):
            print(f"  {n:6,}  {owner}")
            continue
        med = sorted(margins[owner])[len(margins[owner]) // 2]
        print(f"  {n:6,}  {owner:<18} (ortanca fark {med:+.1f} bit)")
    if len(taken) > 15:
        print(f"  ... ve {len(taken) - 15} tip daha")
    print()
    print("  Bu KESIN sonuc: yeni profilin skoru ile her adayin bugunku en")
    print("  yuksek skoru dogrudan karsilastirildi. Kucuk bit farklariyla")
    print("  kazanilan uyeler kirilgandir -- referans seti bir daha")
    print("  degistiginde geri gidebilir.")


def apply_changes(args, manifest, seq, cluster):
    """Dosyalari yaz. Her birinin yedegi alinir."""
    for path in (args.chemistry, args.refs):
        shutil.copy2(path, path + ".bak")
    row = {f: str(manifest.get(f, "") or "") for f in CHEMISTRY_FIELDS}
    row["cluster"] = cluster
    with open(args.chemistry) as fh:
        rows = list(csv.DictReader(fh))
        fields = list(rows[0].keys()) if rows else CHEMISTRY_FIELDS
    rows.append({f: row.get(f, "") for f in fields})
    rows.sort(key=lambda r: r["cluster"])
    with open(args.chemistry, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)
    gene = manifest.get("gene") or cluster.split("_")[-1]
    with open(args.refs, "a") as fh:
        fh.write(f">{cluster}_{gene}_pro\n{seq}\n")
    print(f"[yazildi] {args.chemistry} (+1 satir), {args.refs} (+1 dizi)")
    print(f"[yedek]   {args.chemistry}.bak, {args.refs}.bak")
    print()
    print("Simdi sirayla:")
    print("  python3 build_reference_models.py      # hizalama, profil, motifler")
    print("  python3 run_ro_filter.py               # atama (yarisma yeniden)")
    print("  bash run_all.sh                        # geri kalan pipeline")
    print()
    print("Atama yeniden kosmadan sitedeki sayilar DEGISMEZ ve yeni tip bos")
    print("gorunur.")


def export_manifests(outdir, chemistry, refs):
    """Mevcut 71 tipi bildirim dosyasina cevirir -- gecis icin."""
    os.makedirs(outdir, exist_ok=True)
    written = missing_seq = 0
    for cluster, row in sorted(chemistry.items()):
        seq = refs.get(cluster, "")
        if not seq:
            missing_seq += 1
        m = CLUSTER_RE.match(cluster)
        lines = [f"cluster: {cluster}",
                 f"gene: {m.group(3) if m else ''}",
                 f"group: {m.group(1) if m else ''}"]
        for f in CHEMISTRY_FIELDS[1:]:
            v = (row.get(f) or "").replace("\n", " ")
            lines.append(f'{f}: "{v}"' if ('"' not in v and (":" in v or "#" in v))
                         else f"{f}: {v}")
        lines.append("added_by: (migrated from chemistry.csv)")
        lines.append("added_on: 2026-10-07")
        lines.append(f"sequence: {seq}" if seq else "sequence:   # EKSIK")
        with open(os.path.join(outdir, cluster + ".yaml"), "w") as fh:
            fh.write("\n".join(lines) + "\n")
        written += 1
    print(f"[yazildi] {written} bildirim -> {outdir}")
    if missing_seq:
        print(f"  {missing_seq} tipin referans dizisi bulunamadi "
              "(bildirimde sequence alani bos birakildi)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
