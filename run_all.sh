#!/usr/bin/env bash
#
# ROAR-DB pipeline'ini uctan uca calistirir.
#
# Kullanim:  bash run_all.sh [cpu]
#   cpu: hmmsearch/hmmalign icin cekirdek sayisi (varsayilan 8)
#
# Her adim bir oncekinin ciktisina bagli. Uzun adimlar (hmmsearch, extraction)
# ciktilari varsa atlar -- yeniden calistirmak guvenli.
#
# Not: adim 1 (189k protein filtresi) ve adim 2 (17k gbk cikarimi) uzundur
# (~30-60 dk toplam, CPU'ya bagli). Digerleri dakikalar surer.

set -euo pipefail
cd "$(dirname "$0")"

CPU="${1:-8}"
OUT="analysis_out"
GC="genomic_context"

step() { echo -e "\n\033[1;36m==> $*\033[0m"; }

# --- Gereksinim kontrolu
for tool in hmmsearch hmmalign hmmbuild hmmpress cd-hit python3; do
  command -v "$tool" >/dev/null 2>&1 || { echo "EKSIK: $tool"; exit 1; }
done
python3 -c "import Bio, pandas, numpy" || { echo "EKSIK python paketi (biopython/pandas/numpy)"; exit 1; }

# --- Temiz referans modeller (yoksa kur)
if [ ! -f "ROs_71_Clean/ROmotif71.hmm" ] || [ ! -f "RieskeDB71.hmm" ]; then
  step "Temiz referans modeller kuruluyor (71 referans)"
  echo "  HATA: ROs_71_Clean/ROmotif71.hmm ve RieskeDB71.hmm hazir olmali."
  echo "  Bunlar 72 referanstan IsoMO+CdnD cikarilarak uretildi (bkz. README)."
  exit 1
fi

# --- 1. Protein duzeyi filtre (UZUN)
step "1/14  Protein filtresi (189k PF00355 -> RO alpha)"
[ -f ro_alpha.csv ] && echo "  [atlandi] ro_alpha.csv var" || \
  python3 run_ro_filter.py --fasta combined_pfam.fasta --hmm RieskeDB71.hmm --cpu "$CPU"

# --- 2. Genomik baglam (UZUN)
step "2/14  Genomik baglam cikarimi (17k gbk)"
[ -f "$GC/neighbors.csv" ] && echo "  [atlandi] $GC/neighbors.csv var" || \
  python3 extract_genomic_context.py --gbk-dir gbk_files --out-dir "$GC" --processes "$CPU"

# --- 3. Veritabani
step "3/14  SQLite kuruluyor"
python3 build_db.py --context-dir "$GC" --db roar.sqlite

# --- 4. Dogrulama (hmmsearch + hmmalign)
step "4/14  RO adaylari dogrulaniyor (HMM + motif)"
python3 annotate_ro.py --context-dir "$GC" --db roar.sqlite --cpu "$CPU"

# --- 4b. Operonlar (komsu protein dizileri + Pfam HMM; bkz. build_operons.py)
step "4b/14 Operon bilesenleri ve operon turetme"
[ -f component_hmms/components.hmm ] || { echo "EKSIK: component_hmms/components.hmm (bkz. README)"; exit 1; }
python3 build_operons.py --db roar.sqlite --gbk-dir gbk_files --cpu "$CPU"

# --- 4c. Kanit duzeyi, ETC tipi, regulasyon mimarisi
step "4c/14 Kanit duzeyi (kuratorlu referanslara kimlik)"
python3 evidence_tiers.py --db roar.sqlite --out-dir "$OUT" --threads "$CPU"
step "4d/14 Elektron tasima zinciri bilesen tipleri"
python3 etc_types.py --db roar.sqlite --out-dir "$OUT"
step "4e/14 Operon 5' ucu, promotor bolgesi ve regulator"
python3 analyze_regulation.py --db roar.sqlite --out-dir "$OUT"

# --- 5-6. Temel analizler + null model
step "5/14  Genomik baglam analizleri"
python3 analyze.py --db roar.sqlite --out-dir "$OUT"
step "6/14  Null model (arka plan normalizasyonu)"
python3 null_model.py --db roar.sqlite --gbk-dir gbk_files --processes "$CPU" \
  --out "$OUT/null_model.csv"

# --- 7-9. Varyasyon, alt-aile, ozyinelemeli homojenizasyon
step "7/14  Kume ici varyasyon + SDP"
python3 analyze_variants.py --db roar.sqlite --out-dir "$OUT"
step "8/14  Alt-aile kesfi + novel adaylar"
python3 discover_subfamilies.py --db roar.sqlite --out-dir "$OUT"
step "9/14  Ozyinelemeli homojenizasyon (yapraklar)"
python3 recursive_homogenize.py --db roar.sqlite --out-dir "$OUT"

# --- 10-12. Varyant karakterizasyonu, temsilciler, domain
step "10/14  Yaprak/varyant karakterizasyonu"
python3 characterize_leaves.py --db roar.sqlite --out-dir "$OUT"
step "11/14  Temsilci diziler + katalitik dogrulama"
python3 extract_representatives.py --db roar.sqlite --out-dir "$OUT"
step "12/14  Yasam alani siniflandirmasi"
python3 classify_domains.py --db roar.sqlite --out-dir "$OUT"

# --- 13. Ekoloji
step "13/14  Ekoloji hipotez testi"
python3 analyze_ecology.py --db roar.sqlite --out "$OUT/cluster_ecology_stats.csv"

# --- 12b. Varyant kalinti imzasi (ayirt edici kolonlar)
step "12b/14 Varyant kalinti imzasi"
python3 variant_signature.py --db roar.sqlite --out-dir "$OUT"

# --- 12c. Korunmus merkez istatistigi (web ana sayfasi bunu okur)
step "12c/14 Korunmus merkez istatistigi"
python3 motif_stats.py --db roar.sqlite --out "$OUT/motif_stats.json"

# --- 12d. Arama indeksi (FTS5 + filtre alanlari)
step "12d/14 Arama indeksi"
python3 build_search_index.py --db roar.sqlite --ecology cluster_ecology.csv --chemistry chemistry.csv

# --- 13b. Filogeni + dizi benzerlik agi (FastTree gerekir)
step "13b/14 Filogeni ve dizi benzerlik agi"
if command -v FastTree >/dev/null 2>&1 || [ -x ./bin/FastTree ]; then
  FT=$(command -v FastTree || echo ./bin/FastTree)
  python3 build_phylogeny.py --db roar.sqlite --out-dir "$OUT" --threads "$CPU" --fasttree "$FT"
else
  echo "  [atlandi] FastTree yok; agac ve ag uretilmedi (web Atlas filogeni sayfasi bos kalir)"
fi

# --- 13c. Istatistiksel cikarimlar
step "13c/14 Istatistiksel cikarimlar"
python3 stats_overview.py --db roar.sqlite --out "$OUT/stats.json"

# --- 14. HTML sayfalar
step "14/14  HTML sayfalar uretiliyor"
python3 make_report.py   --db roar.sqlite -o report.html
python3 make_explorer.py --db roar.sqlite -o explorer.html
python3 make_hub.py      --db roar.sqlite --out-dir "$OUT" -o hub.html

step "TAMAM"
echo "  Veritabani : roar.sqlite"
echo "  Analizler  : $OUT/"
echo "  Sayfalar   : hub.html (ana kapi), report.html, explorer.html"
echo "  Web        : cd webapp && uvicorn app:app --port 8000   (bkz. webapp/README_DEPLOY.md)"
