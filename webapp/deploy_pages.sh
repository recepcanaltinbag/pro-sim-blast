#!/usr/bin/env bash
# webapp/site/ icerigini gh-pages dalina ARTIMLI olarak gonderir (GitHub Pages).
#
# Kullanim: bash deploy_pages.sh [repo-url] [branch]
#   varsayilan repo: https://github.com/recepcanaltinbag/pro-sim-blast.git
#   (Pages adresi: recepcanaltinbag.github.io/pro-sim-blast/)
#   site/ bu repo icin --base /pro-sim-blast ile uretilmis olmali.
#
# NEDEN ARTIMLI: ilk surum her yayinda gecici bir dizinde `git init` yapip tum
# agaci force-push ediyordu. Site 26.000 dosya / 370 MB'a cikinca bu, her
# seferinde tum icerigi yeniden yuklemek anlamina geldi ve GitHub Pages
# derlemesi geride kaldi (yayindaki kopya push'tan bir saat eski kaldi).
# Simdi kalici bir klon tutulur, yalnizca DEGISEN dosyalar commit edilir ve
# degisiklik yoksa push hic yapilmaz, boylece gereksiz derleme tetiklenmez.
#
# NEDEN TEK COMMIT: gh-pages her yayinda bir commit daha biriktiriyordu ve dal
# sinirsiz buyuyordu -- olculdu, depo GitHub'da 673 MB'a ciktI ve bunun 7 MB'i
# gercekten izlenen icerik; gerisi gecmis. Site main'den bastan uretilebildigi
# icin gh-pages gecmisinin hicbir degeri yok. Bu yuzden her yayin dali TEK bir
# commit olarak yeniden kurar ve force-push eder.
#
# Bu, her dosyayi yeniden YUKLEMEK anlamina GELMEZ: kalici klon onceki commit'i
# bildigi icin git yalnizca uzak tarafta bulunmayan nesneleri gonderir --
# degismeyen dosyalarin blob'lari ve degismeyen dizinlerin agac nesneleri
# zaten orada. Ilk surumun yavas olmasinin sebebi her seferinde gecici bir
# dizinde `git init` yapmasiydi; ortak gecmis olmadigi icin pazarlik
# yapilamiyordu.
set -euo pipefail
cd "$(dirname "$0")"
URL="${1:-https://github.com/recepcanaltinbag/pro-sim-blast.git}"
BRANCH="${2:-gh-pages}"
CLONE=".gh-pages"
[ -d site ] || { echo "site/ yok -- once: python3 freeze.py --out site --base /pro-sim-blast"; exit 1; }

# Taban onegi KONTROL EDILIR. GitHub proje sayfasi siteyi /pro-sim-blast/
# altinda sunuyor, bu yuzden freeze.py --base verilmeden uretilirse butun
# mutlak baglantilar /about.html gibi alan adinin kokune gider ve 404 doner.
# Bu bir kere yayinlandi ve sitenin gezinmesini kirdi; bir daha gecmemesi
# icin deploy burada durur.
EXPECTED_BASE="$(basename "$URL" .git)"
if ! grep -q "href=\"/${EXPECTED_BASE}/static/style.css\"" site/index.html; then
  echo "HATA: site/index.html taban onegi /${EXPECTED_BASE} ile uretilmemis."
  echo "      Butun mutlak baglantilar 404 donerdi. Dogrusu:"
  echo "      python3 freeze.py --out site --base /${EXPECTED_BASE}"
  exit 1
fi
echo "[kontrol] taban onegi /${EXPECTED_BASE} dogrulandi"

# Derlemenin TAM oldugu dogrulanir. Iki freeze.py ayni klasore ayni anda
# yazinca biri otekinin dosyalarini silerken patladi ve geriye static/ klasoru
# BOS olan bir agac kaldi; o agac yayinlaninca sitenin stil dosyasi 404 verdi
# ve tema tamamen gitti. Eksik bir derlemeyi yayinlamak, yayinlamamaktan
# kotudur: site ayakta gorunur ama kirilmistir.
# TAMLIK OLCUTU: SIHIRLI DOSYA SAYISI DEGIL, BEKLENEN SAYILAR.
#
# Eskiden burada MIN_FILES=5000 vardi ve bu, giris basina bir HTML uretilen
# tasarima gore kalibre edilmisti (tam bir site 13.000-26.000 dosya). Giris
# sayfalari tek bir kabuk + data/ parcalariyla degistirilince tam bir derleme
# 3.950 dosyaya dustu ve esik yanlis alarm verdi: 2026-10-24'te hazir ve
# denetimden gecmis bir derlemenin yayini bu yuzden durdu.
#
# Sabit esigin sorunu, neyin EKSIK oldugunu degil "az" oldugunu olcmesiydi.
# Simdi beklenen sayilar kaynaktan okunuyor: kac tip, kac yaprak, kac giris.
# Eksik bir derlemeyi yayinlamak yayinlamamaktan kotudur -- site ayakta
# gorunur ama kirilmistir.
for required in site/index.html site/static/style.css site/static/table.js site/search_index.json \
                site/entry.html site/data/meta.json; do
  if [ ! -s "$required" ]; then
    echo "HATA: derleme eksik, '$required' yok ya da bos."
    echo "      Once tamamlanmis bir derleme uretin:"
    echo "      python3 freeze.py --out site --base /${EXPECTED_BASE}"
    exit 1
  fi
done
FILE_COUNT=$(find site -type f | wc -l)
SHARDS=$(find site/data/entry -name '*.json' 2>/dev/null | wc -l)
SEQ_SHARDS=$(find site/data/seq -name '*.json' 2>/dev/null | wc -l)
CLUSTERS=$(find site/cluster -name '*.html' 2>/dev/null | wc -l)
LEAVES=$(find site/leaf -name '*.html' 2>/dev/null | wc -l)
INDEXED=$(python3 -c "import json;print(len(json.load(open('site/search_index.json'))))" 2>/dev/null || echo 0)

# Beklenenler kaynaktan: chemistry.csv kac tip, veritabani kac yaprak ve kac
# dogrulanmis giris. Boylece olcut veriyle birlikte buyur, elle guncellenmez.
EXP_CLUSTERS=$(python3 -c "import csv;print(sum(1 for _ in csv.DictReader(open('../chemistry.csv'))))" 2>/dev/null || echo 0)
EXP_LEAVES=$(sqlite3 ../roar.sqlite 'SELECT COUNT(*) FROM leaf' 2>/dev/null || echo 0)
EXP_ENTRIES=$(sqlite3 ../roar.sqlite 'SELECT COUNT(*) FROM ro WHERE is_confirmed=1' 2>/dev/null || echo 0)

# Beklenen PARCA KUMESI arama indeksinden hesaplanir: her kimligin parca
# anahtari md5'inin ilk hanelerinden geliyor, yani hangi parcalarin var olmasi
# gerektigi tam olarak bilinebilir. Ilk surumde burada MIN_SHARDS=200 esigi
# vardi ve 16 parca silindiginde 240 kalip ESIGI GECTI -- yani koruma bos
# cikti. Sihirli sayi yerine kesin kume.
EXP_SHARDS=$(python3 - <<'PYEOF' 2>/dev/null || echo 0
import hashlib, json, re
m = re.compile(r"^(.*)_(\d+)-(\d+)_(-?[01])$")
def unslug(s):
    g = m.match(s)
    return f"{g.group(1)}:{g.group(2)}-{g.group(3)}:{g.group(4)}" if g else s.replace("_", ":")
hexlen = json.load(open("site/data/meta.json")).get("shard_hex", 2)
rows = json.load(open("site/search_index.json"))
print(len({hashlib.md5(unslug(r[0]).encode()).hexdigest()[:hexlen] for r in rows}))
PYEOF
)

fail=0
if [ "$EXP_SHARDS" -gt 0 ] && { [ "$SHARDS" -ne "$EXP_SHARDS" ] || [ "$SEQ_SHARDS" -ne "$EXP_SHARDS" ]; }; then
  echo "HATA: veri parcasi eksik -- entry $SHARDS, seq $SEQ_SHARDS, beklenen $EXP_SHARDS."
  echo "      Giris sayfalari bu parcalardan ciziliyor; eksikse o girisler bos acilir."
  fail=1
fi
if [ "$EXP_CLUSTERS" -gt 0 ] && [ "$CLUSTERS" -lt "$EXP_CLUSTERS" ]; then
  echo "HATA: $CLUSTERS tip sayfasi var, chemistry.csv $EXP_CLUSTERS tip sayiyor."
  fail=1
fi
if [ "$EXP_LEAVES" -gt 0 ] && [ "$LEAVES" -lt "$EXP_LEAVES" ]; then
  echo "HATA: $LEAVES varyant sayfasi var, veritabani $EXP_LEAVES yaprak sayiyor."
  fail=1
fi
if [ "$EXP_ENTRIES" -gt 0 ] && [ "$INDEXED" -lt "$EXP_ENTRIES" ]; then
  echo "HATA: arama indeksinde $INDEXED giris var, veritabani $EXP_ENTRIES dogrulanmis giris sayiyor."
  fail=1
fi
if [ "$fail" -ne 0 ]; then
  echo "      Tekrar uretip oyle yayinlayin:"
  echo "      python3 freeze.py --out site --base /${EXPECTED_BASE}"
  exit 1
fi
echo "[kontrol] derleme tam: $FILE_COUNT dosya, $SHARDS/$EXP_SHARDS parca, "\
     "$CLUSTERS/$EXP_CLUSTERS tip, $LEAVES/$EXP_LEAVES varyant, $INDEXED/$EXP_ENTRIES giris"

if [ ! -d "$CLONE/.git" ]; then
  echo "[kuruluyor] kalici klon: $CLONE"
  rm -rf "$CLONE"
  if git clone --depth 1 --branch "$BRANCH" "$URL" "$CLONE" 2>/dev/null; then
    echo "[bilgi] mevcut $BRANCH dali alindi"
  else
    echo "[bilgi] $BRANCH dali yok, bos bir dal kuruluyor"
    git clone --depth 1 "$URL" "$CLONE"
    git -C "$CLONE" checkout -q --orphan "$BRANCH"
    git -C "$CLONE" rm -rq --cached . 2>/dev/null || true
    find "$CLONE" -mindepth 1 -maxdepth 1 ! -name .git -exec rm -rf {} +
  fi
fi

if git -C "$CLONE" fetch -q --depth 1 origin "$BRANCH" 2>/dev/null; then
  git -C "$CLONE" reset -q --hard FETCH_HEAD
fi

# site/ -> klon: silinen dosyalar da yansitilir, .git korunur
rsync -a --delete --exclude '.git' site/ "$CLONE"/
touch "$CLONE/.nojekyll"          # alt cizgi ile baslayan dosyalari Jekyll yutmasin

cd "$CLONE"
git add -A
if git diff --cached --quiet; then
  echo "[atlandi] degisiklik yok, push yapilmadi"
  exit 0
fi
CHANGED=$(git diff --cached --name-only | wc -l)

# Dali tek commit olarak yeniden kur. Calisma agaci oldugu gibi kalir; yalnizca
# gecmis atilir.
git checkout -q --orphan __deploy
git add -A
git -c user.name="roar-db" -c user.email="roar-db@localhost" \
    commit -q -m "ROAR-DB static site $(date -u +%F' '%T) UTC ($CHANGED files changed)"
git branch -q -M __deploy "$BRANCH"
git push -q --force origin "$BRANCH"
echo "[gonderildi] $CHANGED degisen dosya -> $URL ($BRANCH, tek commit)"
echo "  GitHub > Settings > Pages > Branch: $BRANCH"
