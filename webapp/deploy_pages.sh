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
# normal push yapilir; degisiklik yoksa push hic yapilmaz, boylece gereksiz
# derleme tetiklenmez.
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
MIN_FILES=5000
FILE_COUNT=$(find site -type f | wc -l)
for required in site/index.html site/static/style.css site/static/table.js site/search_index.json; do
  if [ ! -s "$required" ]; then
    echo "HATA: derleme eksik, '$required' yok ya da bos."
    echo "      Once tamamlanmis bir derleme uretin:"
    echo "      python3 freeze.py --out site --base /${EXPECTED_BASE}"
    exit 1
  fi
done
if [ "$FILE_COUNT" -lt "$MIN_FILES" ]; then
  echo "HATA: derlemede yalnizca $FILE_COUNT dosya var, en az $MIN_FILES bekleniyor."
  echo "      Derleme yarida kalmis olabilir; tekrar uretip oyle yayinlayin."
  exit 1
fi
echo "[kontrol] derleme tam: $FILE_COUNT dosya"

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
git -c user.name="roar-db" -c user.email="roar-db@localhost" \
    commit -q -m "ROAR-DB static site $(date -u +%F' '%T) UTC ($CHANGED files)"
git push -q origin "$BRANCH"
echo "[gonderildi] $CHANGED dosya -> $URL ($BRANCH)"
echo "  GitHub > Settings > Pages > Branch: $BRANCH"
