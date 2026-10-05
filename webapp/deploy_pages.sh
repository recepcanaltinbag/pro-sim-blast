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
