#!/usr/bin/env bash
# webapp/site/ icerigini gh-pages dalina gonderir (GitHub Pages).
# Kullanim: bash deploy_pages.sh [repo-url] [branch]
#   varsayilan repo: https://github.com/recepcanaltinbag/pro-sim-blast.git  (Pages: recepcanaltinbag.github.io/pro-sim-blast/)
#   site/ bu repo icin --base /pro-sim-blast ile uretilmis olmali.
set -euo pipefail
cd "$(dirname "$0")"
URL="${1:-https://github.com/recepcanaltinbag/pro-sim-blast.git}"
BRANCH="${2:-gh-pages}"
[ -d site ] || { echo "site/ yok -- once: python3 freeze.py --out site --base /pro-sim-blast"; exit 1; }

TMP=$(mktemp -d)
cp -r site/. "$TMP"/
touch "$TMP/.nojekyll"                 # alt cizgi ile baslayan dosyalari Jekyll yutmasin
cd "$TMP"
git init -q
git checkout -q -b "$BRANCH"
git add -A
git -c user.name="roar-db" -c user.email="roar-db@localhost" commit -q -m "ROAR-DB static site $(date -u +%F)"
git push -f "$URL" "$BRANCH"
cd - >/dev/null
rm -rf "$TMP"
echo "gonderildi: $URL -> $BRANCH. GitHub > Settings > Pages > Branch: $BRANCH"
