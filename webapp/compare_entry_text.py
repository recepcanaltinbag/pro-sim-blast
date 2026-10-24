"""Tarayicida cizilen giris sayfasini Jinja'nin urettigiyle karsilastir.

NEDEN. Giris sayfasi artik iki yolla uretilebiliyor: eski sunucu tarafi
sablon (templates/ro.html) ve yeni istemci tarafi render (static/entry.js).
Ikisi ayrisirsa statik sitede yanlis ya da eksik bilgi yayimlanir ve bunu
kimse fark etmez. Burada GORUNEN METIN karsilastirilir.

Bicim farklari (bosluk, satir sonu) normalize edilir; karsilastirilan sey
anlam tasiyan metindir.

Kullanim:
    python3 compare_entry_text.py --ids 20 --serve-url http://127.0.0.1:8732/pro-sim-blast
"""
import argparse, os, re, subprocess, sqlite3, sys, tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)


def visible_text(html):
    from bs4 import BeautifulSoup
    soup = BeautifulSoup(html, "html.parser")
    for tag in soup(["script", "style", "title", "svg"]):
        tag.decompose()
    main = soup.find("main") or soup
    txt = main.get_text(" ", strip=True)
    return re.sub(r"\s+", " ", txt)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--db", default=os.path.join(os.path.dirname(HERE), "roar.sqlite"))
    ap.add_argument("--serve-url", required=True)
    ap.add_argument("--ids", type=int, default=20)
    ap.add_argument("--show", type=int, default=2)
    a = ap.parse_args()

    os.environ.setdefault("ROAR_BASE", "/pro-sim-blast")
    import app as A                                   # noqa: PLC0415
    # STATIC_MODE ortam degiskeniyle degil, freeze.py tarafindan ELLE set
    # ediliyor. Karsilastirmanin anlamli olmasi icin iki taraf da ayni kipte
    # olmali: aksi halde Sequence bolumu Jinja'da indirme dugmeleri, JS'te
    # statik metin gosterir ve fark gercek bir ayrisma sanilir.
    A.STATIC_MODE = True
    from fastapi.testclient import TestClient         # noqa: PLC0415
    client = TestClient(A.app)

    con = sqlite3.connect(f"file:{a.db}?mode=ro", uri=True)
    ids = [r[0] for r in con.execute(
        "SELECT candidate_id FROM ro WHERE is_confirmed=1 "
        f"ORDER BY candidate_id LIMIT {a.ids}")]
    con.close()

    bad = []
    for n, cid in enumerate(ids, 1):
        slug = A.slug(cid)
        jinja = visible_text(client.get(f"/ro/{slug}").text)
        with tempfile.NamedTemporaryFile("w+", suffix=".html", delete=False) as fh:
            out = fh.name
        r = subprocess.run(
            ["chromium", "--headless", "--disable-gpu", "--no-sandbox",
             "--virtual-time-budget=10000", "--dump-dom",
             f"{a.serve_url}/entry.html?id={slug}"],
            capture_output=True, text=True, timeout=180)
        js = visible_text(r.stdout)
        os.unlink(out)
        if jinja != js:
            bad.append((cid, jinja, js))
        print(f"   {n}/{len(ids)} {cid[:46]:48s} {'FARKLI' if jinja != js else 'ayni'}")

    print()
    if not bad:
        print(f"AYNI: {len(ids)} giriste gorunen metin birebir ortusuyor.")
        return 0
    print(f"FARK: {len(bad)}/{len(ids)}")
    for cid, j, s in bad[:a.show]:
        i = next((k for k in range(min(len(j), len(s))) if j[k] != s[k]), min(len(j), len(s)))
        print(f"\n  {cid}  (ilk fark {i}. karakterde; jinja {len(j)}, js {len(s)})")
        print(f"    jinja: ...{j[max(0,i-70):i+70]}...")
        print(f"    js   : ...{s[max(0,i-70):i+70]}...")
    return 1


if __name__ == "__main__":
    sys.exit(main())
