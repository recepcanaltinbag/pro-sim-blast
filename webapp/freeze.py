"""
Static export of ROAR-DB for free static hosting (GitHub Pages, Netlify, Cloudflare Pages).

Renders every page through the FastAPI app with STATIC_MODE on (links get .html),
writes a client-side search index, and copies downloads. The sequence classifier
needs a server and is omitted from the static site (the Docker/HF Space has it).

    python3 freeze.py --out site --base /pro-sim-blast      # base = repo name for GitHub project pages
    python3 freeze.py --out site --base ""             # user/organisation pages or custom domain

Per-entry and per-variant download files are deliberately NOT exported: they would
add about 13,000 files and 50 MB while the sequence is already shown on the page
and the bulk FASTA covers every entry. GitHub Pages builds a site of this size
slowly, so file count is worth spending carefully.
"""

import argparse
import json
import os
import shutil
import csv
import sqlite3
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", default=os.path.join(HERE, "site"))
    ap.add_argument("--base", default="", help="URL prefix, e.g. /pro-sim-blast for a GitHub project page")
    ap.add_argument("--limit", type=int, default=0, help="only N entry pages (for testing)")
    ap.add_argument("--pages-only", action="store_true",
                    help="refresh the overview, type and variant pages plus static files and "
                         "bulk downloads; leave the per-entry pages in place (they are the slow "
                         "part and change only when the database itself changes)")
    args = ap.parse_args()

    os.environ["ROAR_BASE"] = args.base
    import app as A
    A.STATIC_MODE = True
    A.PAGE_SIZE = 10 ** 6      # static pages carry the full member list (no ?page= files)
    A.BASE = args.base.rstrip("/")
    from fastapi.testclient import TestClient
    client = TestClient(A.app)

    out = args.out
    if os.path.exists(out) and not args.pages_only:
        shutil.rmtree(out)
    os.makedirs(out, exist_ok=True)
    static_dst = os.path.join(out, "static")
    if os.path.exists(static_dst):
        shutil.rmtree(static_dst)
    shutil.copytree(os.path.join(HERE, "static"), static_dst)

    def save(url, path, binary=False):
        r = client.get(url)
        if r.status_code != 200:
            print(f"  !! {r.status_code} {url}")
            return
        full = os.path.join(out, path.lstrip("/"))
        os.makedirs(os.path.dirname(full), exist_ok=True)
        with open(full, "wb") as fh:
            fh.write(r.content)

    con = sqlite3.connect(f"file:{A.DB_PATH}?mode=ro", uri=True)
    # Tip listesi, UYESI OLANLARLA sinirli kalamaz. Kuratorlu 71 tipin 10'u bu
    # derlemede tek bir dogrulanmis uye toplamiyor (ornegin NahAc ve CarAa,
    # birebir ikizleri varken profil atamayi digerine kaydiriyor). Yine de
    # giris sayfalari "en yakin referans" olarak onlara BAGLANIYOR, ve bu
    # sayfalar uretilmezse statik sitede 90 baglanti 404 donuyordu. Dinamik
    # uygulama bu tipler icin bos bir sayfa veriyor; ihracat da onu almali.
    clusters = {r[0] for r in con.execute(
        "SELECT DISTINCT ro_cluster FROM ro WHERE is_confirmed=1") if r[0]}
    chem_path = os.path.join(A.PARENT, "chemistry.csv")
    if os.path.exists(chem_path):
        with open(chem_path, newline="", encoding="utf-8") as fh:
            clusters |= {row["cluster"] for row in csv.DictReader(fh) if row.get("cluster")}
    clusters = sorted(clusters)
    leaves = [r[0] for r in con.execute("SELECT leaf_id FROM leaf")]
    entries = [r[0] for r in con.execute(
        "SELECT candidate_id FROM ro WHERE is_confirmed=1 ORDER BY candidate_id")]
    if args.limit:
        entries = entries[:args.limit]

    print("[pages] fixed")
    save("/", "index.html"); save("/clusters", "clusters.html"); save("/about", "about.html")
    print("[pages] atlas")
    for p in ("atlas", "atlas/phylogeny", "atlas/network", "atlas/taxonomy", "atlas/operons",
              "atlas/regulation", "atlas/ecology", "atlas/evidence", "atlas/novel",
              "atlas/cooccurrence", "atlas/statistics", "atlas/quality",
              "atlas/closeness", "atlas/structure", "atlas/control", "atlas/origin", "atlas/learning",
              "atlas/geography", "atlas/elements", "atlas/disagreements"):
        save("/" + p, p + ".html")
    save("/download/tree_all.nwk", "download/tree_all.nwk")
    save("/download/novel_candidates.fasta", "download/novel_candidates.fasta")
    # Liste atlas.PROVENANCE'tan turetilir. Daha once burada ve indirme
    # rotasinda iki ayri kopya vardi; biri guncellenip oteki unutulunca dosya
    # ya sitede gorunup indirilemiyor ya da indirilip belgelenmemis oluyordu.
    import atlas as _atlas
    for name in sorted(_atlas.ANALYSIS_FILES):
        if name.endswith(".nwk"):
            continue          # agac /download/tree_all.nwk yolundan gidiyor
        save("/download/analysis/" + name, "download/analysis/" + name)
    print(f"[pages] {len(clusters)} clusters")
    for c in clusters:
        save(f"/cluster/{c}", f"cluster/{c}.html")
        save(f"/download/cluster/{c}.fasta", f"download/cluster/{c}.fasta")
        save(f"/download/cluster/{c}.csv", f"download/cluster/{c}.csv")
    # Duzenleyici ve IS ailelerinin sayfalari. Adlar JSON'dan geliyor, yani
    # yeni bir aile esigi gectiginde sayfasi kendiliginden ihrac ediliyor.
    families = _atlas.control_element_names(
        os.path.join(_atlas._DEFAULT_ANALYSIS_DIR, "control_elements.json"))
    if families:
        print(f"[pages] {len(families)} regulator and IS family pages")
        for fam in families:
            slug = _atlas.element_slug(fam)
            save(f"/element/{slug}", f"element/{slug}.html")

    print(f"[pages] {len(leaves)} leaves")
    for l in leaves:
        save(f"/leaf/{A.leaf_slug(l)}", f"leaf/{A.leaf_slug(l)}.html")
    if args.pages_only:
        print("[search] index")
        write_search_index(con, out, A)
        con.close()
        print(f"[done] {out} (entry pages left untouched)")
        return

    print(f"[pages] {len(entries)} entries")
    for i, e in enumerate(entries, 1):
        s = A.slug(e)
        save(f"/ro/{s}", f"ro/{s}.html")
        if i % 1000 == 0:
            print(f"   {i}/{len(entries)}")
    print("[downloads] bulk")
    save("/download/all_confirmed.fasta", "download/all_confirmed.fasta")
    save("/download/all_confirmed.csv", "download/all_confirmed.csv")
    save("/download/operons.csv", "download/operons.csv")

    # client-side search index
    print("[search] index")
    write_search_index(con, out, A)
    con.close()
    print(f"[done] {out}")



def write_fallback_pages(out, A, con):
    """Adresi elle yazan ya da eski bir baglantiyi takip eden okuyucu icin.

    NEDEN. Yayindaki site olculdu: `/classify` ve `/classify.html` 404
    donuyordu, cunku dizi siniflandirici sunucu gerektiriyor ve statik ihracata
    hic girmiyor. Gezinme cubugu baglantiyi gizliyor ama ADRES tahmin edilebilir
    ve eski derlemelerde vardi, yani yer imi olan biri bos bir GitHub 404
    sayfasina dusuyor. Ayni sekilde `/atlas/` (sondaki egik cizgiyle) 404
    donuyordu, `/atlas` ise 200: Pages, egik cizgili yolda `index.html` ariyor.
    Ve sitenin HIC 404 sayfasi yoktu, dolayisiyla yanlis bir adres okuyucuyu
    siteden tamamen cikariyordu.

    Uc sayfa yazilir: bir 404, siteye donus yollariyla; bir `classify` sayfasi,
    neden burada calismadigini ve nerede calistigini soyleyen; ve her atlas
    sayfasi icin egik cizgili surumun karsiligi.
    """
    def page(title, heading, body):
        return (f'<!doctype html><html lang="en"><head><meta charset="utf-8">'
                f'<meta name="viewport" content="width=device-width, initial-scale=1">'
                f'<title>{title} · ROAR-DB</title>'
                f'<link rel="stylesheet" href="{A.BASE}/static/style.css"></head><body>'
                f'<header class="top"><a class="brand" href="{A.BASE}/index.html">'
                f'<span class="brand__mark">RO</span> ROAR-DB</a>'
                f'<nav><a href="{A.BASE}/atlas.html">Atlas</a>'
                f'<a href="{A.BASE}/clusters.html">Types</a>'
                f'<a href="{A.BASE}/search.html">Search</a>'
                f'<a href="{A.BASE}/about.html">Methods</a></nav></header>'
                f'<main><h1>{heading}</h1>{body}</main></body></html>')

    links = (f'<ul><li><a href="{A.BASE}/index.html">Home: the shared architecture, the '
             f'reactions and the enzyme types grouped by chemistry</a></li>'
             f'<li><a href="{A.BASE}/search.html">Search by organism, protein identifier, '
             f'enzyme type or a chemical</a></li>'
             f'<li><a href="{A.BASE}/clusters.html">Every enzyme type</a></li>'
             f'<li><a href="{A.BASE}/atlas.html">Atlas: phylogeny, sequence space, ecology, '
             f'operons, regulation, evidence and the statistics</a></li>'
             f'<li><a href="{A.BASE}/about.html">Methods, with every threshold and its '
             f'justification</a></li></ul>')

    with open(os.path.join(out, "404.html"), "w") as fh:
        fh.write(page("Page not found", "That page is not here",
                      '<p class="lede">The address does not match any page of this database. '
                      'Two common reasons: a link from an older version of the site, or a '
                      'trailing slash. Everything below is a working entry point.</p>'
                      + links +
                      '<p class="note">If you reached this from a link inside the site, that is '
                      'a fault worth reporting, because every internal link is checked '
                      'automatically before each release.</p>'))

    # Siniflandirici: 404 yerine NEDEN burada olmadigini soyleyen bir sayfa.
    with open(os.path.join(out, "classify.html"), "w") as fh:
        fh.write(page("Classify a sequence", "Classifying your own sequence",
                      '<p class="lede">This page needs a server and the published site does not '
                      'have one, so the classifier is not available here.</p>'
                      '<p>Classification runs the same two tests every entry in this database '
                      'passed: profile coverage against the 71 reference models, and a direct '
                      'check of the eight catalytic-centre columns. Both require HMMER, which '
                      'cannot run in a browser.</p>'
                      '<h2>Where it does run</h2>'
                      '<p>The full application, including the classifier and the faceted '
                      'full-text search, is in the repository with a Dockerfile. Running it '
                      'locally gives the complete feature set over the same database.</p>'
                      '<pre>git clone https://github.com/recepcanaltinbag/pro-sim-blast\n'
                      'cd pro-sim-blast/webapp\n'
                      'docker build -t roar-db . &amp;&amp; docker run -p 8000:8000 roar-db</pre>'
                      '<h2>What you can do here instead</h2>'
                      f'<p>If you have a protein identifier or an organism, the '
                      f'<a href="{A.BASE}/search.html">search</a> will find it among the '
                      f'11,422 confirmed entries. If you want to know what defines membership, '
                      f'the <a href="{A.BASE}/about.html">methods page</a> gives every '
                      f'threshold and the measurement behind it.</p>'))

    # Egik cizgili atlas adresleri: Pages bu yolda index.html ariyor.
    for name in ("atlas", "cluster", "leaf", "ro", "download"):
        folder = os.path.join(out, name)
        if not os.path.isdir(folder):
            continue
        target = f"{A.BASE}/{name}.html" if name == "atlas" else f"{A.BASE}/index.html"
        with open(os.path.join(folder, "index.html"), "w") as fh:
            fh.write('<!doctype html><html lang="en"><head><meta charset="utf-8">'
                     f'<meta http-equiv="refresh" content="0; url={target}">'
                     f'<link rel="canonical" href="{target}">'
                     '<title>Redirecting · ROAR-DB</title></head><body>'
                     f'<p>Redirecting to <a href="{target}">{target}</a>.</p>'
                     '</body></html>')
    # Yeniden adlandirilmis tiplerin ESKI adresleri. Statik sitede sunucu
    # yonlendirmesi yok, bu yuzden meta-refresh tasiyan kucuk sayfalar konur.
    # Tip sayfasi TEK basina yetmiyor: ayni ad varyant sayfalarinin ve indirme
    # dosyalarinin yolunda da geciyordu, ve onlar da disaridan baglanmis
    # olabilir. Hepsi icin yonlendirme yazilir.
    renamed = getattr(A, "RENAMED_TYPES", {}) or {}
    leaf_ids = [r[0] for r in con.execute("SELECT leaf_id FROM leaf")]
    redirects = 0

    def write_redirect(path, target, label):
        nonlocal redirects
        full = os.path.join(out, path)
        os.makedirs(os.path.dirname(full), exist_ok=True)
        with open(full, "w") as fh:
            fh.write('<!doctype html><html lang="en"><head><meta charset="utf-8">'
                     f'<meta http-equiv="refresh" content="0; url={target}">'
                     f'<link rel="canonical" href="{target}">'
                     f'<title>Renamed · ROAR-DB</title></head><body>'
                     f'<p>{label} was renamed. Continuing to '
                     f'<a href="{target}">{target}</a>.</p></body></html>')
        redirects += 1

    for old_name, new_name in renamed.items():
        write_redirect(os.path.join("cluster", old_name + ".html"),
                       f"{A.BASE}/cluster/{new_name}.html", "This type")
        # Varyant sayfalari: yaprak kimligi tip adini tasiyor.
        # OLCULDU: bu dongu `leaves` adini kullaniyordu ama o degisken baska bir
        # fonksiyonun yerelindeydi; derleme son adimda NameError ile duruyor ve
        # 404 sayfasi, classify uyarisi ve BUTUN yeniden adlandirma
        # yonlendirmeleri hic yazilmiyordu. Liste artik burada sorgulanir.
        for leaf_id in leaf_ids:
            if leaf_id.startswith(new_name + "#"):
                suffix = leaf_id.split("#", 1)[1]
                write_redirect(os.path.join("leaf", f"{old_name}-{suffix}.html"),
                               f"{A.BASE}/leaf/{new_name}-{suffix}.html", "This variant")
        # Indirme dosyasi: eski yolu da calissin.
        write_redirect(os.path.join("download", "cluster", old_name + ".fasta.html"),
                       f"{A.BASE}/download/cluster/{new_name}.fasta", "This download")
    if redirects:
        print(f"[pages] {redirects} redirect(s) for renamed types")
    print("[pages] 404, classify notice and directory redirects")

def write_search_index(con, out, A):
    """Client-side index plus the static search page.

    The static site has no server, so the facets of the hosted application are
    reproduced in the browser.

    SERBEST METIN de eslesmek zorunda. Bir donem burada yalnizca giris basina
    alanlar yaziliyordu (organizma, urun, tip, varyant, aile) ve KIMYA yoktu;
    sonuc olarak yayinlanan arama bir kimyasali bulamiyordu. Olculdu:
    "benzalkonium" uygulamada 95 giris donuyordu, sitede 0; terephthalate
    273'e 0; caffeine 447'ye 0. Sitenin butun duzeni "enzimleri etkiledikleri
    kimyasallara gore grupla" oldugu icin bu kucuk bir eksik degildi.

    Substrat giris basina DEGIL tip basina bir ozellik, bu yuzden 11.422 satira
    tekrar tekrar yazilmaz: tip basina kucuk bir sozluk (71 kayit) sayfaya
    gomulur ve tarayici aranan metni oradan tamamlar.
    """
    has_search = con.execute(
        "SELECT name FROM sqlite_master WHERE name='ro_search'").fetchone()
    if has_search:
        rows = con.execute("""
            SELECT candidate_id, protein_id, organism, product, cluster, leaf_id,
                   tier, domain, family, is_plasmid, has_partner, has_regulator
            FROM ro_search""").fetchall()
        data = [[A.slug(r[0]), r[1] or "", r[2] or "", r[3] or "", r[4] or "", r[5] or "",
                 r[6] or "", r[7] or "", r[8] or "", r[9] or 0, r[10] or 0, r[11] or 0]
                for r in rows]
    else:
        rows = con.execute("""
            SELECT r.candidate_id, r.protein_id, p.organism, r.product, r.ro_cluster,
                   '', '', '', '', p.is_plasmid, 0, 0
            FROM ro r JOIN replicon p USING(nucleotide_id) WHERE r.is_confirmed=1""").fetchall()
        data = [[A.slug(r[0])] + list(r[1:]) for r in rows]
    # Tip basina kimya: gosterilecek kisa ad, substrat, reaksiyon sinifi.
    # Kisa ad `A.gene` ile alinir, cunku uc kisa ad cakisiyor ve statik arama
    # sayfasinin da tablolarla AYNI etiketi gostermesi gerekiyor.
    chem = {}
    chem_path = os.path.join(A.PARENT, "chemistry.csv")
    if os.path.exists(chem_path):
        with open(chem_path, newline="", encoding="utf-8") as fh:
            for row in csv.DictReader(fh):
                cluster = row.get("cluster")
                if not cluster:
                    continue
                chem[cluster] = [A.gene(cluster),
                                 row.get("substrate_en") or "",
                                 (row.get("reaction_class") or "").replace("_", " ")]
    for cluster in {r[4] for r in data if r[4]}:
        chem.setdefault(cluster, [A.gene(cluster), "", ""])

    with open(os.path.join(out, "search_index.json"), "w") as fh:
        json.dump(data, fh, separators=(",", ":"))
    with open(os.path.join(out, "search.html"), "w") as fh:
        fh.write(SEARCH_PAGE.replace("__BASE__", A.BASE)
                            .replace("__CHEM__", json.dumps(chem, separators=(",", ":"))))
    print(f"   {len(data)} entries indexed, {len(chem)} types with chemistry")
    write_fallback_pages(out, A, con)


SEARCH_PAGE = """<!doctype html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1"><title>Search · ROAR-DB</title>
<link rel="stylesheet" href="__BASE__/static/style.css"></head><body>
<header class="top"><a class="brand" href="__BASE__/index.html"><span class="brand__mark">RO</span> ROAR-DB</a>
<nav><a href="__BASE__/atlas.html">Atlas</a><a href="__BASE__/clusters.html">Types</a>
<a href="__BASE__/search.html">Search</a><a href="__BASE__/about.html">Methods</a></nav></header>
<main>
<h1>Search</h1>
<form class="bigsearch" onsubmit="return false">
  <input id="q" type="search" placeholder="organism, product, protein identifier, enzyme type, or a chemical such as benzalkonium" autofocus>
  <button onclick="run()">Search</button>
</form>
<p class="note">Several words are combined with AND. A trailing asterisk matches a prefix.
Filters on the left narrow the result without a new search.
<span id="status">loading index…</span></p>
<div class="searchlayout">
  <aside class="facets" id="facets"></aside>
  <div class="results">
    <div id="typeblock"></div>
    <h2>Entries <span class="note" id="count"></span></h2>
    <div class="tablewrap"><table class="data">
      <thead><tr><th>Entry</th><th>Organism</th><th>Product</th><th>Type</th>
      <th>Substrate of the reference</th><th>Evidence</th><th>Variant</th></tr></thead>
      <tbody id="rows"></tbody></table></div>
    <nav class="pager" id="pager"></nav>
  </div>
</div>
</main>
<script>
const BASE = "__BASE__";
const TIER_COLORS = {characterized:'#1f7a4d', close_homolog:'#7fb069', family_member:'#f2c14e',
                     distant:'#f78154', novel:'#8e44ad'};
/*  Tip basina kimya. Arama metni buradan tamamlanir, yoksa bir kimyasal adi
    hicbir giriste yazili olmadigi icin bulunamaz.                           */
const CHEM = __CHEM__;
const gene = c => (CHEM[c] && CHEM[c][0]) || c.split('_').slice(2).join('_');
const F = {q:'', tier:'', domain:'', family:'', cluster:'', plasmid:0, partner:0, regulator:0, page:1};
const PAGE = 100;
let IDX = [], VIEW = [];
/*  Etiket metni arka plandan HESAPLANIR, sabit beyaz degil. Olculdu: alti
    etikette beyaz metin WCAG oranini gecmiyordu, en kotusu 1,68:1. Sunucu
    tarafinda `atlas.chip_text` ayni kurali uyguluyor.                      */
function chipText(bg) {
  var v = String(bg).replace('#', '');
  if (v.length !== 6) return '#ffffff';
  var ch = [0, 2, 4].map(function (i) {
    var c = parseInt(v.substr(i, 2), 16) / 255;
    return c <= 0.03928 ? c / 12.92 : Math.pow((c + 0.055) / 1.055, 2.4);
  });
  var L = 0.2126 * ch[0] + 0.7152 * ch[1] + 0.0722 * ch[2];
  var withDark = (L + 0.05) / (0.0114 + 0.05);
  var withLight = (1.05) / (L + 0.05);
  return withDark > withLight ? '#10161c' : '#ffffff';
}
const esc = s => String(s==null?'':s).replace(/[&<>"']/g, c => ({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
// index columns: 0 slug 1 protein 2 organism 3 product 4 cluster 5 leaf 6 tier 7 domain 8 family 9 plasmid 10 partner 11 regulator
const TEXT = r => {
  const c = CHEM[r[4]] || ['', '', ''];
  return (r[0]+' '+r[1]+' '+r[2]+' '+r[3]+' '+r[4]+' '+r[5]+' '+r[8]+' '+
          c[0]+' '+c[1]+' '+c[2]).toLowerCase();
};

fetch(BASE + '/search_index.json').then(r => r.json()).then(d => {
  IDX = d.map(r => { r.push(TEXT(r)); return r; });
  document.getElementById('status').textContent = d.length.toLocaleString() + ' entries indexed.';
  const p = new URLSearchParams(location.search);
  ['q','tier','domain','family','cluster'].forEach(k => { if (p.get(k)) F[k] = p.get(k); });
  ['plasmid','partner','regulator'].forEach(k => { if (p.get(k)) F[k] = 1; });
  document.getElementById('q').value = F.q;
  run();
});

function matches(r) {
  if (F.tier && r[6] !== F.tier) return false;
  if (F.domain && r[7] !== F.domain) return false;
  if (F.family && r[8] !== F.family) return false;
  if (F.cluster && r[4] !== F.cluster) return false;
  if (F.plasmid && !r[9]) return false;
  if (F.partner && !r[10]) return false;
  if (F.regulator && !r[11]) return false;
  if (!F.q) return true;
  const hay = r[12];
  return F.q.toLowerCase().split(/\s+/).filter(Boolean).every(tok =>
    tok.endsWith('*') ? hay.includes(tok.slice(0, -1)) : hay.includes(tok));
}

function run() {
  F.q = document.getElementById('q').value.trim();
  F.page = 1;
  render();
}

function setFilter(k, v) { F[k] = (F[k] === v) ? '' : v; F.page = 1; render(); }
function toggle(k) { F[k] = F[k] ? 0 : 1; F.page = 1; render(); }
function goPage(n) { F.page = n; render(); window.scrollTo(0, 0); }

/*  Sorguya uyan ENZIM TIPLERI ayrica gosterilir. Uye sayisi sifir olan 10 tip
    var ve bunlarin substratlari giris indeksinde hic gecmiyor; yalnizca
    girisleri arayan bir sayfada "xylene" aramasi bos donuyordu, oysa o
    kimyasalin kuratorlu bir tipi VAR. Barindirilan surumde bu bolum zaten
    vardi; statik surumde eksikti.                                          */
function renderTypes() {
  const block = document.getElementById('typeblock');
  const q = F.q.toLowerCase().split(/\s+/).filter(Boolean);
  if (!q.length) { block.innerHTML = ''; return; }
  const counts = {};
  IDX.forEach(r => { counts[r[4]] = (counts[r[4]] || 0) + 1; });
  const hits = Object.keys(CHEM).filter(c => {
    const v = CHEM[c];
    const hay = (c + ' ' + v[0] + ' ' + v[1] + ' ' + v[2]).toLowerCase();
    return q.every(tok => tok.endsWith('*') ? hay.includes(tok.slice(0, -1)) : hay.includes(tok));
  }).sort((a, b) => (counts[b] || 0) - (counts[a] || 0));
  if (!hits.length) { block.innerHTML = ''; return; }
  block.innerHTML = '<h2>Enzyme types matching the query <span class="note">' + hits.length +
    '</span></h2><div class="tablewrap"><table class="data compact"><thead><tr>' +
    '<th>Type</th><th>Substrate of the reference</th><th>Reaction</th>' +
    '<th class="num">Members</th></tr></thead><tbody>' +
    hits.map(c => {
      const v = CHEM[c], n = counts[c] || 0;
      const note = n ? '' : ' <span class="note">no confirmed member in this build</span>';
      return `<tr><td><a href="${BASE}/cluster/${esc(c)}.html">${esc(v[0])}</a>${note}</td>` +
             `<td>${esc(v[1])}</td><td class="small">${esc(v[2])}</td>` +
             `<td class="num" data-v="${n}">${n.toLocaleString()}</td></tr>`;
    }).join('') + '</tbody></table></div>';
}

function render() {
  renderTypes();
  VIEW = IDX.filter(matches);
  document.getElementById('count').textContent = VIEW.length.toLocaleString();
  const start = (F.page - 1) * PAGE;
  document.getElementById('rows').innerHTML = VIEW.slice(start, start + PAGE).map(r => {
    const tier = r[6] ? `<span class="chip" style="background:${TIER_COLORS[r[6]]||'#888'};color:${chipText(TIER_COLORS[r[6]]||'#888')};border:0">${esc(r[6].replace('_',' '))}</span>` : '';
    const leaf = r[5] ? `<a href="${BASE}/leaf/${esc(r[5].replace('#','-'))}.html">${esc(r[5].split('#').pop())}</a>` : '–';
    const dom = (r[7] && r[7] !== 'Bacteria') ? ` <span class="chip">${esc(r[7])}</span>` : '';
    const pl = r[9] ? ' <span class="chip chip--plasmid">plasmid</span>' : '';
    const sub = (CHEM[r[4]] && CHEM[r[4]][1]) || '';
    return `<tr><td><a href="${BASE}/ro/${esc(r[0])}.html">${esc(r[1] || r[0])}</a>${pl}</td>` +
           `<td><i>${esc(r[2])}</i>${dom}</td><td class="small">${esc(r[3])}</td>` +
           `<td><a href="${BASE}/cluster/${esc(r[4])}.html">${esc(gene(r[4]))}</a></td>` +
           `<td class="small">${esc(sub)}</td>` +
           `<td>${tier}</td><td>${leaf}</td></tr>`;
  }).join('') || '<tr><td colspan="7" class="empty">No entry matches. Try fewer words or remove a filter.</td></tr>';

  const pages = Math.ceil(VIEW.length / PAGE);
  document.getElementById('pager').innerHTML = pages > 1
    ? (F.page > 1 ? `<a href="#" onclick="goPage(${F.page-1});return false">previous</a>` : '') +
      `<span>page ${F.page} of ${pages}</span>` +
      (F.page < pages ? `<a href="#" onclick="goPage(${F.page+1});return false">next</a>` : '')
    : '';
  renderFacets();
}

function renderFacets() {
  const defs = [['tier', 'Evidence level', 6], ['domain', 'Domain of life', 7],
                ['family', 'Chemical family', 8], ['cluster', 'Enzyme type', 4]];
  let html = '';
  for (const [key, label, col] of defs) {
    const counts = new Map();
    for (const r of IDX) {
      const saved = F[key]; F[key] = '';
      if (matches(r) && r[col]) counts.set(r[col], (counts.get(r[col]) || 0) + 1);
      F[key] = saved;
    }
    const items = [...counts.entries()].sort((a, b) => b[1] - a[1]).slice(0, 12);
    if (!items.length) continue;
    html += `<div class="facet"><h3>${label}</h3><ul>` + items.map(([v, n]) =>
      `<li><a href="#" onclick="setFilter('${key}', ${JSON.stringify(v)});return false">` +
      (F[key] === v ? `<b>${esc(v)}</b>` : esc(v)) + `</a><span class="note">${n.toLocaleString()}</span></li>`).join('') + '</ul></div>';
  }
  html += '<div class="facet"><h3>Genomic context</h3><ul>' +
    [['plasmid', 'on a plasmid'], ['partner', 'verified operon partner'], ['regulator', 'divergent regulator']]
      .map(([k, l]) => `<li><a href="#" onclick="toggle('${k}');return false">${F[k] ? '<b>' + l + '</b>' : l}</a></li>`).join('') +
    '</ul></div>';
  document.getElementById('facets').innerHTML = html;
}
document.getElementById('q').addEventListener('keydown', e => { if (e.key === 'Enter') run(); });
</script></body></html>"""


if __name__ == "__main__":
    main()
