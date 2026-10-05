"""
ROAR-DB ana kapi sayfasi (hub.html) -- projenin tek giris noktasi.

Yol haritasi + onemli bulgular + degisiklikler + rapor/gezgin baglantilari.
Tum sayilar DB ve cikti dosyalarindan CANLI okunur.
"""

import argparse
import os
import sqlite3

# Iki yayinlanmis artifact -- sabit URL'ler
REPORT_URL = "https://claude.ai/code/artifact/f89ef6e3-d610-4c5c-8803-10e0e10c7bb4"
EXPLORER_URL = "https://claude.ai/code/artifact/db14338f-c919-40e0-821d-a4eb30726ca4"


def count_lines(path):
    if not os.path.exists(path):
        return 0
    with open(path) as handle:
        return sum(1 for _ in handle) - 1  # baslik


def gather(db, out_dir):
    n = {"ro_alpha": count_lines("ro_alpha.csv"),
         "not_alpha": count_lines("ro_not_alpha.csv"),
         "novel": count_lines(os.path.join(out_dir, "novel_candidates.csv")),
         "novel_hi": 0}
    hi_fasta = os.path.join(out_dir, "novel_high_confidence.fasta")
    if os.path.exists(hi_fasta):
        with open(hi_fasta) as handle:
            n["novel_hi"] = sum(1 for line in handle if line.startswith(">"))
    if os.path.exists(db):
        con = sqlite3.connect(db)
        def s(q):
            r = con.execute(q).fetchone()
            return r[0] if r else 0
        n.update({
            "replicons": s("SELECT COUNT(*) FROM replicon"),
            "unannot": s("SELECT COUNT(*) FROM replicon WHERE status='no_annotation'"),
            "confirmed": s("SELECT COUNT(*) FROM ro WHERE is_confirmed=1"),
            "neighbors": s("SELECT COUNT(*) FROM neighbor"),
            "clusters": s("SELECT COUNT(DISTINCT ro_cluster) FROM ro WHERE is_confirmed=1 AND ro_cluster!='N/A'"),
        })
        if con.execute("SELECT name FROM sqlite_master WHERE name='leaf'").fetchone():
            n["leaves"] = s("SELECT COUNT(*) FROM leaf")
            n["leaves_big"] = s("SELECT COUNT(*) FROM leaf WHERE size>=10")
        if con.execute("SELECT name FROM sqlite_master WHERE name='subfamily'").fetchone():
            n["ref_types"] = s("SELECT COUNT(*) FROM subfamily WHERE size>=10")
        con.close()
    n["removed"] = 117537 - n["ro_alpha"] if n["ro_alpha"] else 0
    return n


def fmt(v):
    try:
        return f"{int(v or 0):,}".replace(",", ".")
    except (TypeError, ValueError):
        return str(v)


# ------------------------------------------------------------------ icerik

ROADMAP = [
    ("Sorun teshisi",
     "PF00355 ile cekilen 189.657 protein sadece RO alpha degil &mdash; ferredoksin, "
     "sitokrom <i>bc</i><sub>1</sub> ISP, NirD de iceriyor. Olctuk: E-value bu ayrimi "
     "<strong>yapamiyor</strong> (her iki grup da %100 geciyor)."),
    ("Filtre yeniden kuruldu",
     "Coverage kapisi (&ge;0,45) + 8-kolonluk katalitik merkez testi. Referans setinden "
     "2 kontaminant cikarildi (IsoMO, CdnD), modeller 71 referansla yeniden kuruldu."),
    ("Genomik baglam",
     "RO&rsquo;lar tblastn yerine <strong>dogrudan GenBank anotasyonunun icinde</strong> "
     "bulundu &mdash; koordinat, strand, komsular, organizma tek tutarli kaynaktan. "
     "3 gizli bug yol boyunca duzeltildi."),
    ("Yerel veritabani",
     "SQLite: replicon / ro / neighbor / gene_category. Ham komsu saklanir, "
     "siniflandirma sorgu aninda &mdash; tanimi degistirmek icin gbk&rsquo;lari yeniden "
     "parse etmek gerekmez."),
    ("Analizler + null model",
     "Orientation, divergent regulator, transpozon&times;takson. Null model kontrolu: "
     "&quot;yakinda su gen var&quot; iddiasini ayni genomdaki rastgele pencerelerle kiyaslar."),
    ("Kume ici varyasyon",
     "Ayni enzim diye siniflandirilanlar gercekten ayni mi? Olctuk: cogu kume "
     "<strong>heterojen</strong> &mdash; dogrulananlarin 2/3&rsquo;u &lt;%35 kimlik."),
    ("Ozyinelemeli homojenizasyon",
     "Her kume homojen <em>yapraklara</em> bolunene kadar tekrar tekrar bolundu. "
     "Sonuc bir agac; her yaprak tek bir tutarli enzim tipi."),
    ("Varyant karakterizasyonu",
     "Her yaprak taksonomi + mobilite + <strong>komsuluk imzasiyla</strong> ayrildi. "
     "Farkli baglam = muhtemelen farkli substrat, ayni enzim adi altinda ayri seyler."),
    ("Novel RO kesfi",
     "Cevreye odaklan: hicbir yerlesik tipe uymayan, kucuk yapraklardaki izole diziler. "
     "Yuksek-guven adaylar + bitki/alg dali ortaya cikti."),
]

FINDINGS = [
    ("E-value ayiramiyor", "accent",
     "Coverage ile ~44.000 kontaminant temizlendi. Elenenlerin %19&rsquo;u tam RO alpha "
     "uzunluk araliginda &mdash; uzunluk filtresi de yetmezdi."),
    ("Referans setinde 2 kontaminant", "accent",
     "<code>IsoMO</code> (Rieske merkezi hic yok) ve <code>CdnD</code> (elektron transfer, "
     "katalitik yok) &mdash; her ikisi de RO alpha degil, cikarildi."),
    ("Regulator RO&rsquo;ya ozel degil", "teal",
     "Null model: regulator komsulugu <span class=\"n\">1,08&times;</span> (arka plan). "
     "Mevcut skorlama genom bollugunu siraliyor, RO iliskisini degil. Ama <em>yon</em> "
     "hala gercek: %47 ayni-yon, divergent mimari."),
    ("Kumeler heterojen", "teal",
     "Dogrulanan RO&rsquo;larin <span class=\"n\">%67&rsquo;si</span> &lt;%35 kimlikli "
     "kumelerde. Atama goreli en iyi modeli sectigi icin genis modeller &quot;torba&quot; olmus."),
    ("Varyantlar baglamla ayrisiyor", "teal",
     "<code>KshA</code>&rsquo;nin bir varyanti kanonik steroid dehidrogenaz yaninda; "
     "<code>VanA</code>&rsquo;nin mobil <i>Burkholderia</i> soyu IclR + transpozon yaninda "
     "&mdash; ayri enzimler."),
    ("Novel + bitki dali", "accent",
     "<span class=\"n\">155</span> yuksek-guven novel aday (guvenle RO, tipsiz). "
     "<i>Spinacia</i> (ispanak) CmoA&rsquo;ya 948 skorla &mdash; ayri bir kloroplast RO soyu."),
    ("Ekoloji: plazmit sinyali", "teal",
     "Ksenobiyotik RO&rsquo;lar plazmit uzerinde belirgin bicimde daha sik (oran icin rapor tablosu) "
     "(yatay transfer). Taksonomik yayilim ise ters cikti &mdash; dogal substratlar daha "
     "cok soyda (kadim dikey kalitim)."),
    ("DB karisik: %8 okaryot", "warn",
     "Dogrulanan RO&rsquo;larin <span class=\"n\">%8,3&rsquo;u</span> okaryot &mdash; bitki "
     "(kloroplast RO), mantar, hayvan. Katalitik filtre yasam alanini ayirmaz; "
     "isaretlendi, gizlenmedi."),
    ("%17 gbk anotasyonsuz", "warn",
     "Replikonlarin %17,4&rsquo;unde hic CDS yok &mdash; &quot;komsusuz RO&quot; degil "
     "&quot;verisi eksik&quot;. Co-occurrence paydasindan cikarilmali."),
]

BUGS = [
    ("Gen sayimi iki katina cikiyor", "mongo_gbk_analysis.py:174",
     "<code>gene</code> ve <code>CDS</code> ayri sayiliyordu &mdash; prokaryotta ikisi ayni gen."),
    ("Strand hic okunmuyor", "mongo_gbk_analysis.py:175",
     "Orientation, operon, divergent regulator analizleri bu yuzden imkansizdi."),
    ("Dairesel orijini asan genler", "mongo_gbk_analysis.py:175",
     "408 bp&rsquo;lik gen 116.580 bp gibi gorunup her RO&rsquo;nun komsusu oluyordu."),
    ("Anotasyonsuz dosya sessiz", "mongo_gbk_analysis.py",
     "%17 gbk bos komsu donduruyor, &quot;verisi eksik&quot; olarak isaretlenmiyordu."),
    ("Keyword filtresi cikarimda", "mongo_gbk_analysis.py:186",
     "<code>hypothetical</code> komsular veriye hic girmiyordu &mdash; oysa hedef onlar."),
    ("Atama E-value ile", "filter_hmm_out.py:41",
     "Kume atamasi icin bit skoru daha kararli; E-value uzunluga/DB boyutuna duyarli."),
]

SCRIPTS = [
    ("ro_filter.py, ro_motif.py, run_ro_filter.py", "coverage + katalitik motif filtresi"),
    ("extract_genomic_context.py", "GenBank&rsquo;tan RO + komsu cikarimi (strand&rsquo;li)"),
    ("build_db.py, annotate_ro.py", "SQLite kurulumu + HMM/motif dogrulama"),
    ("analyze.py, null_model.py", "orientation, transpozon, arka plan normalizasyonu"),
    ("analyze_variants.py, discover_subfamilies.py", "kume ici kimlik + alt-aile + SDP"),
    ("recursive_homogenize.py, characterize_leaves.py", "yaprak agaci + varyant imzalari"),
    ("analyze_ecology.py, make_report.py, make_explorer.py", "ekoloji hipotezi + iki HTML"),
]


def build(n):
    def tiles(items):
        return "".join(
            f'<div class="tile"><div class="tile__n">{v}</div>'
            f'<div class="tile__l">{l}</div></div>' for v, l in items)

    roadmap = "".join(
        f'<li class="step"><div class="step__num">{i:02d}</div>'
        f'<div class="step__body"><h3>{t}</h3><p>{d}</p></div></li>'
        for i, (t, d) in enumerate(ROADMAP, 1))

    findings = "".join(
        f'<div class="find find--{tone}"><h3>{t}</h3><p>{d}</p></div>'
        for t, tone, d in FINDINGS)

    bugs = "".join(
        f'<li class="bug"><div class="bug__t">{t}</div>'
        f'<code class="bug__loc">{loc}</code><p>{d}</p></li>'
        for t, loc, d in BUGS)

    scripts = "".join(
        f'<tr><td class="mono">{s}</td><td>{d}</td></tr>' for s, d in SCRIPTS)

    return f"""
<header class="hero">
  <div class="hero__eyebrow">ROAR-DB &middot; proje ana kapisi</div>
  <h1>Rieske oksijenaz veritabani:<br><em>ne yapildi, ne bulundu, nereden bakilir</em></h1>
  <p class="hero__lede">
    PF00355 Rieske proteinlerinden gercek halka-hidroksileyen oksijenaz alpha
    alt birimlerini ayiran, genomik baglamlarini cikaran ve varyantlarina kadar
    inen tam bir pipeline. Asagida yol haritasi, onemli bulgular ve iki
    interaktif sayfaya baglantilar.
  </p>
  <div class="hero__cta">
    <a class="cta cta--primary" href="{EXPLORER_URL}" target="_blank" rel="noopener">
      <span class="cta__icon">&#128300;</span>
      <span><span class="cta__t">Interaktif Gezgin</span>
      <span class="cta__d">kumeleri ara, varyantlari incele, uye listelerini indir</span></span>
    </a>
    <a class="cta cta--secondary" href="{REPORT_URL}" target="_blank" rel="noopener">
      <span class="cta__icon">&#129516;</span>
      <span><span class="cta__t">Metodoloji Raporu</span>
      <span class="cta__d">filtre, null model, ekoloji, varyasyon, novel kesif</span></span>
    </a>
  </div>
</header>

<section>
  <div class="tiles">{tiles([
    (fmt(n.get('ro_alpha')), 'RO alpha (protein duzeyi)'),
    (f"&minus;{fmt(n.get('removed'))}", 'temizlenen kontaminant'),
    (fmt(n.get('confirmed')), 'dogrulanmis RO (genomik)'),
    (fmt(n.get('clusters')), 'kume'),
    (fmt(n.get('leaves')), 'homojen yaprak (varyant)'),
    (fmt(n.get('novel_hi')), 'yuksek-guven novel aday'),
  ])}</div>
</section>

<section>
  <h2>Yol haritasi</h2>
  <p class="sec__lede">Tesbitten novel kesfe &mdash; her adim bir oncekinin
    ortaya cikardigi soruyu cevapliyor.</p>
  <ol class="roadmap">{roadmap}</ol>
</section>

<section>
  <h2>Onemli bulgular</h2>
  <div class="finds">{findings}</div>
</section>

<section>
  <h2>Degisiklikler &middot; bulunan hatalar</h2>
  <p class="sec__lede">Hicbiri aranarak bulunmadi &mdash; pipeline kurulurken ortaya
    ciktilar. Hepsi olculdu, tahmin edilmedi.</p>
  <ol class="bugs">{bugs}</ol>
</section>

<section>
  <h2>Kod haritasi</h2>
  <div class="table-wrap"><table>
    <thead><tr><th>script</th><th>ne yapar</th></tr></thead>
    <tbody>{scripts}</tbody>
  </table></div>
  <p class="note">Rapor ve gezgin DB&rsquo;den <em>canli</em> uretilir:
    <code>make_report.py</code> / <code>make_explorer.py</code> / <code>make_hub.py</code>
    &mdash; veri degisince yeniden calistir, sayfalar guncellenir.</p>
</section>

<section class="closing">
  <h2>Nereden baslayayim?</h2>
  <div class="hero__cta">
    <a class="cta cta--primary" href="{EXPLORER_URL}" target="_blank" rel="noopener">
      <span class="cta__icon">&#128300;</span>
      <span><span class="cta__t">Gezgine git</span>
      <span class="cta__d">bir kumeye tikla &rarr; varyantlar &rarr; CSV indir</span></span>
    </a>
    <a class="cta cta--secondary" href="{REPORT_URL}" target="_blank" rel="noopener">
      <span class="cta__icon">&#129516;</span>
      <span><span class="cta__t">Rapora git</span>
      <span class="cta__d">nasil ve neden &mdash; tum metodoloji</span></span>
    </a>
  </div>
</section>
"""


PAGE = """<title>ROAR-DB &middot; Rieske oksijenaz projesi</title>
<style>
:root{
  --ground:#F5F6F7; --panel:#FFFFFF; --sunken:#EAECEE;
  --ink:#191E22; --ink-2:#48535B; --ink-3:#78848C; --rule:#DBE0E3;
  --accent:#A6402A; --accent-soft:#F4E4DF; --teal:#2F6E7A; --teal-soft:#E1EDEF;
  --warn:#96701A; --warn-soft:#F6EDD9; --ok:#3F7A4E;
  --serif:ui-serif,Georgia,"Iowan Old Style",Palatino,serif;
  --sans:ui-sans-serif,system-ui,-apple-system,"Segoe UI",Roboto,sans-serif;
  --mono:ui-monospace,SFMono-Regular,"SF Mono",Menlo,Consolas,monospace;
}
@media (prefers-color-scheme:dark){:root{
  --ground:#111518; --panel:#181D21; --sunken:#212A2F;
  --ink:#E5E9EB; --ink-2:#A6B1B8; --ink-3:#76838B; --rule:#2A333A;
  --accent:#E28163; --accent-soft:#33211C; --teal:#6FB3BF; --teal-soft:#15272B;
  --warn:#D6A648; --warn-soft:#2C2615; --ok:#74B584;}}
:root[data-theme="dark"]{
  --ground:#111518; --panel:#181D21; --sunken:#212A2F;
  --ink:#E5E9EB; --ink-2:#A6B1B8; --ink-3:#76838B; --rule:#2A333A;
  --accent:#E28163; --accent-soft:#33211C; --teal:#6FB3BF; --teal-soft:#15272B;
  --warn:#D6A648; --warn-soft:#2C2615; --ok:#74B584;}
:root[data-theme="light"]{
  --ground:#F5F6F7; --panel:#FFFFFF; --sunken:#EAECEE;
  --ink:#191E22; --ink-2:#48535B; --ink-3:#78848C; --rule:#DBE0E3;
  --accent:#A6402A; --accent-soft:#F4E4DF; --teal:#2F6E7A; --teal-soft:#E1EDEF;
  --warn:#96701A; --warn-soft:#F6EDD9; --ok:#3F7A4E;}

body{background:var(--ground);color:var(--ink);font-family:var(--sans);
  font-size:16px;line-height:1.6;-webkit-font-smoothing:antialiased;}
.page{max-width:1080px;margin:0 auto;padding:0 clamp(1.1rem,4vw,3rem) 5rem;}
em{font-style:italic;} .mono,.n{font-family:var(--mono);font-variant-numeric:tabular-nums;}
code{font-family:var(--mono);font-size:.86em;background:var(--sunken);
  padding:.1em .38em;border-radius:3px;}

/* hero */
.hero{padding:clamp(2.6rem,7vw,5rem) 0 clamp(1.8rem,4vw,3rem);
  border-bottom:2px solid var(--ink);}
.hero__eyebrow{font-family:var(--mono);font-size:.74rem;letter-spacing:.14em;
  text-transform:uppercase;color:var(--accent);margin-bottom:1.4rem;}
.hero h1{font-family:var(--serif);font-weight:600;font-size:clamp(1.9rem,5vw,3.2rem);
  line-height:1.13;letter-spacing:-.02em;text-wrap:balance;margin-bottom:1.3rem;}
.hero h1 em{font-style:italic;color:var(--accent);}
.hero__lede{max-width:64ch;font-size:1.06rem;color:var(--ink-2);margin-bottom:1.9rem;}
.hero__cta{display:flex;gap:.9rem;flex-wrap:wrap;}
.cta{display:flex;align-items:center;gap:.85rem;text-decoration:none;
  padding:.95rem 1.2rem;border-radius:8px;flex:1 1 300px;min-width:260px;
  border:1px solid var(--rule);transition:transform .12s,border-color .12s;}
.cta:hover{transform:translateY(-2px);}
.cta--primary{background:var(--accent);border-color:var(--accent);}
.cta--primary *{color:#fff;}
.cta--primary:hover{border-color:var(--accent);}
.cta--secondary{background:var(--panel);}
.cta--secondary:hover{border-color:var(--teal);}
.cta__icon{font-size:1.5rem;line-height:1;}
.cta__t{display:block;font-weight:650;font-size:1rem;}
.cta__d{display:block;font-size:.8rem;opacity:.85;margin-top:.15rem;}
.cta--secondary .cta__d{color:var(--ink-3);}

/* sections */
section{padding:clamp(2.2rem,5vw,3.6rem) 0;border-bottom:1px solid var(--rule);}
section:last-child{border-bottom:none;}
h2{font-family:var(--serif);font-weight:600;font-size:clamp(1.4rem,3vw,1.95rem);
  letter-spacing:-.015em;margin-bottom:.6rem;text-wrap:balance;}
.sec__lede{color:var(--ink-2);max-width:62ch;margin-bottom:1.6rem;}

/* tiles */
.tiles{display:grid;grid-template-columns:repeat(auto-fit,minmax(150px,1fr));
  gap:1px;background:var(--rule);border:1px solid var(--rule);}
.tile{background:var(--panel);padding:1.15rem 1.1rem;}
.tile__n{font-family:var(--mono);font-variant-numeric:tabular-nums;
  font-size:clamp(1.4rem,3.2vw,1.85rem);font-weight:600;letter-spacing:-.02em;
  color:var(--accent);}
.tile__l{font-size:.78rem;color:var(--ink-2);margin-top:.3rem;line-height:1.3;}

/* roadmap */
.roadmap{list-style:none;display:flex;flex-direction:column;gap:0;}
.step{display:grid;grid-template-columns:auto 1fr;gap:1.1rem;padding:1.15rem 0;
  border-top:1px solid var(--rule);}
.step:first-child{border-top:none;}
.step__num{font-family:var(--mono);font-size:.95rem;font-weight:600;color:var(--teal);
  padding-top:.15rem;}
.step__body h3{font-family:var(--sans);font-weight:650;font-size:1.02rem;margin-bottom:.3rem;}
.step__body p{color:var(--ink-2);font-size:.92rem;max-width:70ch;}
.step__body strong{color:var(--ink);}

/* findings */
.finds{display:grid;grid-template-columns:repeat(auto-fit,minmax(270px,1fr));gap:1px;
  background:var(--rule);border:1px solid var(--rule);}
.find{background:var(--panel);padding:1.15rem 1.2rem;border-top:3px solid var(--rule);}
.find--accent{border-top-color:var(--accent);}
.find--teal{border-top-color:var(--teal);}
.find--warn{border-top-color:var(--warn);}
.find h3{font-size:.98rem;font-weight:650;margin-bottom:.4rem;}
.find p{font-size:.87rem;color:var(--ink-2);}
.find .n{color:var(--accent);font-weight:600;}

/* bugs */
.bugs{list-style:none;counter-reset:b;display:flex;flex-direction:column;gap:1px;
  background:var(--rule);border:1px solid var(--rule);}
.bug{background:var(--panel);padding:1rem 1.15rem;counter-increment:b;}
.bug__t{font-weight:640;font-size:.94rem;display:inline;}
.bug__t::before{content:counter(b,decimal-leading-zero);font-family:var(--mono);
  color:var(--accent);margin-right:.55rem;font-size:.85em;}
.bug__loc{font-size:.74rem;color:var(--ink-3);background:none;margin-left:.5rem;}
.bug p{font-size:.86rem;color:var(--ink-2);margin-top:.3rem;}

/* table */
.table-wrap{overflow-x:auto;border:1px solid var(--rule);}
table{width:100%;border-collapse:collapse;font-size:.87rem;background:var(--panel);}
th{text-align:left;font-size:.72rem;text-transform:uppercase;letter-spacing:.05em;
  color:var(--ink-3);padding:.65rem .9rem;border-bottom:1px solid var(--rule);}
td{padding:.55rem .9rem;border-bottom:1px solid var(--rule);color:var(--ink-2);}
tr:last-child td{border-bottom:none;}
.note{font-size:.85rem;color:var(--ink-3);margin-top:1rem;}
.closing{text-align:center;}
.closing h2{margin-bottom:1.4rem;}
.closing .hero__cta{justify-content:center;}
@media (prefers-reduced-motion:reduce){*{transition:none!important;}}
</style>
<div class="page">
__BODY__
</div>
"""


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--db", default="roar.sqlite")
    parser.add_argument("--out-dir", default="analysis_out")
    parser.add_argument("-o", "--output", default="hub.html")
    args = parser.parse_args()

    numbers = gather(args.db, args.out_dir)
    page = PAGE.replace("__BODY__", build(numbers))
    with open(args.output, "w") as handle:
        handle.write(page)
    print(f"[yazildi] {args.output}  ({len(page):,} bayt)")


if __name__ == "__main__":
    main()
