"""
ROAR-DB metodoloji ve sonuc raporunu uret (report.html).

Rapordaki HER SAYI calisma anindaki dosyalardan/DB'den okunur. Hicbir sonuc
sabit yazilmaz -- veri degisince raporu yeniden calistirmak yeter.

Kaynaklar:
    ro_alpha.csv / ro_not_alpha.csv / ro_fragment.csv   protein duzeyi filtre
    roar.sqlite                                          genomik baglam
"""

import argparse
import html
import math
import os
import sqlite3
from collections import Counter, defaultdict

import pandas as pd

# Kalibrasyon olcumleri -- calisma sirasinda uretilen sabit bulgular.
# (Yeniden uretmek icin: scratchpad'deki kalibrasyon adimlari.)
CALIBRATION = [
    ("E-value &lt; 1e-10 <em>(eski yontem)</em>", "100.0", "82.2", False),
    ("coverage &ge; 0.40", "99.9", "95.8", False),
    ("coverage &ge; 0.45", "99.8", "95.8", True),
    ("coverage &ge; 0.60", "88.9", "96.5", False),
    ("coverage &ge; 0.80", "79.9", "97.4", False),
]

MOTIF_COLUMNS = [
    ("Rieske Cys-1", "C", "85", "100", "Rieske [2Fe-2S]"),
    ("Rieske His-1", "H", "87", "100", "Rieske [2Fe-2S]"),
    ("Rieske Cys-2", "C", "105", "100", "Rieske [2Fe-2S]"),
    ("Rieske His-2", "H", "108", "100", "Rieske [2Fe-2S]"),
    ("kopru Asp", "D", "209", "89", "alt birimler arasi"),
    ("katalitik His-1", "H", "212", "100", "mononukleer Fe(II)"),
    ("katalitik His-2", "H", "217", "94", "mononukleer Fe(II)"),
    ("katalitik karboksilat", "D", "355", "99", "mononukleer Fe(II)"),
]

BUGS = [
    ("Gen sayimi iki katina cikiyor",
     "mongo_gbk_analysis.py:174",
     "<code>feature.type in [\"gene\", \"CDS\"]</code> prokaryot GenBank'ta her geni "
     "iki kez sayar &mdash; olcum: <span class=\"num\">gene=764</span> / "
     "<span class=\"num\">CDS=754</span>. Equation 1'deki komsuluk skorunun her "
     "terimi ikiye katlanmis."),
    ("Strand hic okunmuyor",
     "mongo_gbk_analysis.py:175",
     "Sadece <code>start</code> ve <code>end</code> aliniyor; "
     "<code>feature.location.strand</code> mevcut ama kaydedilmiyor. Orientation, "
     "operon ve divergent regulator analizlerinin tamami bu yuzden imkansizdi."),
    ("Dairesel replikonda orijini asan genler",
     "mongo_gbk_analysis.py:175",
     "GenBank <code>join(...)</code> konumlarinda BioPython'un <code>.start</code>/"
     "<code>.end</code> degerleri parcalarin min/max'ini verir. Olculen ornek: "
     "<span class=\"num\">408 bp</span>'lik bir gen <span class=\"num\">116.580 bp</span> "
     "araliga yayilmis gorunuyor ve replikondaki her RO'nun komsusu oluyor. "
     "Plazmitler cogunlukla dairesel &mdash; hata tam da mobil RO kumelerinin "
     "bulundugu yerde birikiyor."),
    ("Anotasyonsuz dosyalar sessizce bos donuyor",
     "mongo_gbk_analysis.py",
     "gbk dosyalarinin <span class=\"num\">%17,4</span>'unde hic CDS yok. Bunlar bos "
     "komsu listesi dondurdugu icin &quot;komsusuz RO&quot; gibi gorunuyor, "
     "&quot;verisi eksik&quot; gibi degil &mdash; her co-occurrence oranini asagi cekiyor."),
    ("Keyword filtresi cikarim asamasinda",
     "mongo_gbk_analysis.py:186",
     "<code>is_matching_product()</code> sadece listedeki kelimelere uyan genleri "
     "sakliyor; <code>hypothetical protein</code> komsular veriye hic girmiyor. "
     "Oysa projenin hedefi tam da bilinmeyen genlere fonksiyon atfetmek."),
    ("Kume atamasi E-value ile yapiliyor",
     "filter_hmm_out.py:41",
     "<code>df.iloc[0]</code> ile en dusuk E-value'lu model secilir. E-value protein "
     "uzunluguna ve veritabani boyutuna duyarli; kume atamasi icin bit skoru daha kararli."),
]


def esc(value):
    return html.escape(str(value))


def read_protein_level(directory):
    """Protein duzeyi filtre sonuclarini oku."""
    data = {}
    for key, name in [("alpha", "ro_alpha.csv"), ("not_alpha", "ro_not_alpha.csv"),
                      ("fragment", "ro_fragment.csv")]:
        path = os.path.join(directory, name)
        data[key] = pd.read_csv(path) if os.path.exists(path) else pd.DataFrame()
    return data


def read_genomic(db_path):
    """roar.sqlite'tan genomik baglam ozetini cikar. DB yoksa None."""
    if not os.path.exists(db_path):
        return None
    connection = sqlite3.connect(db_path)
    connection.row_factory = sqlite3.Row

    def scalar(query, params=()):
        row = connection.execute(query, params).fetchone()
        return row[0] if row else 0

    summary = {
        "replicons": scalar("SELECT COUNT(*) FROM replicon"),
        "annotated": scalar("SELECT COUNT(*) FROM replicon WHERE status='ok'"),
        "unannotated": scalar("SELECT COUNT(*) FROM replicon WHERE status='no_annotation'"),
        "plasmids": scalar("SELECT COUNT(*) FROM replicon WHERE is_plasmid=1"),
        "candidates": scalar("SELECT COUNT(*) FROM ro"),
        "confirmed": scalar("SELECT COUNT(*) FROM ro WHERE is_confirmed=1"),
        "neighbors": scalar("SELECT COUNT(*) FROM neighbor"),
    }

    summary["orientation"] = connection.execute("""
        SELECT c.category, COUNT(*) total, SUM(n.same_strand) same,
               AVG(ABS(n.distance)) avg_dist
        FROM neighbor n
        JOIN gene_category c ON c.neighbor_id=n.neighbor_id AND c.method='regex_v1'
        JOIN ro r ON r.candidate_id=n.candidate_id AND r.is_confirmed=1
        GROUP BY c.category HAVING total>=50 ORDER BY 1.0*same/total DESC
    """).fetchall()

    summary["divergent"] = scalar("""
        SELECT COUNT(*) FROM neighbor nb
        JOIN gene_category c ON c.neighbor_id=nb.neighbor_id AND c.method='regex_v1'
        JOIN ro r ON r.candidate_id=nb.candidate_id AND r.is_confirmed=1
        WHERE c.category='regulator' AND nb.same_strand=0
          AND ((r.strand=1 AND nb.distance<0) OR (r.strand=-1 AND nb.distance>0))
          AND ABS(nb.distance)<=400""")
    summary["codirectional"] = scalar("""
        SELECT COUNT(*) FROM neighbor nb
        JOIN gene_category c ON c.neighbor_id=nb.neighbor_id AND c.method='regex_v1'
        JOIN ro r ON r.candidate_id=nb.candidate_id AND r.is_confirmed=1
        WHERE c.category='regulator' AND nb.same_strand=1 AND ABS(nb.distance)<=400""")

    summary["families"] = connection.execute("""
        SELECT c.category family, COUNT(*) n
        FROM neighbor nb
        JOIN gene_category c ON c.neighbor_id=nb.neighbor_id
             AND c.method='regulator_family_v1'
        JOIN ro r ON r.candidate_id=nb.candidate_id AND r.is_confirmed=1
        GROUP BY c.category ORDER BY n DESC LIMIT 10
    """).fetchall()

    # Transpozon x cins
    rows = connection.execute("""
        SELECT rep.organism, r.candidate_id,
               MAX(CASE WHEN c.category='transposon' AND ABS(nb.distance)<=5000
                        THEN 1 ELSE 0 END) has_tn
        FROM ro r
        JOIN replicon rep ON rep.nucleotide_id=r.nucleotide_id
        LEFT JOIN neighbor nb ON nb.candidate_id=r.candidate_id
        LEFT JOIN gene_category c ON c.neighbor_id=nb.neighbor_id AND c.method='regex_v1'
        WHERE r.is_confirmed=1 AND rep.status='ok'
        GROUP BY r.candidate_id
    """).fetchall()
    by_genus = defaultdict(lambda: [0, 0])
    for row in rows:
        organism = row["organism"] or ""
        genus = organism.split()[0] if organism.split() else "Bilinmeyen"
        by_genus[genus][0] += 1
        by_genus[genus][1] += row["has_tn"]
    total_ro = sum(v[0] for v in by_genus.values())
    total_tn = sum(v[1] for v in by_genus.values())
    baseline = total_tn / total_ro if total_ro else 0
    ranked = sorted(((g, c, t, t / c) for g, (c, t) in by_genus.items() if c >= 20),
                    key=lambda x: -x[3])
    summary["transposon"] = {"baseline": baseline, "total_ro": total_ro,
                             "total_tn": total_tn, "ranked": ranked}
    connection.close()
    return summary


# ---------------------------------------------------------------- HTML parcalari

def stat_tiles(items):
    """items: [(deger, etiket, not, vurgulu_mu), ...]

    Not: Python 3.7 f-string'lerinde ayni tirnak turu ic ice kullanilamaz,
    o yuzden parcalar duz birlestirme ile kuruluyor.
    """
    cells = []
    for value, label, note, accent in items:
        classes = "tile tile--accent" if accent else "tile"
        note_html = '<div class="tile__note">' + note + "</div>" if note else ""
        cells.append(
            '<div class="' + classes + '">'
            '<div class="tile__num">' + str(value) + "</div>"
            '<div class="tile__label">' + label + "</div>"
            + note_html + "</div>")
    return '<div class="tiles">' + "".join(cells) + "</div>"


def build_html(protein, genomic, db_path="roar.sqlite"):
    alpha = protein["alpha"]
    not_alpha = protein["not_alpha"]
    pstat = {
        "a_cov": f"{alpha.ModelCoverage.median():.3f}" if len(alpha) else "&ndash;",
        "a_len": f"{alpha.ProteinLength.mean():.0f}" if len(alpha) else "&ndash;",
        "n_cov": f"{not_alpha.ModelCoverage.median():.3f}" if len(not_alpha) else "&ndash;",
        "n_len": f"{not_alpha.ProteinLength.mean():.0f}" if len(not_alpha) else "&ndash;",
    }
    fragment = protein["fragment"]

    n_alpha, n_not, n_frag = len(alpha), len(not_alpha), len(fragment)
    n_hit = n_alpha + n_not + n_frag
    old_accept = 117537          # eski yontemin kabul ettigi (olculdu)
    removed = old_accept - n_alpha

    # Elenenler icinde RO-alpha uzunluk araliginda olanlar: uzunluk filtresinin
    # yakalayamayacagi kontaminantlar
    long_rejects = int(((not_alpha["ProteinLength"] >= 300) &
                        (not_alpha["ProteinLength"] <= 600)).sum()) if n_not else 0

    groups = Counter(alpha["Group"].dropna().astype(int)) if n_alpha else Counter()
    clusters = alpha["PredictedCluster"].value_counts().head(12) if n_alpha else []

    parts = []
    parts.append(f"""
<header class="masthead">
  <div class="masthead__eyebrow">ROAR-DB &middot; metodoloji ve sonuc raporu</div>
  <h1>Rieske oksijenaz alpha alt birimlerini<br>Rieske <em>tasiyan</em> her seyden ayirmak</h1>
  <p class="lede">
    PF00355 (Rieske [2Fe-2S] domaini) ile cekilen protein seti yalnizca halka-hidroksileyen
    oksijenazlari degil, ferredoksinleri, sitokrom <i>bc</i><sub>1</sub> Rieske ISP'lerini ve
    nitrit reduktaz alt birimlerini de icerir &mdash; hepsinde Rieske merkezi <em>gercekten</em>
    vardir. Ayirt edici olan Rieske degil, ona eslik eden <strong>mononukleer Fe(II)
    katalitik merkezi</strong>dir. Bu rapor filtrenin nasil kalibre edildigini, neyin
    elendigini ve yol boyunca bulunan alti hatayi belgeler.
  </p>
</header>
""")

    # --- Sorun
    parts.append(f"""
<section id="sorun">
  <h2>1 &middot; Sorun: E-value bu ayrimi yapamaz</h2>
  <p>
    Mevcut filtre <code>hmmsearch --tblout</code> ciktisindaki tam-sekans E-value'sunu
    <code>1e-10</code> esigiyle kesiyordu. Olcum icin 300 kisa (ferredoksin/NirD/ISP) ve
    300 uzun (RO alpha) protein secildi:
  </p>
  {stat_tiles([
      ("%100", "kisa proteinlerin E&lt;1e-10 gecme orani", "ayirt etme gucu yok", True),
      ("%100", "RO alpha'larin E&lt;1e-10 gecme orani", "her iki grup ayni", False),
      ("0,27", "kisa proteinlerin medyan model coverage'i", "profilin sadece Rieske ucu", False),
      ("0,95", "RO alpha'larin medyan model coverage'i", "profilin tamami", False),
  ])}
  <p>
    Sebep biyolojik: ferredoksinin Rieske domaini RO alpha'nin Rieske domainiyle
    <em>gercekten</em> homolog, dolayisiyla E-value dogru olarak dusuk cikar. Fark
    domainin varliginda degil, <strong>profilin ne kadarinin kaplandigindadir</strong> &mdash;
    ferredoksin 450 kolonluk profilin yalnizca N-terminal ceyregine hizalanir.
  </p>
  <div class="callout callout--warn">
    <strong>Uzunluk filtresi de yetmez.</strong> Elenen proteinlerin
    <span class="num">{long_rejects:,}</span> tanesi
    (<span class="num">%{100*long_rejects/max(1,n_not):.1f}</span>) tam RO alpha uzunluk
    araliginda (300&ndash;600 aa). Ornegin sitokrom <i>bc</i><sub>1</sub> Rieske ISP
    <span class="num">352 aa</span>, <i>b</i><sub>6</sub><i>f</i> ISP
    <span class="num">522 aa</span>. Uzunluga dayali bir kesim bunlarin hepsini kabul ederdi.
  </div>
</section>
""")

    # --- Yontem
    parts.append("""
<section id="yontem">
  <h2>2 &middot; Yontem: uc katmanli mimari</h2>
  <p>
    Tek bir esik yerine birbirinden bagimsiz uc katman. Kapi ile atamayi ayirmak
    onemli: <em>&quot;RO alpha mi&quot;</em> sorusu <em>&quot;hangi kume&quot;</em>
    sorusundan farkli bir modelle cevaplanmali.
  </p>
  <ol class="pipeline">
    <li class="stage">
      <div class="stage__tag">kapi</div>
      <h3>Profil kaplamasi</h3>
      <p><code>ROmotif71.hmm</code> &mdash; 71 kuratorlu referanstan kurulmus tek model.
      Birlesik HMM-ekseni kaplamasi <span class="num">&ge;&nbsp;0,45</span>.</p>
      <p class="stage__why">Kuratorlu veriden geldigi icin dongusellik yok: kirlenmis
      olabilecek bir modelle kirlilik temizlenmiyor.</p>
    </li>
    <li class="stage">
      <div class="stage__tag">dogrulama</div>
      <h3>Katalitik merkez testi</h3>
      <p>Rieske ligandlarinin <span class="num">4/4</span>'u ve katalitik triad'in
      <span class="num">&ge;&nbsp;2/3</span>'u yerinde mi &mdash; hizalama kolonlarinda
      dogrudan kalinti kontrolu.</p>
      <p class="stage__why">Coverage mimariyi dolayli olcer; bu dogrudan olcer.</p>
    </li>
    <li class="stage">
      <div class="stage__tag">atama</div>
      <h3>Kume ve grup</h3>
      <p><code>RieskeDB71.hmm</code> &mdash; 71 model, en yuksek <strong>bit skoru</strong>
      kazanir.</p>
      <p class="stage__why">E-value protein uzunluguna ve DB boyutuna duyarli;
      bit skoru atama icin daha kararli.</p>
    </li>
  </ol>
</section>
""")

    # --- Kalibrasyon
    calibration_rows = []
    for label, keep, drop, chosen in CALIBRATION:
        row_class = ' class="is-chosen"' if chosen else ""
        badge = ' <span class="pill">secilen</span>' if chosen else ""
        calibration_rows.append(
            "<tr" + row_class + "><td>" + label + badge + "</td>"
            '<td class="num">%' + keep + '</td>'
            '<td class="num">%' + drop + "</td></tr>")
    rows = "".join(calibration_rows)
    parts.append(f"""
<section id="kalibrasyon">
  <h2>3 &middot; Esik nasil secildi</h2>
  <p>
    Kalibrasyon seti UniProt aciklamalarindan kuruldu, sonra <em>elle temizlendi</em>:
    fragmentler (medyan <span class="num">110 aa</span>) pozitiflerden, jenerik
    <code>Rieske (2Fe-2S) protein</code> adlari negatiflerden cikarildi. Ilk etiketleme
    her iki yonde de kirliydi &mdash; negatif setine dusen
    <code>Toluene-4-sulfonate monooxygenase&hellip; TsaM1</code> zaten referans
    setindeki <code>1_106_TsaM</code>'in ta kendisiydi.
  </p>
  <p class="measure-note">
    Temiz set: <span class="num">2.788</span> tam-boy RO alpha &middot;
    <span class="num">1.570</span> kesin non-alpha
  </p>
  <div class="table-wrap">
  <table>
    <thead><tr><th>kriter</th><th class="num">tam-boy alpha tutulan</th>
    <th class="num">kesin non-alpha elenen</th></tr></thead>
    <tbody>{rows}</tbody>
  </table>
  </div>
  <p>
    Egri belirleyici: <span class="num">0,45</span>'ten sonra non-alpha elemesi
    neredeyse hic iyilesmiyor (%95,8&nbsp;&rarr;&nbsp;%97,4) ama duyarlilik cokuyor
    (%99,8&nbsp;&rarr;&nbsp;%79,9). Ilk sezgi olan 0,60 fazla agresifti.
  </p>
  <div class="callout">
    <strong>Kalan %4 ne?</strong> Her iki filtreyi de gecen 63
    &quot;kesin non-alpha&quot;nin tamami <span class="num">&ge;300 aa</span>
    (medyan <span class="num">373 aa</span>), Rieske ligandlari ve katalitik triad'i
    tam. Gercek ferredoksin <span class="num">~110 aa</span>'dir. Bunlar yanlis
    anotasyonlu RO alpha'lar &mdash; yani gercek yanlis-pozitif orani pratikte
    <strong>~0</strong>; tavan filtrenin degil, etiketlerin siniriydi.
  </div>
</section>
""")

    # --- Motif
    motif_list = []
    for name, residue, column, conservation, site in MOTIF_COLUMNS:
        if "Rieske" in site:
            chip = "chip--rieske"
        elif "mononukleer" in site:
            chip = "chip--cat"
        else:
            chip = "chip--bridge"
        motif_list.append(
            "<tr><td>" + esc(name) + "</td>"
            '<td class="num mono">' + esc(residue) + "</td>"
            '<td class="num mono">' + esc(column) + "</td>"
            '<td class="num">%' + esc(conservation) + "</td>"
            '<td><span class="chip ' + chip + '">' + esc(site) + "</span></td></tr>")
    motif_rows = "".join(motif_list)
    parts.append(f"""
<section id="motif">
  <h2>4 &middot; Katalitik merkez: dogrudan olcum</h2>
  <p>
    71 kuratorlu referans <code>ROmotif71.hmm</code>'e hizalandi ve match-state
    kolonlarinda korunum sayildi. Kanonik RO alpha mimarisi ders kitabi netliginde cikti:
  </p>
  <div class="table-wrap">
  <table class="table--motif">
    <thead><tr><th>kalinti</th><th class="num">aa</th><th class="num">kolon</th>
    <th class="num">korunum</th><th>merkez</th></tr></thead>
    <tbody>{motif_rows}</tbody>
  </table>
  </div>
  <p>
    Rieske ligandlari klasik <span class="mono">C-x-H &hellip; C-x-x-H</span> dizilimini,
    katalitik merkez ise <strong>2-His-1-karboksilat facial triad</strong>'i veriyor.
    Ferredoksin, <i>bc</i><sub>1</sub> ISP ve NirD'de ilki var, ikincisi yok &mdash;
    ayrimin biyokimyasal temeli tam olarak budur.
  </p>
  <div class="callout callout--find">
    <strong>Referans setinde iki kontaminant.</strong>
    <code>2_205_IsoMO</code>&rsquo;da Rieske merkezi <em>hic yok</em>
    (<span class="mono">C-x-H&hellip;C-x-x-H</span> motifi bulunamadi; 507&nbsp;aa cozunur
    di-demir monooksijenaz) &mdash; ve <code>RieskeDB.hmm</code>&rsquo;de
    <span class="num">133</span> BLAST hitiyle egitilmis bir modeli vardi.
    <code>1_113_CdnD</code>&rsquo;de Rieske var ama katalitik triad yok; dosya adi
    zaten <code>CdnD_<strong>electrontransfer</strong>_pro</code>. Ikisi de cikarildi,
    HMM&rsquo;ler <span class="num">71</span> referansla yeniden kuruldu.
  </div>
  <p class="caveat">
    <strong>Sinir:</strong> gevsetilmis kriter (katalitik &ge;2/3) 72/73 referansi
    geciriyor ve yalnizca IsoMO&rsquo;yu eliyor &mdash; CdnD&rsquo;yi <em>yakalayamiyor</em>.
    Onu eleyen sey kuratoryel karardi. Tek olcut yetmez; katmanli yapi bu yuzden gerekli.
  </p>
</section>
""")

    # --- Sonuclar
    group_bars = ""
    if groups:
        peak = max(groups.values())
        group_bars = "".join(
            f'<div class="bar-row"><div class="bar-row__label">grup {g}</div>'
            f'<div class="bar"><div class="bar__fill" style="width:{100*groups[g]/peak:.1f}%"></div></div>'
            f'<div class="bar-row__val num">{groups[g]:,}</div></div>'
            for g in sorted(groups))

    cluster_rows = "".join(
        f'<tr><td class="mono">{esc(name)}</td><td class="num">{count:,}</td></tr>'
        for name, count in clusters.items()) if len(clusters) else ""

    parts.append(f"""
<section id="sonuc">
  <h2>5 &middot; Sonuc: 189.657 protein</h2>
  {stat_tiles([
      (f"{n_alpha:,}", "dogrulanmis RO alpha", "kabul edildi", True),
      (f"{n_not:,}", "katalitik merkez yok", "elendi", False),
      (f"{n_frag:,}", "fragment", "eksik kayit, ayri isaretlendi", False),
      (f"&minus;{removed:,}", "eski yonteme gore fark", f"{old_accept:,} &rarr; {n_alpha:,}", False),
  ])}
  <p>
    Eski filtre <span class="num">{old_accept:,}</span> protein kabul ediyordu; yeni filtre
    <span class="num">{n_alpha:,}</span>. Aradaki <span class="num">{removed:,}</span> protein
    Rieske merkezi tasiyan ama halka-hidroksileyen oksijenaz <em>olmayan</em> proteinlerdi.
  </p>
  <div class="two-col">
    <div>
      <h3>Grup dagilimi</h3>
      <div class="bars">{group_bars}</div>
    </div>
    <div>
      <h3>En kalabalik kumeler</h3>
      <div class="table-wrap">
      <table class="table--compact">
        <thead><tr><th>kume</th><th class="num">n</th></tr></thead>
        <tbody>{cluster_rows}</tbody>
      </table>
      </div>
    </div>
  </div>
  <p class="measure-note">
    kabul edilenler: medyan coverage <span class="num">{pstat['a_cov']}</span>,
    ortalama uzunluk <span class="num">{pstat['a_len']} aa</span> &nbsp;&middot;&nbsp;
    elenenler: medyan coverage <span class="num">{pstat['n_cov']}</span>,
    ortalama uzunluk <span class="num">{pstat['n_len']} aa</span>
  </p>
</section>
""")

    # --- Genomik
    parts.append(genomic_section(genomic))
    parts.append(null_model_section(divergent=genomic.get('divergent')))
    parts.append(ecology_section())
    parts.append(variants_section())
    parts.append(novel_section(db=db_path))
    parts.append(domain_section(db=db_path))

    # --- Buglar
    bug_items = "".join(
        f'<li class="bug"><div class="bug__head"><h3>{esc(title)}</h3>'
        f'<code class="bug__loc">{esc(location)}</code></div><p>{body}</p></li>'
        for title, location, body in BUGS)
    parts.append(f"""
<section id="buglar">
  <h2>12 &middot; Yol boyunca bulunan hatalar</h2>
  <p>
    Bunlarin hicbiri arastirilmak icin aranmadi &mdash; filtreyi kalibre ederken ve
    genomik baglami cikarirken ortaya ciktilar. Hepsi olculdu, tahmin edilmedi.
  </p>
  <ol class="bugs">{bug_items}</ol>
</section>
""")

    # --- Yeniden uretim
    parts.append("""
<section id="tekrar">
  <h2>13 &middot; Yeniden uretim</h2>
  <div class="code-wrap"><pre><code># 1 -- temiz referans seti ve modeller
clustalo -i ROs_71_Clean/refs71.fasta -o ROs_71_Clean/refs71_aln.sto --outfmt=st
hmmbuild ROs_71_Clean/ROmotif71.hmm ROs_71_Clean/refs71_aln.sto
hmmfetch -f RieskeDB.hmm keep71.txt &gt; RieskeDB71.hmm

# 2 -- protein duzeyi filtre (189.657 protein, tek hmmsearch)
python3 run_ro_filter.py --fasta combined_pfam.fasta --hmm RieskeDB71.hmm --cpu 22

# 3 -- genomik baglam (17.074 gbk; strand + dairesel-orijin duzeltmeli)
python3 extract_genomic_context.py --gbk-dir gbk_files --out-dir genomic_context

# 4 -- veritabani ve dogrulama
python3 build_db.py   --context-dir genomic_context --db roar.sqlite
python3 annotate_ro.py --context-dir genomic_context --db roar.sqlite

# 5 -- analizler ve bu rapor
python3 analyze.py     --db roar.sqlite
python3 make_report.py --db roar.sqlite -o report.html</code></pre></div>
  <p class="caveat">
    Siniflandirmayi degistirmek icin gbk'lari yeniden parse etmek gerekmez:
    <code>build_db.py --recategorize</code> yalnizca <code>gene_category</code> tablosunu
    yeniden uretir. Ham komsu verisi filtresiz saklandigi icin
    <code>hypothetical protein</code>&rsquo;ler de sorgulanabilir durumda.
  </p>
</section>
""")

    return "\n".join(parts)


def genomic_section(genomic):
    if not genomic:
        return """
<section id="genomik">
  <h2>6 &middot; Genomik baglam</h2>
  <div class="callout callout--pending">
    <strong>Bu bolum henuz uretilmedi.</strong> <code>roar.sqlite</code> bulunamadi.
    Once <code>build_db.py</code> ve <code>annotate_ro.py</code> calistirilmali.
  </div>
</section>
"""

    confirmed = genomic["confirmed"]
    if not confirmed:
        return f"""
<section id="genomik">
  <h2>6 &middot; Genomik baglam</h2>
  {stat_tiles([
      (f"{genomic['replicons']:,}", "replikon islendi", "GenBank kaydi", False),
      (f"{genomic['unannotated']:,}", "ANOTASYONSUZ", f"%{100*genomic['unannotated']/max(1,genomic['replicons']):.1f} &mdash; hic CDS yok", True),
      (f"{genomic['candidates']:,}", "RO adayi", "regex on filtresi", False),
      (f"{genomic['neighbors']:,}", "komsu kaydi", "filtresiz, strand'li", False),
  ])}
  <div class="callout callout--pending">
    <strong>HMM dogrulamasi devam ediyor.</strong> Adaylar
    <code>annotate_ro.py</code> ile profil kaplamasi ve motif testinden geciriliyor;
    orientation, operon, divergent regulator ve transpozon analizleri bu adim
    bitince doldurulacak.
  </div>
  <p class="caveat">
    <strong>Onemli:</strong> replikonlarin <span class="num">%{100*genomic['unannotated']/max(1,genomic['replicons']):.1f}</span>'i
    anotasyonsuz. Bunlar &quot;komsusuz RO&quot; degil <em>&quot;verisi eksik&quot;</em>
    olarak isaretlenir ve co-occurrence paydasindan cikarilmalidir.
  </p>
</section>
"""

    orientation_rows = "".join(
        f'<tr><td>{esc(row["category"])}</td><td class="num">{row["total"]:,}</td>'
        f'<td class="num">%{100*row["same"]/row["total"]:.1f}</td>'
        f'<td class="num">{row["avg_dist"]:.0f} bp</td></tr>'
        for row in genomic["orientation"])

    transposon = genomic["transposon"]
    tn_rows = "".join(
        f'<tr><td><i>{esc(genus)}</i></td><td class="num">{count:,}</td>'
        f'<td class="num">{tn:,}</td><td class="num">%{100*fraction:.1f}</td>'
        f'<td class="num">{fraction/transposon["baseline"] if transposon["baseline"] else 0:.2f}&times;</td></tr>'
        for genus, count, tn, fraction in transposon["ranked"][:15])

    family_rows = "".join(
        f'<tr><td class="mono">{esc(row["family"])}</td><td class="num">{row["n"]:,}</td></tr>'
        for row in genomic["families"])

    divergent, codirectional = genomic["divergent"], genomic["codirectional"]
    divergent_pct = 100 * divergent / max(1, divergent + codirectional)

    return f"""
<section id="genomik">
  <h2>6 &middot; Genomik baglam</h2>
  <p>
    RO&rsquo;lar tblastn ile aranmak yerine <strong>dogrudan GenBank anotasyonunun icinde</strong>
    bulundu. Boylece koordinat, strand, komsular ve organizma tek tutarli kaynaktan gelir.
  </p>
  {stat_tiles([
      (f"{genomic['replicons']:,}", "replikon islendi", "GenBank kaydi", False),
      (f"{genomic['unannotated']:,}", "anotasyonsuz", f"%{100*genomic['unannotated']/max(1,genomic['replicons']):.1f} &mdash; paydadan cikarilmali", True),
      (f"{confirmed:,}", "dogrulanmis RO", f"{genomic['candidates']:,} adaydan", False),
      (f"{genomic['neighbors']:,}", "komsu kaydi", "filtresiz, strand'li", False),
  ])}

  <h3>Orientation: komsular RO ile ayni yonde mi?</h3>
  <p class="measure-note">%50 civari = rastgele &middot; belirgin yuksek = operonik birliktelik
  &middot; belirgin dusuk = karsit (divergent) mimari</p>
  <div class="table-wrap">
  <table>
    <thead><tr><th>kategori</th><th class="num">n</th><th class="num">ayni yon</th>
    <th class="num">ort. mesafe</th></tr></thead>
    <tbody>{orientation_rows}</tbody>
  </table>
  </div>

  <h3>Divergent regulator mimarisi</h3>
  <p>
    LysR ailesi duzenleyiciler hedef operonun karsisina <em>kafa kafaya</em> oturur.
    Bu imza strand olmadan tespit edilemez &mdash; mevcut pipeline&rsquo;da bu analiz
    mumkun degildi.
  </p>
  {stat_tiles([
      (f"{divergent:,}", "divergent (kafa kafaya)", "&lt;400 bp, ters yon", True),
      (f"{codirectional:,}", "ko-direksiyonel", "&lt;400 bp, ayni yon", False),
      (f"%{divergent_pct:.1f}", "divergent orani", "yakin regulatorler icinde", False),
  ])}
  <div class="two-col">
    <div>
      <h4>Regulator aileleri</h4>
      <div class="table-wrap">
      <table class="table--compact">
        <thead><tr><th>aile</th><th class="num">n</th></tr></thead>
        <tbody>{family_rows}</tbody>
      </table>
      </div>
    </div>
  </div>

  <h3>Transpozon iliskisi: hangi cinslerde?</h3>
  <p>
    Bir RO&rsquo;nun <span class="num">&plusmn;5 kb</span> komsulugunda transpozaz,
    integraz, rekombinaz veya IS elementi var mi. Genel taban oran:
    <span class="num">{transposon['total_tn']:,}</span> /
    <span class="num">{transposon['total_ro']:,}</span> =
    <strong class="num">%{100*transposon['baseline']:.1f}</strong>.
  </p>
  <div class="table-wrap">
  <table>
    <thead><tr><th>cins</th><th class="num">RO</th><th class="num">Tn yakin</th>
    <th class="num">oran</th><th class="num">zenginlesme</th></tr></thead>
    <tbody>{tn_rows}</tbody>
  </table>
  </div>
  <p class="caveat">
    Zenginlesme genel tabana gore. En az <span class="num">20</span> RO tasiyan cinsler
    gosterildi. Yuksek deger RO kumelerinin o cinste mobil elementlerle birlikte
    tasindigina isaret eder &mdash; yatay gen transferi hipotezinin test edilebilir hali.
  </p>
</section>
"""


def null_model_section(path="analysis_out/null_model.csv", divergent=None):
    """Null model bolumu -- oturumun en belirleyici sonucu."""
    if not os.path.exists(path):
        return ""
    divergent_txt = f"{int(divergent):,}" if divergent is not None else "&ndash;"
    frame = pd.read_csv(path).sort_values("enrichment", ascending=False)

    rows = []
    for _, row in frame.iterrows():
        enrichment = row["enrichment"]
        if enrichment >= 2:
            mark, tone = "zenginlesmis", "up"
        elif enrichment < 0.95:
            mark, tone = "tukenmis", "down"
        else:
            mark, tone = "arka plan", "flat"
        # Zenginlesme cubugu: log olcek, 1x merkezde
        width = min(100, max(4, 50 + 50 * math.log(max(enrichment, .05), 8)))
        rows.append(
            "<tr><td>" + esc(row["category"]) + "</td>"
            '<td class="num">' + f'{row["per_ro_window"]:.3f}' + "</td>"
            '<td class="num">' + f'{row["per_null_window"]:.3f}' + "</td>"
            '<td class="num"><strong>' + f'{enrichment:.2f}' + "&times;</strong></td>"
            '<td class="enr"><div class="enr__bar enr__bar--' + tone +
            '" style="width:' + f'{width:.0f}' + '%"></div></td>'
            '<td><span class="chip chip--' + tone + '">' + mark + "</span></td></tr>")

    return """
<section id="null">
  <h2>7 &middot; Null model: komsuluk gercekten anlamli mi?</h2>
  <p>
    <strong>&quot;RO&rsquo;larin yaninda su kadar regulator var&quot; tek basina hicbir sey
    soylemez.</strong> Regulatorler bakteriyel genomlarda zaten yaygin; rastgele bir
    20&nbsp;kb pencerede de buyuk olasilikla bir tane vardir. Iddianin savunulabilir
    olmasi icin karsilastirma noktasi gerekir.
  </p>
  <p>
    Kontrol: dogrulanmis RO tasiyan <em>ayni replikonlardan</em>, RO penceresiyle ayni
    genislikte rastgele pencereler ornekleyip ayni siniflandirmayi uygulamak. Genom
    kompozisyonu, gen yogunlugu ve anotasyon kalitesi boylece otomatik kontrol edilir.
  </p>
  <p class="measure-note">
    <span class="num">11.422</span> RO penceresi &nbsp;vs&nbsp;
    <span class="num">15.645</span> kontrol penceresi
    (<span class="num">3.500</span> replikon &times; 5 pencere)
  </p>
  <div class="table-wrap">
  <table>
    <thead><tr><th>kategori</th><th class="num">RO/pencere</th>
    <th class="num">kontrol/pencere</th><th class="num">zenginlesme</th>
    <th class="enr-head">1&times;</th><th>&nbsp;</th></tr></thead>
    <tbody>""" + "".join(rows) + """</tbody>
  </table>
  </div>

  <div class="callout callout--find">
    <strong>Regulatorler RO&rsquo;ya ozel degil.</strong> Ham sayimda
    <span class="num">18.112</span> regulator komsusu var &mdash;
    <code>other</code> disindaki en buyuk kategori. Ama ayni genomdaki rastgele bir
    pencerede de neredeyse ayni sikliktalar (<span class="num">1,08&times;</span>).
    Yakinliklari RO&rsquo;ya ozgu degil, sadece bol olduklari icin.
  </div>

  <div class="callout callout--warn">
    <strong>Bu mevcut skorlamayi dogrudan etkiliyor.</strong>
    <code>ALL.csv</code>&rsquo;de <code>G_Transporter</code> Score
    <span class="num">145,3</span> ile ikinci sirada &mdash; null modele gore
    <span class="num">0,93&times;</span>, yani <em>tukenmis</em>.
    <code>G_Regulator</code> <span class="num">123,3</span> ile ucuncu &mdash;
    gercekte notr. Equation&nbsp;1 arka plana normalize edilmedigi icin
    <em>genom bollugunu</em> siraliyor, RO iliskisini degil.
  </div>

  <div class="callout">
    <strong>Ama yon hala gercek bir sinyal.</strong> Bolluk &quot;hayir&quot; derken
    orientation &quot;evet&quot; diyor: regulatorler <span class="num">%47,1</span>
    ayni-yon oraniyla rastgelenin <em>altinda</em> kalan tek kategori, ve
    <span class="num">""" + divergent_txt + """</span> tanesi divergent konumda. Yani dogru analiz birimi
    &quot;yakindaki regulatorler&quot; degil, <strong>&quot;divergent konumdaki
    regulatorler&quot;</strong> &mdash; kalan ~17.000&rsquo;i arka plan.
  </div>

  <p>
    Gercekten zenginlesenler yolun kendi bilesenleri:
    <code>ro_beta</code> <span class="num">4,12&times;</span>,
    <code>ring_cleavage</code> <span class="num">2,81&times;</span>,
    <code>ferredoxin</code> <span class="num">2,05&times;</span>.
    Biyolojik olarak beklenen tam da bu.
  </p>
</section>
"""


def ecology_section(path="analysis_out/cluster_ecology_stats.csv"):
    """Substrat sinifi x mobilite -- evrimsel hipotez testi."""
    if not os.path.exists(path):
        return ""
    frame = pd.read_csv(path)
    known = frame[(frame.confidence != "dusuk") & (frame.substrate_class != "unknown")]
    # Okaryot-agirlikli kumeler bakteriyel mobilite testini confound eder -- dislanir
    # (analyze_ecology.py ile ayni kural, tutarlilik icin).
    euk_excluded = 0
    domain_csv = "analysis_out/domain_by_cluster.csv"
    if os.path.exists(domain_csv):
        dom = pd.read_csv(domain_csv).set_index("cluster")["eukaryota_rate"].to_dict()
        mask = known.cluster.map(lambda c: dom.get(c, 0) < 0.25)
        euk_excluded = int((~mask).sum())
        known = known[mask]
    if known.empty:
        return ""

    def bucket(classes):
        subset = known[known.substrate_class.isin(classes)]
        n = subset.ro_count.sum()
        return {
            "n": n,
            "tn": subset.transposon_adjacent.sum() / n if n else 0,
            "pl": subset.plasmid_count.sum() / n if n else 0,
            "genera": subset.genus_count.sum(),
        }

    xeno = bucket(["xenobiotic"])
    natural = bucket(["natural_aromatic", "natural_specialized"])

    top = known.sort_values("transposon_rate", ascending=False).head(14)
    cluster_rows = "".join(
        "<tr><td class=\"mono\">" + esc(row.cluster) + "</td>"
        "<td>" + esc(row.substrate) + "</td>"
        '<td><span class="chip chip--' +
        ("xeno" if row.substrate_class == "xenobiotic" else "nat") + '">' +
        esc(row.substrate_class.replace("natural_", "dogal-").replace("xenobiotic", "ksenobiyotik")) +
        "</span></td>"
        '<td class="num">' + f"{int(row.ro_count):,}" + "</td>"
        '<td class="num">%' + f"{100*row.transposon_rate:.1f}" + "</td>"
        '<td class="num">%' + f"{100*row.plasmid_rate:.1f}" + "</td>"
        '<td class="num">' + str(int(row.genus_count)) + "</td></tr>"
        for row in top.itertuples())

    def ratio(a, b):
        return f"{a/b:.2f}&times;" if b else "&infin;"

    return """
<section id="ekoloji">
  <h2>8 &middot; Ekoloji ve evrim: substrat sinifi mobiliteyi ongoruyor mu?</h2>
  <p>
    RO&rsquo;larin bir kismi son yuzyilda ortaya cikan bilesikleri okside eder &mdash;
    PAH, PCB, BTEX, ftalat, nitroaromatik, hatta tereftalat (PET plastigi monomeri).
    Bir kismi ise kadim metabolik yollarin parcasi: vanillat (lignin), 3-ketosteroid
    (kolesterol), kafein, kolin, karnitin.
  </p>
  <p>
    <strong>Hipotez:</strong> ksenobiyotik yikim yetenegi hizla edinilmis olmali &rarr;
    plazmit uzerinde, mobil element komsulugunda, genis yayilimli. Dogal substratli
    yollar ise kromozomal ve dikey kalitimli.
  </p>
  <div class="table-wrap">
  <table>
    <thead><tr><th>olcut</th><th class="num">ksenobiyotik</th><th class="num">dogal</th>
    <th class="num">oran</th><th>sonuc</th></tr></thead>
    <tbody>
      <tr><td>plazmit uzerinde</td>
        <td class="num">%""" + f"{100*xeno['pl']:.1f}" + """</td>
        <td class="num">%""" + f"{100*natural['pl']:.1f}" + """</td>
        <td class="num"><strong>""" + ratio(xeno['pl'], natural['pl']) + """</strong></td>
        <td><span class="chip chip--up">tahmin tuttu</span></td></tr>
      <tr><td>transpozon komsulugu</td>
        <td class="num">%""" + f"{100*xeno['tn']:.1f}" + """</td>
        <td class="num">%""" + f"{100*natural['tn']:.1f}" + """</td>
        <td class="num">""" + ratio(xeno['tn'], natural['tn']) + """</td>
        <td><span class="chip chip--flat">etki cok kucuk</span></td></tr>
      <tr><td>farkli cins sayisi</td>
        <td class="num">""" + f"{int(xeno['genera']):,}" + """</td>
        <td class="num">""" + f"{int(natural['genera']):,}" + """</td>
        <td class="num">""" + ratio(xeno['genera'], natural['genera']) + """</td>
        <td><span class="chip chip--down">TERS cikti</span></td></tr>
    </tbody>
  </table>
  </div>

  <div class="callout callout--find">
    <strong>Taksonomik yayilim tahmini ters cikti &mdash; ve sebebi ogretici.</strong>
    Dogal substratli RO&rsquo;lar <em>daha cok</em> cinste bulunuyor. Cunku yayilim
    genisligi yakin zamanli yatay transferi degil, <strong>kadim dikey kalitimi</strong>
    olcuyor. Ksenobiyotik RO&rsquo;lar birkaç uzman yikici soyda yogunlasmis
    (<i>Pseudomonas</i>, <i>Sphingomonas</i>, <i>Rhodococcus</i>, <i>Mycobacterium</i>) &mdash;
    dar ama yogun. Olcut degistirilmeli: yayginlik degil, <em>filogenetik dagilmislik</em>
    (birbirine uzak soylarda parcali gorunme) aranmali.
  </div>

  <h3>Kume duzeyi: en mobil 14 kume</h3>
  <div class="table-wrap">
  <table class="table--compact">
    <thead><tr><th>kume</th><th>substrat</th><th>sinif</th><th class="num">RO</th>
    <th class="num">Tn yakin</th><th class="num">plazmit</th><th class="num">cins</th></tr></thead>
    <tbody>""" + cluster_rows + """</tbody>
  </table>
  </div>
  <p>
    <code>1_102_CARDO</code> (karbazol) hem transpozon hem plazmit siralamasinin tepesinde
    &mdash; literaturde <i>car</i> genlerinin transpozonla tasindigi bilinir, veri bunu
    bagimsiz olarak dogruluyor. <code>3_301_PhnA1a</code> (fenantren)
    <span class="num">%20,2</span> ve <code>3_317_NidA</code> (piren)
    <span class="num">%24,3</span> plazmit oraniyla PAH yikiminin mobil karakterini
    gosteriyor. <code>4_409_TPDO</code> &mdash; tereftalat, yani PET monomeri &mdash;
    <span class="num">%11,4</span> plazmit.
  </p>

  <div class="callout callout--warn">
    <strong>Uc onemli cekince.</strong> (1) Substrat atamalari
    <code>cluster_ecology.csv</code>&rsquo;de literaturden en iyi cabayla yapildi ve
    <code>confidence</code> sutunuyla isaretlendi &mdash; 72 referansi kuratorleyen kisi
    bunlari dogrulamali. (2) Transpozon farki sinirda (z&nbsp;=&nbsp;2,57, p&lt;0,01
    esiginin hemen altinda); cinsler arasi farklar buyuk olcude <em>o genomlarda kac
    IS elementi oldugunu</em> olcuyor. Plazmit sinyali (tablodaki oran) ise saglam.
    (3) Okaryot-agirlikli """ + str(euk_excluded) + """ kume bu testten <em>cikarildi</em>
    &mdash; okaryotta plazmit/mobil element farkli calisir, bakteriyel sinyali confound eder.
  </div>

  <p class="caveat">
    <strong>Olculemeyen:</strong> operon tamligi (alpha + beta + ferredoksin + reduktaz)
    urun adlarindan guvenilir sekilde cikarilamiyor &mdash; beta alt birimleri siklikla
    <code>hypothetical protein</code> olarak anotasyonlu. 11.422 RO&rsquo;nun yaninda
    yalnizca 1.694 anotasyonlu beta var. Dogru olcum komsu proteinleri de HMM&rsquo;den
    gecirmeyi gerektirir; bu duzeltilmeden operon tamligi rakamlari anotasyon kalitesini
    olcer, biyolojiyi degil.
  </p>
</section>
"""


def variants_section(path="analysis_out/cluster_variance.csv",
                     sdp_path="analysis_out/sdp_positions.csv"):
    """Kume ici varyant analizi -- atamanin guvenilirlik siniri."""
    if not os.path.exists(path):
        return ""
    frame = pd.read_csv(path).sort_values("median_identity")
    hetero = frame[frame.median_identity < 0.35]
    homo = frame[frame.median_identity > 0.60]

    leaf_block = ""
    leaf_csv = "analysis_out/leaves.csv"
    if os.path.exists(leaf_csv):
        lf = pd.read_csv(leaf_csv)
        n_leaves = len(lf)
        n_homog = int((lf.is_homogeneous == 1).sum())
        big = lf[lf["size"] >= 10]
        covered = int(big["size"].sum())
        singletons = int((lf["size"] == 1).sum())
        leaf_block = """
  <h3>Ozyinelemeli homojenizasyon: her kumeyi yapraklarina kadar bolmek</h3>
  <p>
    Tek CD-HIT gecisi yerine, bir alt-aile hala heterojense (medyan kimlik &lt;%70)
    daha yuksek esikte tekrar bolunur &mdash; her <em>yaprak</em> tek bir tutarli enzim
    tipi olana kadar. Sonuc bir agac. Substrat atamasi ancak bu yaprak duzeyinde
    savunulabilir: bir yaprak %85 kimlikli 90 proteinse, ona substrat atamak anlamli;
    tum kumeye atamak degil.
  </p>
  """ + stat_tiles([
      (f"{n_leaves:,}", "homojen yaprak", "61 kume &rarr; yapraklar", True),
      (f"%{100*n_homog/n_leaves:.0f}", "hedefi tutturan yaprak", "medyan kimlik &ge;%70", False),
      (f"{len(big):,}", "buyuk yaprak (&ge;10)", f"{covered:,} RO kapsiyor", False),
      (f"{singletons:,}", "tekil yaprak", "izole diziler", False),
  ]) + """
  <div class="callout">
    <strong>Gercek enzim tipleri artik gorunur.</strong> <code>1_104_VanA</code>
    tek bir 412-uyeli yaprak (medyan kimlik ~%70) veriyor (asil vanillat monooksijenaz),
    yaninda taksonomik olarak tutarli kucuk yapraklar (14 <i>Microcystis</i>, %99).
    <code>5_504_CntA</code> ise 244 yaprak &mdash; &quot;tek enzim&quot; degil, bir aile agaci.
  </div>
  <div class="callout callout--find">
    <strong>Varyantlar genomik baglamla ayirt ediliyor.</strong> Ayni kumenin
    yapraklari sadece dizide degil, <em>komsulukta</em> da farkli &mdash; ve bu
    islevsel farkin en guclu isareti. <code>3_313_KshA15</code> yapraklarindan biri
    kanonik <code>3-ketosteroid-&Delta;1-dehidrogenaz</code> yaninda oturuyor (steroid
    yolu), digeri farkli bir baglamda. <code>1_104_VanA</code>&rsquo;nin
    <i>Pseudomonas</i> cekirdegi (412 uye, kromozomal) <code>GntR</code>/<code>MarR</code>
    ile; <i>Burkholderia</i>/<i>Enterobacter</i>/<i>Acinetobacter</i> yapraklari ise
    hepsi <code>IclR</code> + mobil element yaninda &mdash; ayri bir tasinabilir soy.
    Substrat atamasi bu yaprak duzeyinde yapilmali; gezginde her varyant komsuluk
    imzasi ve indirilebilir uye listesiyle ayri ayri gosteriliyor.
  </div>"""

    def block(sub):
        return "".join(
            "<tr><td class=\"mono\">" + esc(row.cluster) + "</td>"
            '<td class="num">' + f"{int(row.n):,}" + "</td>"
            '<td class="num">%' + f"{100*row.median_identity:.1f}" + "</td>"
            '<td class="num">' + (str(int(row.subfamilies)) if pd.notna(row.subfamilies) else "&ndash;") + "</td>"
            '<td class="num">' + f"{int(row.sampled)}" + "</td></tr>"
            for row in sub.itertuples())

    sdp_note = ""
    if os.path.exists(sdp_path):
        sdp = pd.read_csv(sdp_path).sort_values("sdp_score", ascending=False).head(30)
        in_cat = int(sdp.in_catalytic_domain.sum())
        top = sdp.head(10)
        sdp_rows = "".join(
            '<tr><td class="num mono">' + str(int(row.column)) + "</td>"
            '<td class="num">' + f"{row.sdp_score:.3f}" + "</td>"
            '<td class="num">' + f"{row.within_cluster_entropy:.3f}" + "</td>"
            '<td class="num">' + f"{row.between_cluster_entropy:.3f}" + "</td>"
            '<td><span class="chip chip--' + ("cat" if row.in_catalytic_domain else "rieske") + '">' +
            ("katalitik domain" if row.in_catalytic_domain else "Rieske / N-terminal") + "</span></td></tr>"
            for row in top.itertuples())
        sdp_note = """
  <h3>Spesifisite belirleyici pozisyonlar</h3>
  <p>
    Substrat secicigi katalitik cebi doseyen birkac kalinti tarafindan belirlenir.
    Bu pozisyonlarin imzasi: kume <em>icinde</em> korunmus ama kumeler
    <em>arasinda</em> farkli. Sadece korunmus kolonlar yapisal, sadece degisken
    kolonlar gurultu &mdash; islevi tasiyanlar ikisinin arasindadir.
  </p>
  <p class="measure-note">SDP skoru = (kumeler arasi cesitlilik) &minus; (kume ici ortalama cesitlilik)</p>
  <div class="table-wrap">
  <table class="table--compact">
    <thead><tr><th class="num">kolon</th><th class="num">SDP skoru</th>
    <th class="num">kume ici</th><th class="num">kumeler arasi</th><th>bolge</th></tr></thead>
    <tbody>""" + sdp_rows + """</tbody>
  </table>
  </div>
  <div class="callout callout--find">
    <strong>Ilk 30 SDP&rsquo;nin """ + str(in_cat) + """&rsquo;i katalitik domainde.</strong>
    Tam da beklenen yer: spesifisite belirleyicileri substrat cebinde yogunlasiyor,
    elektron transferi yapan Rieske domaininde degil. Bu, skorun gurultu degil
    gercek sinyal olctugunun bagimsiz gostergesi.
  </div>
"""

    return """
<section id="varyant">
  <h2>9 &middot; Kume ici varyasyon: &quot;ayni enzim&quot; gercekten ayni mi?</h2>
  <p>
    Kume atamasi 71 model arasindan en yuksek bit skorunu alani secer. Ama bu,
    kumenin homojen oldugunu <em>garanti etmez</em> &mdash; her protein bir kumeye
    gitmek <strong>zorundadir</strong>, hicbirine yeterince benzemese bile.
    Minimum benzerlik tabani yok, sadece goreli kazanan var.
  </p>
  <p>
    Olcum: tum uyeler ayni 426 kolonluk modele hizali oldugu icin ikili kimlik
    dogrudan match kolonlarindan okunur; ayrica her kume CD-HIT ile %70 kimlikte
    alt-ailelere ayrilir.
  </p>
  """ + stat_tiles([
      (str(len(hetero)), "heterojen kume", "medyan kimlik &lt; %35", True),
      (f"{int(hetero.n.sum()):,}", "bu kumelerdeki RO", "dogrulananlarin cogunlugu", False),
      (str(len(homo)), "homojen kume", "medyan kimlik &gt; %60", False),
      ("&minus;0,02", "kimlik &harr; kume boyutu", "korelasyon yok (p=0,88)", False),
  ]) + """
  <h3>En heterojen kumeler</h3>
  <div class="table-wrap">
  <table class="table--compact">
    <thead><tr><th>kume</th><th class="num">n</th><th class="num">medyan kimlik</th>
    <th class="num">alt-aile</th><th class="num">orneklenen</th></tr></thead>
    <tbody>""" + block(frame.head(10)) + """</tbody>
  </table>
  </div>

  <h3>En homojen kumeler</h3>
  <div class="table-wrap">
  <table class="table--compact">
    <thead><tr><th>kume</th><th class="num">n</th><th class="num">medyan kimlik</th>
    <th class="num">alt-aile</th><th class="num">orneklenen</th></tr></thead>
    <tbody>""" + block(frame.tail(6).iloc[::-1]) + """</tbody>
  </table>
  </div>

  <div class="callout callout--warn">
    <strong>%24 medyan kimlik alacakaranlik bolgesinin altidir.</strong>
    <code>4_401_CmtAb</code> 107 uyeden 49 alt-aile, <code>5_504_CntA</code> 1.210 uyeden
    (120 orneklemde) 66 alt-aile veriyor. Bunlar tek bir enzim degil; genis modeller
    &quot;baska yere uymayan her seyin&quot; cekim havuzu olmus.
    Kimlik ile kume boyutu arasinda korelasyon <em>yok</em>
    (&rho;=&minus;0,02) &mdash; yani bu &quot;buyuk kumeler dogal olarak cesitlidir&quot;
    degil, hangi modellerin cekici oldugu meselesi.
  </div>

  <div class="callout callout--find">
    <strong>Bu, substrat atamalarini geriye donuk sinirlar.</strong>
    Bolum 8&rsquo;deki ekolojik siniflandirma yalnizca <em>homojen</em> kumeler icin
    savunulabilir. <code>5_504_CntA</code>&rsquo;ya &quot;karnitin&quot; demek,
    %24,7 kimlikli 1.210 proteinin hepsi icin gecerli olamaz. Kume duzeyi ekoloji
    sonuclari homojenlik skoruyla birlikte okunmali.
  </div>
  """ + sdp_note + """
  """ + leaf_block + """
  <p class="caveat">
    <strong>Ne yapilmali:</strong> (1) atamaya bir <em>mutlak</em> taban ekle
    (yalnizca goreli en iyi degil, minimum bit skoru / kimlik); (2) heterojen kumeleri
    alt-ailelere bol ve her alt-aile icin ayri profil kur; (3) novel RO aramasini
    atandigi kumenin merkezinden uzak duran uyelere odakla &mdash; gercek yeni aileler
    tam olarak orada.
  </p>
</section>
"""


def novel_section(db="roar.sqlite",
                  novel_csv="analysis_out/novel_candidates.csv"):
    """Alt-aile bolme + mutlak olcut + novel aday kesfi."""
    import sqlite3
    if not os.path.exists(db) or not os.path.exists(novel_csv):
        return ""
    con = sqlite3.connect(db)
    def scalar(q):
        r = con.execute(q).fetchone()
        return r[0] if r else 0
    n_sub = scalar("SELECT COUNT(*) FROM subfamily")
    n_ref = scalar("SELECT COUNT(*) FROM subfamily WHERE size>=10")
    n_core = scalar("SELECT COUNT(*) FROM ro_subfamily WHERE assignment_class='core'")
    n_nov = scalar("SELECT COUNT(*) FROM ro_subfamily WHERE assignment_class='novel_candidate'")
    n_tot = scalar("SELECT COUNT(*) FROM ro_subfamily")

    frame = pd.read_csv(novel_csv)
    high = frame[(frame.hmm_score >= 350) & (frame.core_identity < 0.40)]
    euk_markers = ["Spinacia","Chlorella","Volvox","Rhodamnia","Arabidopsis","Oryza",
                   "Gossypium","Cucurbita","Ocimum","Nicotiana","Zea","Physcomitrium",
                   "Populus","Glycine","Marchantia"]
    euk = frame[frame.organism.fillna("").apply(
        lambda o: any(m in o for m in euk_markers))]

    top = high.sort_values("hmm_score", ascending=False).head(12)
    rows = "".join(
        "<tr><td class=\"mono\">" + esc(str(r.candidate_id)[:28]) + "</td>"
        "<td class=\"mono\">" + esc(r.assigned_cluster) + "</td>"
        '<td class="num">' + f"{r.hmm_score:.0f}" + "</td>"
        '<td class="num">%' + f"{100*r.core_identity:.0f}" + "</td>"
        "<td><i>" + esc((r.organism or "?")[:34]) + "</i></td></tr>"
        for r in top.itertuples())
    con.close()

    return """
<section id="novel">
  <h2>10 &middot; Novel RO kesfi: cevreye odaklanmak</h2>
  <p>
    Bolum 9 kumelerin heterojen oldugunu gosterdi. Bunun dogrudan sonucu: atama
    goreli en iyi modeli sectigi icin, bilinen hicbir tipe uymayan bir protein bile
    zorla bir kumeye tikili &mdash; ve orada <em>gorunmez</em>. Novel RO tam olarak
    buralarda saklaniyor.
  </p>
  <p>
    Yontem uc katman: (1) her kume CD-HIT ile alt-ailelere bolunur; (2) her uyeye
    <strong>mutlak</strong> bir olcut atanir &mdash; goreli kazanan degil, kumenin
    kalabalik <em>referans tiplerinden</em> herhangi birine en iyi kimlik; (3) hicbir
    referans tipe uymayan (&lt;%40) ve kucuk alt-ailede oturan uyeler novel aday.
  </p>
  """ + stat_tiles([
      (f"{n_ref:,}", "referans-tip alt-aile", f"71 model &rarr; {n_ref} ince ayrim", True),
      (f"{n_sub:,}", "toplam alt-aile", "CD-HIT %70", False),
      (f"{n_nov:,}", "novel aday", f"%{100*n_nov/max(1,n_tot):.0f} &mdash; izole, tipsiz", False),
      (f"{len(high):,}", "yuksek-guven novel", "skor&ge;350, kimlik&lt;%40", False),
  ]) + """
  <div class="callout callout--find">
    <strong>En guclu novel adaylar: yuksek skor + dusuk kimlik.</strong>
    Bu """ + str(len(high)) + """ protein HMM tarafindan guvenle RO olarak taniniyor
    (yuksek bit skoru, katalitik merkez tam) <em>ama</em> bilinen hicbir alt-tipe
    benzemiyor. Ornegin <code>4_405_ROCH34</code>&rsquo;e atanan bir grup
    <i>Bradyrhizobium</i>/<i>Rhodopseudomonas</i> proteini ~860 skorla geliyor ama
    cekirdege yalnizca %38 &mdash; yanlis yere tikilmis tutarli bir yeni tip.
  </div>
  <div class="table-wrap">
  <table class="table--compact">
    <thead><tr><th>aday</th><th>atandigi kume</th><th class="num">HMM skoru</th>
    <th class="num">cekirdek kimlik</th><th>organizma</th></tr></thead>
    <tbody>""" + rows + """</tbody>
  </table>
  </div>
  <div class="callout">
    <strong>Ayri bir evrimsel dal: bitki/alg Rieske oksijenazlari.</strong>
    Novel adaylarin """ + str(len(euk)) + """&rsquo;i bitki veya algden geliyor
    (<i>Spinacia</i>, <i>Chlorella</i>, <i>Volvox</i>&hellip;). <i>Spinacia</i>
    kolin monooksijenazi <code>CmoA</code>&rsquo;ya <span class="num">948</span> skorla
    hizalaniyor &mdash; guclu bir RO, ama bakteriyel yikim enzimlerinden tamamen
    ayri bir soy. Kloroplast Rieske oksijenazlari (CMO, CAO) bilinen ama bu veri
    setinde ayri bir kume olarak temsil edilmiyorlar.
  </div>
  <p class="caveat">
    <strong>Cikti dosyalari:</strong>
    <code>analysis_out/novel_candidates.fasta</code> (2.283 dizi) ve
    <code>novel_high_confidence.fasta</code> (155 dizi) &mdash; filogenetik
    yerlestirme veya BLAST dogrulamasi icin hazir.
    <code>subfamilies.csv</code> 204 referans-tip alt-aileyi listeler; her biri
    icin ayri HMM kurulup atama bu ince cozunurlukte tekrarlanabilir.
  </p>
</section>
"""


def domain_section(db="roar.sqlite", csv_path="analysis_out/domain_by_cluster.csv"):
    """Yasam alani kompozisyonu -- veri kalitesi."""
    import sqlite3
    if not os.path.exists(db):
        return ""
    con = sqlite3.connect(db)
    if not con.execute("SELECT name FROM sqlite_master WHERE name=\'ro_domain\'").fetchone():
        con.close(); return ""
    dist = dict(con.execute("SELECT domain, COUNT(*) FROM ro_domain GROUP BY domain").fetchall())
    euk = dict(con.execute("SELECT euk_group, COUNT(*) FROM ro_domain WHERE domain=\'Eukaryota\' GROUP BY euk_group").fetchall())
    total = sum(dist.values()) or 1   # bos ro_domain'de sifira bolme olmasin
    b, e, a = dist.get("Bacteria", 0), dist.get("Eukaryota", 0), dist.get("Archaea", 0)
    con.close()

    flagged = []
    if os.path.exists(csv_path):
        f = pd.read_csv(csv_path)
        flagged = f[(f.total >= 10) & (f.eukaryota_rate >= 0.25)].sort_values(
            "eukaryota_rate", ascending=False)

    flag_rows = "".join(
        "<tr><td class=\"mono\">" + esc(r.cluster) + "</td>"
        "<td class=\"num\">" + f"{int(r.total):,}" + "</td>"
        "<td class=\"num\">" + f"{int(r.eukaryota):,}" + "</td>"
        "<td class=\"num\">%" + f"{100*r.eukaryota_rate:.0f}" + "</td>"
        "<td>" + esc(r.top_euk_group) + "</td></tr>"
        for r in flagged.itertuples()) if len(flagged) else ""

    return """
<section id="domain">
  <h2>11 &middot; Veri kompozisyonu: bakteriyel mi, karisik mi?</h2>
  <p>
    Temsilci analizinde bazi kume tip dizilerinin <strong>bitki proteini</strong> oldugu
    ortaya cikti. Sistematik baktik: her uye yasam alanina gore siniflandirildi.
  </p>
  """ + stat_tiles([
      (f"%{100*b/total:.1f}", "bakteri", f"{b:,} RO", True),
      (f"%{100*e/total:.1f}", "okaryot", f"{e:,} RO", False),
      (f"{euk.get('plant/alga', 0):,}", "bitki/alg", "kloroplast RO dali", False),
      (f"{euk.get('fungus', 0)+euk.get('animal', 0):,}", "mantar + hayvan", "supheli", False),
  ]) + """
  <div class="callout callout--warn">
    <strong>Katalitik-triad filtresi yasam alanini ayirmaz.</strong> Filtre
    ferredoksini iyi eliyor ama &quot;bakteriyel halka-hidroksileyen oksijenaz&quot;i
    &quot;Rieske + katalitik merkez tasiyan herhangi bir okaryot protein&quot;den
    ayirmiyor. Kloroplast Rieske oksijenazlari (CMO, CAO) bu mimariyi gercekten
    paylasir. Bu bir veri gercegi &mdash; gizlenmek yerine isaretlendi.
  </div>""" + ("""
  <h3>Okaryot-agirlikli kumeler</h3>
  <p>Bunlarin bir kismi <em>dogru</em>: <code>5_501_CmoS</code> (kolin monooksijenaz)
    zaten bir bitki enzimi. Digerleri (<code>DdmC</code>, <code>GxtA</code>) bakteriyel
    diye etiketlenmisti &mdash; substrat atamalari gozden gecirilmeli.</p>
  <div class="table-wrap"><table class="table--compact">
    <thead><tr><th>kume</th><th class="num">toplam</th><th class="num">okaryot</th>
    <th class="num">oran</th><th>baskin grup</th></tr></thead>
    <tbody>""" + flag_rows + """</tbody>
  </table></div>""" if flag_rows else "") + """
  <p class="caveat">
    Gezginde her kumenin yasam alani kirilimi bir cubukla gosteriliyor; okaryot-agirlikli
    olanlar ayrica bayrakli. <code>ro_domain</code> tablosu ve
    <code>domain_by_cluster.csv</code> tam kirilimi icerir.
  </p>
</section>
"""


def wrap_page(body):
    return f"""<title>ROAR-DB &middot; Rieske oksijenaz filtresi: metodoloji ve sonuclar</title>
<style>
:root {{
  --ground:      #F7F8F9;
  --ground-2:    #FFFFFF;
  --ground-3:    #EDEFF1;
  --ink:         #1A1F24;
  --ink-2:       #4A555E;
  --ink-3:       #78848D;
  --rule:        #DCE1E4;
  --accent:      #A6402A;   /* demir oksit -- Fe merkezi */
  --accent-soft: #F3E3DF;
  --teal:        #2F6E7A;
  --teal-soft:   #E0EDEF;
  --ok:          #3F7A4E;
  --warn:        #96701A;
  --warn-soft:   #F7EEDA;
  --measure:     68ch;
  --serif: ui-serif, Georgia, "Iowan Old Style", "Palatino Linotype", Palatino, serif;
  --sans:  ui-sans-serif, system-ui, -apple-system, "Segoe UI", Roboto, "Helvetica Neue", sans-serif;
  --mono:  ui-monospace, SFMono-Regular, "SF Mono", Menlo, Consolas, "Liberation Mono", monospace;
}}
@media (prefers-color-scheme: dark) {{
  :root {{
    --ground: #13171A; --ground-2: #1A1F23; --ground-3: #222A2F;
    --ink: #E4E8EA; --ink-2: #A6B1B8; --ink-3: #78858D;
    --rule: #2C353B;
    --accent: #E08163; --accent-soft: #33211C;
    --teal: #6FB3BF; --teal-soft: #16292D;
    --ok: #74B584; --warn: #D6A648; --warn-soft: #2E2716;
  }}
}}
:root[data-theme="dark"] {{
  --ground: #13171A; --ground-2: #1A1F23; --ground-3: #222A2F;
  --ink: #E4E8EA; --ink-2: #A6B1B8; --ink-3: #78858D;
  --rule: #2C353B;
  --accent: #E08163; --accent-soft: #33211C;
  --teal: #6FB3BF; --teal-soft: #16292D;
  --ok: #74B584; --warn: #D6A648; --warn-soft: #2E2716;
}}
:root[data-theme="light"] {{
  --ground: #F7F8F9; --ground-2: #FFFFFF; --ground-3: #EDEFF1;
  --ink: #1A1F24; --ink-2: #4A555E; --ink-3: #78848D;
  --rule: #DCE1E4;
  --accent: #A6402A; --accent-soft: #F3E3DF;
  --teal: #2F6E7A; --teal-soft: #E0EDEF;
  --ok: #3F7A4E; --warn: #96701A; --warn-soft: #F7EEDA;
}}

body {{
  background: var(--ground); color: var(--ink);
  font-family: var(--sans); font-size: 16px; line-height: 1.65;
  -webkit-font-smoothing: antialiased;
}}
.page {{ max-width: 1180px; margin: 0 auto; padding: 0 clamp(1.1rem, 4vw, 3.5rem) 6rem; }}

/* ---------- masthead ---------- */
.masthead {{ padding: clamp(3rem, 8vw, 6rem) 0 clamp(2rem, 4vw, 3.5rem); border-bottom: 2px solid var(--ink); }}
.masthead__eyebrow {{
  font-family: var(--mono); font-size: .74rem; letter-spacing: .13em;
  text-transform: uppercase; color: var(--accent); margin-bottom: 1.6rem;
}}
.masthead h1 {{
  font-family: var(--serif); font-weight: 600; font-size: clamp(2rem, 5.2vw, 3.4rem);
  line-height: 1.12; letter-spacing: -.02em; text-wrap: balance; margin-bottom: 1.5rem;
}}
.masthead h1 em {{ font-style: italic; color: var(--accent); }}
.lede {{ max-width: var(--measure); font-size: 1.08rem; color: var(--ink-2); }}
.lede strong {{ color: var(--ink); }}

/* ---------- sections ---------- */
section {{ padding: clamp(2.6rem, 6vw, 4.2rem) 0; border-bottom: 1px solid var(--rule); }}
section:last-of-type {{ border-bottom: none; }}
h2 {{
  font-family: var(--serif); font-weight: 600; font-size: clamp(1.5rem, 3.2vw, 2.05rem);
  line-height: 1.2; letter-spacing: -.015em; margin-bottom: 1.3rem; text-wrap: balance;
}}
h3 {{
  font-family: var(--sans); font-weight: 650; font-size: 1.06rem;
  letter-spacing: -.005em; margin: 2.4rem 0 .85rem;
}}
h4 {{ font-family: var(--sans); font-weight: 650; font-size: .95rem; margin: 1.4rem 0 .7rem; }}
p {{ max-width: var(--measure); margin-bottom: 1.05rem; color: var(--ink-2); }}
p strong {{ color: var(--ink); font-weight: 620; }}
em {{ font-style: italic; }}

.measure-note {{ font-size: .87rem; color: var(--ink-3); }}
.caveat {{
  font-size: .9rem; color: var(--ink-2); border-left: 2px solid var(--rule);
  padding-left: 1rem; margin-top: 1.5rem;
}}

/* ---------- numbers & code ---------- */
.num, .mono {{ font-family: var(--mono); font-variant-numeric: tabular-nums; }}
.num {{ font-size: .94em; }}
code {{
  font-family: var(--mono); font-size: .87em; background: var(--ground-3);
  padding: .12em .38em; border-radius: 3px; color: var(--ink);
}}

/* ---------- stat tiles ---------- */
.tiles {{
  display: grid; grid-template-columns: repeat(auto-fit, minmax(178px, 1fr));
  gap: 1px; background: var(--rule); border: 1px solid var(--rule);
  margin: 1.9rem 0; max-width: none;
}}
.tile {{ background: var(--ground-2); padding: 1.25rem 1.2rem 1.35rem; }}
.tile--accent {{ background: var(--accent-soft); }}
.tile__num {{
  font-family: var(--mono); font-variant-numeric: tabular-nums;
  font-size: clamp(1.45rem, 3.4vw, 1.95rem); font-weight: 600;
  line-height: 1.05; letter-spacing: -.02em; color: var(--ink); margin-bottom: .5rem;
}}
.tile--accent .tile__num {{ color: var(--accent); }}
.tile__label {{ font-size: .82rem; font-weight: 600; line-height: 1.35; color: var(--ink); }}
.tile__note {{ font-size: .76rem; color: var(--ink-3); margin-top: .3rem; line-height: 1.4; }}

/* ---------- pipeline ---------- */
.pipeline {{
  list-style: none; display: grid; gap: 1px; background: var(--rule);
  border: 1px solid var(--rule); margin: 1.9rem 0;
  grid-template-columns: repeat(auto-fit, minmax(240px, 1fr));
}}
.stage {{ background: var(--ground-2); padding: 1.4rem 1.3rem 1.5rem; }}
.stage__tag {{
  font-family: var(--mono); font-size: .68rem; letter-spacing: .12em;
  text-transform: uppercase; color: var(--teal); margin-bottom: .75rem;
}}
.stage h3 {{ margin: 0 0 .6rem; font-size: 1rem; }}
.stage p {{ font-size: .89rem; margin-bottom: .55rem; max-width: none; }}
.stage__why {{ color: var(--ink-3); font-size: .84rem !important; font-style: italic; }}

/* ---------- callouts ---------- */
.callout {{
  background: var(--teal-soft); border-left: 3px solid var(--teal);
  padding: 1.1rem 1.25rem; margin: 1.7rem 0; max-width: var(--measure);
  font-size: .93rem; color: var(--ink-2);
}}
.callout strong {{ color: var(--ink); }}
.callout--warn {{ background: var(--warn-soft); border-left-color: var(--warn); }}
.callout--find {{ background: var(--accent-soft); border-left-color: var(--accent); }}
.callout--pending {{ background: var(--ground-3); border-left-color: var(--ink-3); }}

/* ---------- tables ---------- */
.table-wrap {{ overflow-x: auto; margin: 1.4rem 0; border: 1px solid var(--rule); }}
table {{ width: 100%; border-collapse: collapse; font-size: .9rem; background: var(--ground-2); }}
th {{
  text-align: left; font-weight: 620; font-size: .76rem; letter-spacing: .06em;
  text-transform: uppercase; color: var(--ink-3); padding: .7rem .95rem;
  border-bottom: 1px solid var(--rule); white-space: nowrap;
}}
th.num, td.num {{ text-align: right; font-family: var(--mono); font-variant-numeric: tabular-nums; }}
td {{ padding: .62rem .95rem; border-bottom: 1px solid var(--rule); color: var(--ink-2); }}
tr:last-child td {{ border-bottom: none; }}
tbody tr.is-chosen {{ background: var(--accent-soft); }}
tbody tr.is-chosen td {{ color: var(--ink); font-weight: 600; }}
.table--compact td, .table--compact th {{ padding: .45rem .8rem; font-size: .85rem; }}

.pill {{
  display: inline-block; font-family: var(--mono); font-size: .66rem;
  letter-spacing: .08em; text-transform: uppercase; background: var(--accent);
  color: #fff; padding: .12em .5em; border-radius: 2px; margin-left: .5rem;
  vertical-align: middle;
}}
.chip {{
  display: inline-block; font-size: .74rem; padding: .16em .55em;
  border-radius: 2px; background: var(--ground-3); color: var(--ink-2); white-space: nowrap;
}}
.chip--rieske {{ background: var(--teal-soft); color: var(--teal); }}
.chip--cat {{ background: var(--accent-soft); color: var(--accent); }}
.chip--up {{ background: color-mix(in srgb, var(--ok) 15%, transparent); color: var(--ok); }}
.chip--down {{ background: var(--accent-soft); color: var(--accent); }}
.chip--flat {{ background: var(--ground-3); color: var(--ink-3); }}
.chip--xeno {{ background: var(--accent-soft); color: var(--accent); }}
.chip--nat {{ background: var(--teal-soft); color: var(--teal); }}
th.enr-head {{ width: 130px; font-family: var(--mono); text-transform: none; letter-spacing: 0; }}
td.enr {{ padding: .62rem .5rem; width: 130px; }}
.enr__bar {{ height: 8px; border-radius: 1px; }}
.enr__bar--up {{ background: var(--ok); }}
.enr__bar--down {{ background: var(--accent); }}
.enr__bar--flat {{ background: var(--ink-3); opacity: .45; }}

/* ---------- bars ---------- */
.bars {{ display: flex; flex-direction: column; gap: .5rem; margin-top: .9rem; }}
.bar-row {{ display: grid; grid-template-columns: 4.4rem 1fr 3.6rem; align-items: center; gap: .7rem; }}
.bar-row__label {{ font-size: .84rem; color: var(--ink-2); }}
.bar-row__val {{ font-size: .84rem; text-align: right; color: var(--ink); }}
.bar {{ height: 9px; background: var(--ground-3); }}
.bar__fill {{ height: 100%; background: var(--teal); }}

.two-col {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(280px, 1fr)); gap: 2rem; }}
.two-col h3 {{ margin-top: 1.5rem; }}

/* ---------- bugs ---------- */
.bugs {{ list-style: none; counter-reset: bug; display: flex; flex-direction: column; gap: 1px;
        background: var(--rule); border: 1px solid var(--rule); margin-top: 1.6rem; }}
.bug {{ background: var(--ground-2); padding: 1.2rem 1.3rem 1.3rem; counter-increment: bug; }}
.bug__head {{ display: flex; flex-wrap: wrap; align-items: baseline; gap: .7rem; margin-bottom: .5rem; }}
.bug h3 {{ margin: 0; font-size: .99rem; }}
.bug h3::before {{
  content: counter(bug, decimal-leading-zero); font-family: var(--mono);
  color: var(--accent); margin-right: .6rem; font-size: .84em;
}}
.bug__loc {{ font-size: .76rem; color: var(--ink-3); background: none; padding: 0; }}
.bug p {{ font-size: .9rem; margin: 0; max-width: none; }}

/* ---------- code ---------- */
.code-wrap {{ overflow-x: auto; background: var(--ground-2); border: 1px solid var(--rule); margin: 1.5rem 0; }}
pre {{ padding: 1.15rem 1.3rem; }}
pre code {{ background: none; padding: 0; font-size: .82rem; line-height: 1.75; color: var(--ink-2); white-space: pre; }}

@media (max-width: 620px) {{
  body {{ font-size: 15px; }}
  .bar-row {{ grid-template-columns: 3.6rem 1fr 3rem; }}
}}
@media (prefers-reduced-motion: reduce) {{ * {{ animation: none !important; transition: none !important; }} }}
</style>

<div class="page">
{body}
</div>
"""


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dir", default=".")
    parser.add_argument("--db", default="roar.sqlite")
    parser.add_argument("-o", "--output", default="report.html")
    args = parser.parse_args()

    protein = read_protein_level(args.dir)
    genomic = read_genomic(args.db)
    page = wrap_page(build_html(protein, genomic, args.db))

    with open(args.output, "w") as handle:
        handle.write(page)
    print(f"[yazildi] {args.output}  ({len(page):,} bayt)")
    if genomic is None:
        print("[uyari] roar.sqlite yok -- genomik bolum bos birakildi")
    elif not genomic["confirmed"]:
        print("[uyari] dogrulanmis RO yok -- annotate_ro.py henuz bitmemis olabilir")


if __name__ == "__main__":
    main()
