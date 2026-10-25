/* Giris sayfasini parcali JSON'dan ciz.
 *
 * Bu dosya templates/ro.html'in karsiligidir. Komsuluk SVG'si neighsvg.js'te
 * ve Python asliyla BIREBIR karsilastirilarak sinaniyor (verify_entry_render.py).
 * Operon-regulator SVG'si sunucuda hesaplanip sekil olarak geliyor; burada
 * yalnizca iki yer tutucu (__BP__, __PROD__) dolduruluyor, yani o geometri
 * JS'e port EDILMEDI.
 */
(function () {
  'use strict';

  var root = document.getElementById('entry-root');
  var BASE = root.getAttribute('data-base') || '';
  var esc = window.ROARneigh.esc;

  function slug(id) { return id.replace(/:/g, '_'); }
  function unslug(s) {
    var m = /^(.*)_(\d+)-(\d+)_(-?[01])$/.exec(s);
    return m ? m[1] + ':' + m[2] + '-' + m[3] + ':' + m[4] : s.replace(/_/g, ':');
  }
  function leafSlug(l) { return String(l).replace('#', '-'); }
  function fmt(n) {
    if (n === null || n === undefined || n === '') return '–';
    return Number(n).toLocaleString('en-US');
  }
  function num(v, d) {
    /* Python "%.<d>f" ile ayni yuvarlama -- bkz neighsvg.js fixed(). */
    return (v === null || v === undefined) ? '–' : window.ROARneigh.fixed(Number(v), d);
  }
  function expo(v, d) {
    /* Python "%.1e" ile ayni: US EN AZ IKI HANE olur.
       JS toExponential(1) "0.0e+0" verir, Python "0.0e+00". Olculdu:
       8 giriste 1'i yalnizca bu yuzden ayrisiyordu. */
    var s = Number(v).toExponential(d);
    return s.replace(/e([+-])(\d)$/, 'e$10$2');
  }
  function href(p) { return BASE + p; }

  var md5 = window.ROARmd5;

  function unpack(arr, fields) {
    if (!arr) return null;
    var o = {};
    for (var i = 0; i < fields.length; i++) o[fields[i]] = arr[i];
    return o;
  }

  function getJSON(path) {
    return fetch(href(path)).then(function (r) {
      if (!r.ok) throw new Error(path + ' ' + r.status);
      return r.json();
    });
  }

  function fail(msg) {
    root.innerHTML = '<div class="caveat"><b>Entry not found.</b> ' + esc(msg) +
      ' <a href="' + href('/search.html') + '">Search for an entry</a>.</div>';
  }

  var params = new URLSearchParams(location.search);
  var raw = params.get('id') || '';
  if (!raw) { fail('No entry identifier in the address.'); return; }
  var cid = raw.indexOf(':') >= 0 ? raw : unslug(raw);

  var META = null;
  getJSON('/data/meta.json').then(function (meta) {
    META = meta;
    var sh = md5(cid).slice(0, meta.shard_hex);
    return getJSON('/data/entry/' + sh + '.json').then(function (shard) {
      var rec = shard.e[cid];
      if (!rec) throw new Error('not in shard');
      return { shard: shard, rec: rec, sh: sh };
    });
  }).then(function (ctx) {
    return render(ctx.rec, ctx.shard, ctx.sh);
  }).catch(function (e) {
    fail('Identifier ' + esc(cid) + ' is not in this build (' + esc(e.message) + ').');
  });

  function render(rec, shard, sh) {
    var M = META;
    var ro = unpack(rec.ro, M.ro_fields);
    var ev = unpack(rec.ev, M.ev_fields);
    var etc = unpack(rec.etc, M.etc_fields);
    var reg = unpack(rec.reg, M.reg_fields);
    var op = unpack(rec.op, M.operon_fields);
    var str = shard.s || [];
    function S(i) { return (i === null || i === undefined) ? null : str[i]; }

    var nb = (rec.nb || []).map(function (row) {
      var o = unpack(row, M.nb_fields);
      M.nb_interned.forEach(function (f) { o[f] = S(o[f]); });
      return o;
    });
    var opg = (rec.opg || []).map(function (row) {
      var o = unpack(row, M.opgene_fields);
      M.opgene_interned.forEach(function (f) { o[f] = S(o[f]); });
      return o;
    });

    var gene = (M.cluster_gene || {})[ro.ro_cluster] || ro.ro_cluster;
    var h = [];

    h.push('<nav class="crumbs"><a href="' + href('/index.html') + '">Home</a> › ' +
      '<a href="' + href('/cluster/' + ro.ro_cluster + '.html') + '">' + esc(gene) + '</a>' +
      (ro.leaf_id ? ' › <a href="' + href('/leaf/' + leafSlug(ro.leaf_id) + '.html') + '">variant ' +
        esc(String(ro.leaf_id).split('#').pop()) + '</a>' : '') + ' › entry</nav>');
    h.push('<h1>' + esc(ro.product) + '</h1>');
    h.push('<p class="note mono">' + esc(ro.candidate_id) + '</p>');
    if (!ro.is_confirmed) {
      h.push('<div class="caveat"><b>Not confirmed.</b> This candidate failed the coverage or ' +
        'catalytic-motif test and is shown for completeness only.</div>');
    }

    // --- kanit duzeyi
    h.push('<div class="caveat"><b>What is known about this protein.</b> ');
    if (ev && ev.tier) {
      var tc = (M.tier_colors || {})[ev.tier] || '#888';
      h.push('<span class="chip chip--data" style="--chip:' + tc + '">' +
        esc(ev.tier.replace(/_/g, ' ')) + '</span> ' +
        esc((M.tier_meaning || {})[ev.tier] || '') + '<br>Nearest curated reference: ' +
        '<a href="' + href('/cluster/' + ev.nearest_ref + '.html') + '">' + esc(ev.nearest_ref) +
        '</a> at ' + num(ev.ref_identity, 1) + ' % identity over ' +
        num(ev.ref_qcov, 0) + ' % of the sequence.');
      if (ev.same_as_hmm === 0) {
        h.push(' The highest-scoring profile belongs to a different type, <a href="' +
          href('/cluster/' + ro.ro_cluster + '.html') + '">' + esc(ro.ro_cluster) +
          '</a>; where profile and nearest sequence disagree, neither should be read as a ' +
          'functional assignment.');
      }
      var sub = (M.cluster_substrate || {})[ro.ro_cluster];
      if (sub) {
        h.push('<br>The characterised enzyme of the assigned type acts on <b>' + esc(sub) +
          '</b>. Whether this protein shares that substrate follows from the evidence level ' +
          'above, not from the type name.');
      }
    } else {
      h.push('This protein carries an intact Rieske centre and catalytic triad. Its assigned ' +
        'type is the profile it scores highest against, which names its closest characterised ' +
        'relative rather than its function.');
    }
    h.push('</div>');

    // --- kunye
    h.push('<div class="cols"><dl class="kv">');
    h.push('<dt>Organism</dt><dd><i>' + esc(ro.organism) + '</i>' +
      (ro.is_plasmid ? '<span class="chip chip--plasmid">plasmid</span>' : '') +
      (ro.domain ? '<span class="chip">' + esc(ro.domain) +
        (ro.euk_group ? ' / ' + esc(ro.euk_group) : '') + '</span>' : '') + '</dd>');
    h.push('<dt>Lineage</dt><dd class="small">' + esc(ro.taxonomy) + '</dd>');
    h.push('<dt>Genome record</dt><dd><a href="https://www.ncbi.nlm.nih.gov/nuccore/' +
      encodeURIComponent(ro.nucleotide_id) + '" rel="noopener">' + esc(ro.nucleotide_id) +
      '</a> <span class="note">' + esc(ro.description) + ' · ' + fmt(ro.replicon_length) +
      ' bp · ' + fmt(ro.cds_count) + ' coding sequences' +
      (ro.is_circular ? ' · circular' : '') + '</span></dd>');
    h.push('<dt>Protein</dt><dd>' + (ro.protein_id
      ? '<a href="https://www.ncbi.nlm.nih.gov/protein/' + encodeURIComponent(ro.protein_id) +
        '" rel="noopener">' + esc(ro.protein_id) + '</a>' : '–') +
      (ro.locus_tag ? ' · ' + esc(ro.locus_tag) : '') +
      (ro.gene ? ' · gene <i>' + esc(ro.gene) + '</i>' : '') +
      ' · <span id="seqlen">…</span> residues</dd>');
    h.push('<dt>Position</dt><dd>' + fmt(ro.start) + '–' + fmt(ro.end) + ' on the ' +
      (ro.strand === 1 ? '+' : '−') + ' strand</dd>');
    h.push('</dl><dl class="kv">');
    h.push('<dt>Assigned type</dt><dd><a href="' + href('/cluster/' + ro.ro_cluster + '.html') +
      '"><b>' + esc(gene) + '</b></a> <span class="note">' + esc(ro.ro_cluster) +
      ', RO group ' + esc(ro.ro_group) + '</span></dd>');
    h.push('<dt>Variant</dt><dd>' + (ro.leaf_id
      ? '<a href="' + href('/leaf/' + leafSlug(ro.leaf_id) + '.html') + '">' + esc(ro.leaf_id) +
        '</a> <span class="note">' + fmt(ro.leaf_size) + ' members</span>' : '–') + '</dd>');
    h.push('<dt>Profile evidence</dt><dd>' + num(ro.hmm_score || 0, 1) + ' bits · E = ' +
      expo(ro.hmm_evalue || 0, 1) + ' · model coverage ' +
      num(ro.model_coverage || 0, 2) + '</dd>');
    h.push('<dt>Metal centres</dt><dd><span class="chip ' +
      (ro.rieske_intact ? 'chip--ok' : 'chip--bad') + '">Rieske [2Fe-2S] C-x-H…C-x-x-H ' +
      (ro.rieske_intact ? 'complete' : 'incomplete') + '</span> <span class="chip ' +
      (ro.catalytic_intact ? 'chip--ok' : 'chip--bad') +
      '">mononuclear Fe(II) 2-His-1-carboxylate ' +
      (ro.catalytic_intact ? 'complete' : 'incomplete') + '</span></dd>');
    if (etc && etc.etc_profile) {
      h.push('<dt>Electron transport</dt><dd class="mono">' + esc(etc.etc_profile) +
        ' <span class="note">' + (etc.in_operon ? 'a partner is inside the operon'
          : 'partners only outside the operon') + '</span></dd>');
    }
    h.push('<dt>Assignment class</dt><dd><span class="chip chip--' +
      esc(ro.assignment_class || '') + '">' + esc(ro.assignment_class || '–') + '</span>' +
      (ro.core_identity !== null && ro.core_identity !== undefined
        ? ' <span class="note">alignment identity to nearest reference ' +
          num(ro.core_identity, 2) + '</span>' : '') + '</dd>');
    h.push('</dl></div>');

    // --- komsuluk haritasi
    h.push('<h2>Genomic neighbourhood</h2>');
    h.push('<div class="legend">' +
      '<span><i style="background:#c0392b"></i>this alpha subunit</span>' +
      '<span><i style="background:#e67e22"></i>another alpha</span>' +
      '<span><i style="background:#2980b9"></i>beta subunit</span>' +
      '<span><i style="background:#27ae60"></i>ferredoxin</span>' +
      '<span><i style="background:#f39c12"></i>reductase</span>' +
      '<span><i style="background:#8e44ad"></i>regulator</span>' +
      '<span><i style="background:#2c3e50"></i>transposase or IS</span>' +
      '<span><i style="background:#d35400"></i>ring cleavage</span>' +
      '<span><i style="background:#b2bec3"></i>other</span>' +
      '<span><i class="dash"></i>predicted operon</span></div>');
    h.push('<figure class="fig"><div class="archscroll">' +
      window.ROARneigh.neighborhoodSvg(ro, nb, M.component_colors, M.neighbor_window) +
      '</div><figcaption>Genes within 10 kb on either side, drawn to scale; arrow direction is ' +
      'the coding strand. The dashed outline marks the predicted operon, the run of same-strand ' +
      'genes separated by at most 150 base pairs. Hover any gene for its product, coordinates ' +
      'and the Pfam evidence behind its component label.</figcaption></figure>');

    // --- regulasyon
    if (reg && reg.architecture) {
      var ARCH = {
        divergent_regulator: '<b>Regulator in divergent orientation.</b> The regulator and this ' +
          'operon face each other across a single intergenic region, the arrangement in which ' +
          'one regulatory site can control both.',
        codirectional_regulator: 'A regulator lies upstream on the same strand, so the two may ' +
          'belong to one transcript rather than facing each other.',
        divergent_other: 'The upstream gene is transcribed away from the operon but is not a ' +
          'recognised regulator.',
        codirectional_other: 'The upstream gene is on the same strand and is not a recognised ' +
          'regulator, so the operon may extend further than predicted.'
      };
      h.push('<h2>Operon start and promoter region</h2><dl class="kv">');
      h.push('<dt>Architecture</dt><dd>' +
        (ARCH[reg.architecture] || 'No upstream gene inside the analysed window.') + '</dd>');
      h.push('<dt>Operon 5\' end</dt><dd>' + fmt(reg.operon_5p) + '</dd>');
      h.push('<dt>Intergenic region</dt><dd>' + fmt(reg.intergenic_bp) + ' bp <span class="note">' +
        'the promoter lies in this interval; its sequence is not analysed here</span></dd>');
      h.push('<dt>Upstream gene</dt><dd>' + esc(reg.upstream_product || '–') +
        (reg.upstream_family && reg.upstream_family !== 'regulator_unclassified'
          ? ' <span class="chip">' + esc(reg.upstream_family) + ' family</span>' : '') +
        ' <span class="note">' + (reg.upstream_divergent ? 'opposite strand' : 'same strand') +
        '</span></dd></dl>');
      h.push('<p class="note"><a href="' + href('/atlas/regulation.html') +
        '">How this is measured across the family</a></p>');
    }

    // --- operon
    if (op) {
      h.push('<h2>Predicted operon</h2>');
      h.push('<p>' + fmt(op.n_genes) + ' genes spanning ' + fmt(op.start) + '–' + fmt(op.end) +
        ':' + ['beta', 'ferredoxin', 'reductase'].map(function (k) {
          var on = op['has_' + k];
          return ' <span class="chip ' + (on ? 'chip--ok' : 'chip--bad') + '">' + k + ' ' +
            (on ? 'present' : 'not detected') + '</span>';
        }).join('') + '</p>');
      h.push('<p class="note">Within 10 kb but outside the operon: beta ' +
        (op.nearby_beta ? 'yes' : 'no') + ', ferredoxin ' +
        (op.nearby_ferredoxin ? 'yes' : 'no') + ', reductase ' +
        (op.nearby_reductase ? 'yes' : 'no') + '.</p>');

      /* BIR GEN DIZISI, PAYLASILAN BIR ISLEV DEMEK DEGILDIR.
         Jinja sablonu (templates/ro.html) ile AYNI metni uretir; iki render
         yolu ayrisirsa yayimlanan sayfa ile gelistirme sunucusu farkli sey
         soyler. Olculdu: dogrulanmis girislerin %11,2'sinde operonda hic RO
         bileseni yok, %18,8'inde operon yalnizca alfadan olusuyor. */
      if (!(op.has_beta || op.has_ferredoxin || op.has_reductase)) {
        var near = nb.filter(function (n) {
          return !n.in_operon && Math.abs(n.distance) <= 5000 &&
            ['regulator', 'transporter', 'transposon', 'ring_cleavage']
              .indexOf(n.category) >= 0;
        });
        var box = '<div class="caveat"><b>No electron-transport partner in this operon.</b> ' +
          'The genes below are the run of same-strand genes next to the alpha subunit, ' +
          'separated by at most 150 base pairs. That is a statement about position and ' +
          'strand, not about function: adjacency is how the operon is called, so when no ' +
          'beta subunit, ferredoxin or reductase is present the run may hold genes that ' +
          'share nothing with this enzyme but their location.';
        if (near.length) {
          box += '<br><br><b>What does sit nearby, outside the operon.</b> ' +
            near.map(function (n) {
              return esc(String(n.product).slice(0, 52)) + ' <span class="note">(' +
                esc(n.category) + ', ' + fmt(Math.abs(n.distance)) + ' bp, ' +
                (n.same_strand ? 'same strand' : 'opposite strand') + ')</span>';
            }).join('; ') +
            '. A gene on the opposite strand between the alpha subunit and a further gene ' +
            'rules the two out of one transcript, which is why same-strand neighbours ' +
            'beyond it are not in the operon even though they may belong to the same ' +
            'functional module.';
        }
        h.push(box + '</div>');
      }
      h.push('<div id="opsvg"></div>');
      h.push('<table class="data compact"><thead><tr><th class="num">Position</th>' +
        '<th>Component by sequence</th><th>Category by annotation</th><th>Product</th>' +
        '<th class="num">Gap to previous gene</th></tr></thead><tbody>' +
        opg.map(function (g) {
          return '<tr class="' + (g.position === 0 ? 'hl' : '') + '">' +
            '<td class="num">' + g.position + '</td><td>' +
            (g.component && g.component !== 'none'
              ? '<span class="comp comp--' + esc(g.component) + '">' + esc(g.component) + '</span>'
              : '<span class="note">not detected</span>') +
            '</td><td>' + esc(g.category) + '</td><td class="small">' + esc(g.product) +
            '</td><td class="num">' + (g.position !== 0 ? fmt(g.gap_bp) + ' bp' : '') +
            '</td></tr>';
        }).join('') + '</tbody></table>');
    }

    // --- butun komsular
    h.push('<h2>All neighbouring genes</h2><div class="tablewrap">' +
      '<table class="data compact"><thead><tr><th class="num">Gene offset</th>' +
      '<th class="num">Distance</th><th>Strand</th><th>Component by sequence</th>' +
      '<th>Category</th><th>Product</th><th>Protein</th></tr></thead><tbody>' +
      nb.map(function (n) {
        return '<tr class="' + (n.in_operon ? 'op' : '') + '">' +
          '<td class="num">' + n.gene_offset + '</td>' +
          '<td class="num">' + fmt(n.distance) + ' bp</td>' +
          '<td>' + (n.strand === 1 ? '+' : '−') +
          (n.same_strand ? ' <span class="note">same</span>' : '') + '</td><td>' +
          (n.component && n.component !== 'none'
            ? '<span class="comp comp--' + esc(n.component) + '" title="' + esc(n.evidence) +
              '">' + esc(n.component) + '</span>' : '') +
          '</td><td>' + esc(n.category) + '</td>' +
          '<td class="small">' + esc(n.product) +
          (n.gene ? ' <i>(' + esc(n.gene) + ')</i>' : '') + '</td>' +
          '<td class="small">' + (n.protein_id
            ? '<a href="https://www.ncbi.nlm.nih.gov/protein/' + encodeURIComponent(n.protein_id) +
              '" rel="noopener">' + esc(n.protein_id) + '</a>' : '') + '</td></tr>';
      }).join('') + '</tbody></table></div>');

    // --- dizi (ayri parcadan)
    h.push('<h2>Sequence</h2><p class="note">The sequence is below for copying. Per-entry FASTA ' +
      'files are not published on the static site; <a href="' +
      href('/download/all_confirmed.fasta') + '">the bulk FASTA</a> contains every entry, keyed ' +
      'by the identifier shown above.</p><pre class="seq" id="seqbox">loading…</pre>');

    root.innerHTML = h.join('');
    document.title = (ro.protein_id || ro.locus_tag || ro.candidate_id) + ' · ROAR-DB';

    if (rec.opsvg) {
      getJSON('/data/opsvg/' + rec.opsvg.slice(0, M.shard_hex) + '.json').then(function (b) {
        var svg = b[rec.opsvg];
        if (!svg) return;
        var t = (shard.t || {})[ro.candidate_id] || ['promoter', 'regulator'];
        svg = svg.replace('__BP__', esc(t[0])).replace('__PROD__', esc(t[1]));
        var el = document.getElementById('opsvg');
        if (el) el.innerHTML = svg;
      }).catch(function () { /* sema yoksa sayfa yine dogru */ });
    }

    getJSON('/data/seq/' + sh + '.json').then(function (b) {
      var seq = b[ro.candidate_id] || '';
      var lines = [];
      for (var i = 0; i < seq.length; i += 60) lines.push(seq.slice(i, i + 60));
      document.getElementById('seqbox').textContent =
        '>' + ro.candidate_id + ' ' + (ro.protein_id || '') + ' ' + ro.product +
        ' [' + ro.organism + ']\n' + lines.join('\n') + '\n';
      var sl = document.getElementById('seqlen');
      if (sl) sl.textContent = String(seq.length);
    }).catch(function () {
      document.getElementById('seqbox').textContent = '(sequence not available)';
    });
  }
})();
