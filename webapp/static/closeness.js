/* ROAR-DB · closeness.js
 *
 * Kopruleyen karboksilat iliskisini TARAYICIDA yeniden hesaplar.
 *
 * NEDEN. Yayinlanan tek bir V degeri, hangi uyelerin sayildigina gore 0,74 ile
 * 1,00 arasinda degisiyor. Okuyucunun bunu kendi secimiyle gormesi, bir tabloya
 * bakmasindan daha ikna edici; ve sunucu olmadigi icin hesap burada yapilmak
 * zorunda. Veri satirlari dizgi TEKRARLAMAZ: her sutun bir sozluge indekstir,
 * boylece 11.422 satir 144 KB'a siginiyor.
 */
(function () {
  'use strict';
  var out = document.getElementById('cbx-out');
  if (!out || typeof CBX_URL === 'undefined') return;

  var DATA = null, tierBoxes = [], classBoxes = [];

  function chi2AndV(table) {
    /* Bos satir ve sutunlar ATILIR: tek bir bos hucre bolme hatasi verir ve
       tabakanin hic olculememesine yol acar. */
    table = table.filter(function (r) { return r.reduce(function (a, b) { return a + b; }, 0) > 0; });
    if (table.length < 2) return null;
    var cols = table[0].length, keep = [];
    for (var j = 0; j < cols; j++) {
      var s = 0;
      for (var i = 0; i < table.length; i++) s += table[i][j];
      if (s > 0) keep.push(j);
    }
    if (keep.length < 2) return null;
    table = table.map(function (r) { return keep.map(function (j) { return r[j]; }); });
    var rowSums = table.map(function (r) { return r.reduce(function (a, b) { return a + b; }, 0); });
    var colSums = table[0].map(function (_, j) {
      return table.reduce(function (a, r) { return a + r[j]; }, 0);
    });
    var n = rowSums.reduce(function (a, b) { return a + b; }, 0);
    if (!n) return null;
    var chi2 = 0;
    for (var a = 0; a < table.length; a++) {
      for (var b = 0; b < table[a].length; b++) {
        var e = rowSums[a] * colSums[b] / n;
        if (e > 0) chi2 += Math.pow(table[a][b] - e, 2) / e;
      }
    }
    var k = Math.min(table.length - 1, table[0].length - 1);
    return { chi2: chi2, n: n, v: Math.sqrt(chi2 / (n * k)),
             rows: table.length, cols: table[0].length };
  }

  function selected(boxes) {
    var picked = {};
    boxes.forEach(function (b) { if (b.checked) picked[b.value] = true; });
    return picked;
  }

  function render() {
    if (!DATA) return;
    var tiers = selected(tierBoxes), classes = selected(classBoxes);
    var groups = DATA.groups, nGroups = groups.length;
    /* satir = grup, sutun = Asp / Glu */
    var table = [];
    for (var i = 0; i < nGroups; i++) table.push([0, 0]);
    var kept = 0, glu = 0;
    DATA.rows.forEach(function (r) {
      var residue = DATA.residues[r[2]];
      if (residue !== 'D' && residue !== 'E') return;
      if (!tiers[DATA.tiers[r[3]]]) return;
      if (!classes[DATA.classes[r[4]]]) return;
      table[r[0]][residue === 'D' ? 0 : 1] += 1;
      kept += 1;
      if (residue === 'E') glu += 1;
    });
    var stat = chi2AndV(table);
    var head = '<table class="data compact"><thead><tr><th>RO group</th>' +
               '<th class="num">Aspartate</th><th class="num">Glutamate</th>' +
               '<th class="num">Glutamate share</th></tr></thead><tbody>';
    var body = '';
    for (var g = 0; g < nGroups; g++) {
      var tot = table[g][0] + table[g][1];
      if (!tot) continue;
      body += '<tr><td>group ' + groups[g] + '</td>' +
              '<td class="num">' + table[g][0].toLocaleString() + '</td>' +
              '<td class="num">' + table[g][1].toLocaleString() + '</td>' +
              '<td class="num">' + (100 * table[g][1] / tot).toFixed(1) + ' %</td></tr>';
    }
    var summary;
    if (!stat) {
      summary = '<p class="note"><b>Not computable.</b> The selection leaves fewer than two ' +
                'groups or only one residue, so there is no table to test.</p>';
    } else {
      var strength = stat.v >= 0.9 ? 'near-deterministic' :
                     stat.v >= 0.5 ? 'strong' : stat.v >= 0.25 ? 'moderate' :
                     stat.v >= 0.1 ? 'weak' : 'negligible';
      summary = '<p class="note"><b>Cramér&rsquo;s V = ' + stat.v.toFixed(3) + '</b> (' + strength +
                ') on ' + stat.n.toLocaleString() + ' entries, ' +
                (100 * glu / kept).toFixed(1) + ' % glutamate, χ² = ' + stat.chi2.toFixed(0) +
                ' over a ' + stat.rows + ' by ' + stat.cols + ' table.' +
                (stat.n < 200 ? ' <b>Treat this with care:</b> below a few hundred entries the ' +
                 'figure moves a lot with small changes in the selection.' : '') + '</p>';
    }
    out.innerHTML = summary + '<div class="tablewrap">' + head + body + '</tbody></table></div>';
  }

  function boxes(container, values, labels) {
    return values.map(function (v) {
      var id = container.id + '-' + v;
      var wrap = document.createElement('label');
      wrap.className = 'cbx__item';
      var input = document.createElement('input');
      input.type = 'checkbox'; input.value = v; input.checked = true; input.id = id;
      input.addEventListener('change', render);
      wrap.appendChild(input);
      wrap.appendChild(document.createTextNode(' ' + (labels[v] || v.replace(/_/g, ' '))));
      container.appendChild(wrap);
      return input;
    });
  }

  fetch(CBX_URL).then(function (r) { return r.json(); }).then(function (d) {
    DATA = d;
    tierBoxes = boxes(document.getElementById('cbx-tiers'), d.tiers, {});
    classBoxes = boxes(document.getElementById('cbx-classes'), d.classes, {});
    document.getElementById('cbx-all').addEventListener('click', function () {
      tierBoxes.concat(classBoxes).forEach(function (b) { b.checked = true; });
      render();
    });
    document.getElementById('cbx-close').addEventListener('click', function () {
      tierBoxes.forEach(function (b) {
        b.checked = (b.value === 'characterized' || b.value === 'close_homolog');
      });
      classBoxes.forEach(function (b) { b.checked = true; });
      render();
    });
    render();
  }).catch(function () {
    out.innerHTML = '<p class="note">The per-entry table could not be loaded, so this ' +
                    'section cannot recompute anything. The tables above are unaffected.</p>';
  });
})();
