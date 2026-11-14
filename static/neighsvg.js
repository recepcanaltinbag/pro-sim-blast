/* Komsuluk gen haritasi -- app.py icindeki neighborhood_svg'nin PORTU.
 *
 * NEDEN PORT. Bu SVG giris sayfasinin 6,9 KB'ini tutuyordu ve her girisin
 * koordinatlari farkli oldugu icin sunucuda hesaplanip HTML'e gomuldugunde
 * 19.109 kez tekrarlaniyordu. Burada veriden cizilir.
 *
 * NEDEN GUVENLI. Port, verify_entry_render.py ile sinanir: ayni girdi icin
 * Python'un urettigi dizgi ile buranin urettigi dizgi BIREBIR karsilastirilir,
 * ornekleme degil, butun girisler uzerinde. Ayrismasi halinde sinama kirilir.
 *
 * Python tarafiyla ayni kalmasi gereken ayrintilar:
 *   - sayi bicimi: Python f"{x:.1f}" -> burada fmt1()
 *   - ok geometrisi, renk secimi, etiket kirpmasi (14 karakter)
 *   - HTML kacisi (html.escape): & < > ve " (quote=True varsayilani)
 */
(function (global) {
  'use strict';

  function fixed(x, d) {
    /* Python'un "%.<d>f" bicimiyle BIREBIR ayni cikti.
     *
     * NEDEN toFixed YETMIYOR. Python yarimi CIFTE yuvarlar (round-half-even),
     * JS toFixed ise yukari. Olculdu: komsuluk SVG'sinde 6.000 giriste 10'u
     * (216.25 -> Python "216.2", toFixed "216.3"), giris sayfasinda ise
     * kaplama yuzdesi (98.5 -> Python "98", toFixed "99") ayrisiyordu.
     * Bu degerler double olarak TAM temsil edilebildigi icin esitlik
     * gercekten olusuyor, yuvarlama hatasi degil.
     *
     * Yontem: tam ondalik acilimi toFixed(20) ile al, hedef haneden SONRAKI
     * kisim tam olarak "5" + sifirlar ise cifte yuvarla, degilse normal.
     * Cikti dizgisi elle kurulur ki son adimda yeniden kayan nokta hatasi
     * girmesin.
     */
    var neg = x < 0 || (x === 0 && 1 / x < 0);
    var s = Math.abs(x).toFixed(20);
    var dot = s.indexOf('.');
    var keep = s.slice(0, dot) + s.substr(dot + 1, d);     // d hane tutulur
    var rest = s.slice(dot + 1 + d);
    var n = parseInt(keep, 10);
    if (/^50*$/.test(rest)) {
      if (n % 2 !== 0) n += 1;
    } else if (rest && rest.charCodeAt(0) >= 53) {
      n += 1;
    }
    var str = String(n);
    if (d === 0) return (neg ? '-' : '') + str;
    while (str.length <= d) str = '0' + str;
    return (neg ? '-' : '') + str.slice(0, str.length - d) + '.' + str.slice(str.length - d);
  }

  function fmt1(x) { return fixed(x, 1); }

  function esc(s) {
    // html.escape(s) varsayilani quote=True: & < > " ' kacirir.
    return String(s === null || s === undefined ? '' : s)
      .replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;')
      .replace(/"/g, '&quot;').replace(/'/g, '&#x27;');
  }

  function neighborhoodSvg(ro, genes, colors, window_) {
    var W = 1000, H = 110, Y = 55;
    var left = Math.max(0, ro.start - window_);
    var right = ro.end + window_;
    var drawn = genes.filter(function (g) { return !g.spans_origin; });
    if (drawn.length) {
      drawn.forEach(function (g) {
        if (g.start < left) left = g.start;
        if (g.end > right) right = g.end;
      });
    }
    var span = Math.max(1, right - left);
    function sx(x) { return (x - left) / span * W; }

    function arrow(start, end, strand, color, label, title, bold, box) {
      var x1 = sx(start), x2 = sx(end);
      if (x2 - x1 < 4) x2 = x1 + 4;
      var h = bold ? 18 : 14;
      var head = Math.min(8, (x2 - x1) / 2);
      var pts;
      if (strand === -1) {
        pts = [[x2, Y - h], [x1 + head, Y - h], [x1, Y], [x1 + head, Y + h], [x2, Y + h]];
      } else {
        pts = [[x1, Y - h], [x2 - head, Y - h], [x2, Y], [x2 - head, Y + h], [x1, Y + h]];
      }
      var poly = pts.map(function (p) { return fmt1(p[0]) + ',' + fmt1(p[1]); }).join(' ');
      var stroke = bold ? '#111' : '#555';
      var out = '<polygon points="' + poly + '" fill="' + color + '" stroke="' + stroke +
        '" stroke-width="' + (bold ? 1.5 : 0.8) + '"><title>' + title + '</title></polygon>';
      if (box) {
        out += '<rect x="' + fmt1(x1 - 1) + '" y="' + (Y - h - 5) + '" width="' +
          fmt1(x2 - x1 + 2) + '" height="' + (2 * h + 10) +
          '" fill="none" stroke="#2d3436" stroke-dasharray="3,2" stroke-width="0.8"/>';
      }
      if (label && x2 - x1 > 24) {
        var ty = strand !== -1 ? Y + h + 14 : Y - h - 6;
        out += '<text x="' + fmt1((x1 + x2) / 2) + '" y="' + ty +
          '" font-size="9" text-anchor="middle" fill="#2d3436">' +
          String(label).slice(0, 14) + '</text>';
      }
      return out;
    }

    var parts = ['<svg viewBox="0 0 ' + W + ' ' + H + '" width="100%" preserveAspectRatio="none" ' +
      'xmlns="http://www.w3.org/2000/svg" class="genes">',
      '<line x1="0" y1="' + Y + '" x2="' + W + '" y2="' + Y + '" stroke="#999" stroke-width="1"/>'];

    drawn.forEach(function (g) {
      var comp = g.component || 'none';
      var key = comp !== 'none' ? comp : (g.category || 'other');
      var color = colors[key] || '#b2bec3';
      var label = g.gene || (comp !== 'none' ? comp : '');
      var title = esc(g.product + ' | ' + g.start + '-' + g.end + ' (' +
        (g.strand === 1 ? '+' : '-') + ') | component: ' + comp +
        ' | category: ' + g.category + ' | distance ' + g.distance + ' bp');
      parts.push(arrow(g.start, g.end, g.strand, color, label, title, false, !!g.in_operon));
    });
    parts.push(arrow(ro.start, ro.end, ro.strand, colors.alpha,
      ro.gene || 'RO α', esc(ro.product + ' (this entry)'), true, true));

    var kb = 2000;
    parts.push('<line x1="' + fmt1(W - 10 - kb / span * W) + '" y1="' + (H - 6) + '" x2="' +
      (W - 10) + '" y2="' + (H - 6) + '" stroke="#333" stroke-width="2"/>' +
      '<text x="' + fmt1(W - 10 - kb / span * W / 2) + '" y="' + (H - 9) +
      '" font-size="9" text-anchor="middle">2 kb</text>');
    parts.push('</svg>');
    return parts.join('');
  }

  var api = { neighborhoodSvg: neighborhoodSvg, esc: esc, fmt1: fmt1, fixed: fixed };
  if (typeof module !== 'undefined' && module.exports) module.exports = api;
  else global.ROARneigh = api;
})(typeof window !== 'undefined' ? window : this);
