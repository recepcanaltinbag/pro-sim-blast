/* ROAR-DB · structures.js
 *
 * SMILES -> SVG, tek yerden. Daha once ayni betik index.html ve cluster.html
 * icinde IKI KEZ duruyordu; biri duzeltilip oteki unutulabilirdi.
 *
 * NEDEN DEGISTI. Cizim 190x150 piksellik sabit bir tuvale yapiliyordu ve bag
 * uzunlugu sabitti, bu yuzden buyuk molekullerde atomlar UST USTE biniyordu:
 * 5,5'-dehidrodivanillat (43 karakterlik SMILES), dehidroabietat ve
 * 3-ketosteroid en kotuleriydi. Ayrica cizilen SVG'ye viewBox yazilmiyordu,
 * dolayisiyla CSS onu kucultunce oran bozuluyordu.
 *
 * Cozum iki parcali: molekul GENIS bir tuvale cizilir, sonra SVG'ye viewBox
 * verilip sabit en/boy nitelikleri silinir. Boylece olcegi tarayici yapar,
 * oranlar korunur ve kucuk molekul de buyuk molekul de kutusunu duzgun
 * doldurur.
 */
(function () {
  'use strict';

  var els = document.querySelectorAll('svg.struct[data-smiles]');
  if (!els.length) return;

  /* Kutuphane yuklenmediyse her kutuya aciklama konur. Sessiz bos kare
     birakmak en kotusu: okuyucu cizimin eksik oldugunu anlamaz. */
  if (typeof SmilesDrawer === 'undefined') {
    els.forEach(function (el) { placeholder(el, 'structure library did not load'); });
    return;
  }

  var dark = matchMedia('(prefers-color-scheme: dark)').matches;

  /* Tuval, gosterildigi kutudan iki kat buyuk. Bag uzunlugu da buyudu:
     amac molekule YER acmak, kucuk cizip buyutmek degil. */
  var drawer = new SmilesDrawer.SvgDrawer({
    width: 380, height: 300,
    bondThickness: 1.4, bondLength: 26,
    atomVisualization: 'default', terminalCarbons: false,
    explicitHydrogens: false, compactDrawing: false, padding: 14,
    fontSizeLarge: 7, fontSizeSmall: 5
  });

  function placeholder(el, message) {
    el.insertAdjacentHTML('afterend',
      '<div class="struct struct--none">' + message + '</div>');
    el.remove();
  }

  /*  Cizimden SONRA olceklenebilir hale getirir. SmilesDrawer en/boy
      niteliklerini piksel olarak yaziyor; viewBox olmadan CSS ile kucultmek
      ici bozuyor.                                                          */
  function makeResponsive(el) {
    var w = parseFloat(el.getAttribute('width')) || 380;
    var h = parseFloat(el.getAttribute('height')) || 300;
    if (!el.getAttribute('viewBox')) el.setAttribute('viewBox', '0 0 ' + w + ' ' + h);
    el.setAttribute('preserveAspectRatio', 'xMidYMid meet');
    el.removeAttribute('width');
    el.removeAttribute('height');
  }

  els.forEach(function (el) {
    var smiles = el.dataset.smiles;
    try {
      SmilesDrawer.parse(smiles, function (tree) {
        try {
          drawer.draw(tree, el, dark ? 'dark' : 'light', false);
          makeResponsive(el);
        } catch (e) {
          placeholder(el, 'structure could not be drawn');
        }
      }, function () {
        /* Ayristirma hatasi: SMILES'in kendisi sorunlu, bu bir kuratorluk
           bulgusudur ve gizlenmemeli. */
        placeholder(el, 'SMILES could not be parsed');
      });
    } catch (e) {
      placeholder(el, 'structure could not be drawn');
    }
  });
})();
