/*  Plotly'nin yatay cubuk grafiklerinde sol kenar boslugu SABITTI (l: 180).
    360 px genisliginde bir telefonda bu, grafige 180 px birakiyor ve kategori
    etiketleri kirpiliyordu (gorsel denetim, 2026-11-01). `responsive: true`
    bunu cozmez: genisligi izler, kenar boslugunu izlemez.

    plotMargin(wide) kenar boslugunu ekran genisliginden TURETIR, ve
    bindPlotMargin ekran dondugunde yeniden duzenler. Grafigi kuran sayfa
    yalnizca `margin: plotMargin(180)` yazar.                              */
function plotMargin(wide, opts) {
  const w = window.innerWidth || 1024;
  // Grafigin kendisine en az yarim ekran kalmali; dar ekranda etiketler
  // kisalir ama GORUNUR, uzun etiket hover'da tam okunuyor.
  const l = Math.max(44, Math.min(wide, Math.round(w * 0.38)));
  return Object.assign({l: l, r: 24, t: 10, b: 64}, opts || {});
}

function bindPlotMargin(id, wide, opts) {
  let last = null;
  const apply = () => {
    const el = document.getElementById(id);
    if (!el || typeof Plotly === 'undefined') return;
    const m = plotMargin(wide, opts);
    if (last === m.l) return;          // gereksiz yeniden cizim yok
    last = m.l;
    Plotly.relayout(el, {'margin.l': m.l});
  };
  let timer = null;
  window.addEventListener('resize', () => {
    clearTimeout(timer);
    timer = setTimeout(apply, 150);
  });
  apply();
}
