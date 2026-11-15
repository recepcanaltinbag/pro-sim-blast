/* ROAR-DB · charts.js
 *
 * Grafikler icin ortak yardimcilar. Bu dosya, Plotly'nin betik adresi 404
 * dondugu ve HICBIR grafik gorunmedigi ortaya cikinca yazildi: adres
 * duzeltilince grafikler ilk kez goruldu ve bicimlerinin hic denetlenmemis
 * oldugu anlasildi. Uc kusur vardi ve ucu de burada toparlaniyor.
 *
 * 1. YUKSEKLIK. Uc grafik 12-14 cubuk icin 1500 piksele ayarlanmisti, yani
 *    cubuk basina 110 piksel. Cubuklar devasa, boslugu buyuk, sayfa
 *    gereksiz uzundu. Yukseklik artik cubuk SAYISINDAN hesaplanir.
 * 2. ACIKLAMA KUTUSU. `y:1.015` ile grafigin USTUNE konuyordu ve ust bosluk
 *    yalnizca 10 pikseldi; 13 girisli yatay bir aciklama iki-uc satira
 *    sariyor ve en ustteki cubuklarin uzerine biniyordu. Artik grafigin
 *    ALTINA konuyor ve alt bosluk ona gore ayrilyor.
 * 3. RENKLER. Eski palette dort kirmizi-kahve (#8c2f22, #b5503f, #c2703f) ve
 *    uc koyu mavi-gri (#1f4e79, #2d4356, #5b6670) vardi; yigili cubukta
 *    yan yana dusunce ayirt edilemiyorlardi. Yeni palette ARDIL renkler
 *    birbirinden uzak secildi, ve her dilime ince beyaz bir kenar konuyor:
 *    yigili cubukta sinirlari gormek rengi ayirt etmekten daha guvenilir.
 */
window.ROAR = window.ROAR || {};

/* Ardil ogeler arasinda ton, parlaklik ve doygunluk bakimindan belirgin fark
   olacak sekilde siralanmis 12 renk. Son oge "diger" icin notr gri. */
window.ROAR.PALETTE = [
  '#1f4e79', /* koyu mavi      */
  '#c2703f', /* turuncu-kahve  */
  '#1f7a4d', /* yesil          */
  '#8e5ea2', /* mor            */
  '#d9a13b', /* hardal         */
  '#2a8f87', /* turkuaz        */
  '#8c2f22', /* koyu kirmizi   */
  '#7fa650', /* zeytin yesili  */
  '#4a6fa5', /* orta mavi      */
  '#b5503f', /* kiremit        */
  '#5b6670', /* kursun grisi   */
  '#9e6b7d'  /* gul kurusu     */
];
window.ROAR.OTHER_COLOR = '#b9b4ab';

/* Yigili cubukta dilim siniri. Renk yakinligi kalirsa bile kenar ayirir. */
window.ROAR.sliceLine = function () {
  var dark = matchMedia('(prefers-color-scheme: dark)').matches;
  return { color: dark ? '#14181d' : '#ffffff', width: 0.8 };
};

/* Yatay cubuk grafigi icin yukseklik: cubuk basina sabit yer, arti eksen ve
   aciklama icin pay. Alt ve ust sinir, tek cubuklu ve cok cubuklu uc
   durumlarda grafigin okunur kalmasi icin. */
window.ROAR.barHeight = function (bars, legendRows) {
  var perBar = 30;
  var chrome = 70 + 20 * (legendRows || 1);
  return Math.max(260, Math.min(900, bars * perBar + chrome));
};

/* Aciklama kutusunu grafigin ALTINA koyar ve gereken alt boslugu dondurur. */
window.ROAR.legendBelow = function (entries, perRow) {
  var rows = Math.ceil(entries / (perRow || 6));
  return {
    legend: { orientation: 'h', y: -0.02, yanchor: 'top', x: 0, xanchor: 'left',
              font: { size: 10 } },
    bottomMargin: 52 + 18 * rows
  };
};
