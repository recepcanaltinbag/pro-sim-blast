/*  Uzun aciklamalari KIRPAR, silmez.  (2026-11-02)

    SORUN, OLCUMLE. Kullanici: "sayfayi ilk acan direkt bisuru yaziyla
    karsilasiyor ya tablonun altin bi yazi var inanilmaz uzun aciklama vs".
    Olculdu: site genelinde 45 kelimeden uzun 131 <figcaption>/<p class="note">
    blogu var ve bunlar toplam 11.481 kelime; en uzunu 283 kelime. Bu metnin
    HEPSI gerekli -- her biri bir olcumun sinirini ya da bir caveat'i yaziyor
    ve bu projede "yayinlanan her iddianin kosulu yazili" kurali var. Yani
    cozum silmek DEGIL, ilk bakista 3 satir gostermek ve gerisini istege
    bagli acmak.

    NEDEN line-clamp DEGIL. `-webkit-line-clamp` icin elemana
    `display:-webkit-box` vermek gerekiyor; bu bloklarin bir kismi blok
    seviyeli cocuk tasiyor (orn. reaksiyon semasindaki .rxnote) ve box
    modeli onlari bozuyor. `max-height` + `overflow:hidden` hangi cocuk
    olursa olsun calisir.

    NEDEN DOM CERRAHISI DEGIL. Ilk tasarim metni ilk cumleden kesip gerisini
    yeni bir elemana tasiyordu. O yol icerideki <a>, <b> ve <code>
    etiketlerini ortasindan bolebiliyor ve bozuk HTML uretiyor. Burada
    HICBIR dugum tasinmiyor: yalnizca bir sinif ve bir dugme ekleniyor.

    JS YOKSA metin TAM gorunur -- yani bozulma yonunde degil, eski davranisa
    dogru bozuluyor.  */
(function () {
  "use strict";
  var MIN_WORDS = 45;

  function words(el) {
    return (el.textContent || "").trim().split(/\s+/).length;
  }

  function clamp(el, i) {
    if (el.closest(".longnote")) return;          // ic ice sarma yok
    el.classList.add("longnote");
    var btn = document.createElement("button");
    btn.type = "button";
    btn.className = "longnote__more";
    var id = el.id || ("longnote-" + i);
    el.id = id;
    btn.setAttribute("aria-controls", id);
    btn.setAttribute("aria-expanded", "false");
    btn.textContent = "Show the full note";
    btn.addEventListener("click", function () {
      var open = el.classList.toggle("is-open");
      btn.setAttribute("aria-expanded", open ? "true" : "false");
      btn.textContent = open ? "Show less" : "Show the full note";
    });
    el.insertAdjacentElement("afterend", btn);
  }

  function run() {
    var sel = "figcaption, p.note, .groupnote p, figcaption + .note";
    var seen = 0;
    document.querySelectorAll(sel).forEach(function (el) {
      if (el.dataset.noclamp !== undefined) return;
      if (words(el) < MIN_WORDS) return;
      clamp(el, seen++);
    });
  }

  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", run);
  } else {
    run();
  }
})();
