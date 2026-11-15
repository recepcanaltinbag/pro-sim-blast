/* md5 -- YALNIZCA parca anahtari icin (guvenlik amaci YOK).
 *
 * build_entry_data.py parcalari md5(candidate_id)'nin ilk iki hex hanesine
 * gore ayiriyor. Tarayici ayni anahtari hesaplayamazsa giris bulunamaz, yani
 * bu fonksiyonun Python hashlib.md5 ile BIREBIR ayni olmasi sart.
 * verify_shard_key.py bunu 19.109 gercek kimlik uzerinde siniyor.
 */
(function (global) {
  'use strict';
  var md5 = (function () {
    function s(x, c) { return (x << c) | (x >>> (32 - c)); }
    function a(q, a_, b, x, s_, t) {
      a_ = (((a_ + q) | 0) + ((x + t) | 0)) | 0; return (((s(a_, s_)) + b) | 0);
    }
    function ff(a_, b, c, d, x, s_, t) { return a(((b & c) | (~b & d)), a_, b, x, s_, t); }
    function gg(a_, b, c, d, x, s_, t) { return a(((b & d) | (c & ~d)), a_, b, x, s_, t); }
    function hh(a_, b, c, d, x, s_, t) { return a((b ^ c ^ d), a_, b, x, s_, t); }
    function ii(a_, b, c, d, x, s_, t) { return a((c ^ (b | ~d)), a_, b, x, s_, t); }
    function toBlocks(str) {
      var utf = unescape(encodeURIComponent(str));
      var n = utf.length, blks = new Array(((n + 8) >> 6) + 1).fill(0).map(function () { return new Array(16).fill(0); });
      for (var i = 0; i < n; i++) blks[i >> 6][(i % 64) >> 2] |= utf.charCodeAt(i) << ((i % 4) * 8);
      blks[n >> 6][(n % 64) >> 2] |= 0x80 << ((n % 4) * 8);
      blks[blks.length - 1][14] = n * 8;
      return blks;
    }
    return function (str) {
      var blks = toBlocks(str), A = 1732584193, B = -271733879, C = -1732584194, D = 271733878;
      var K = [7, 12, 17, 22, 5, 9, 14, 20, 4, 11, 16, 23, 6, 10, 15, 21];
      for (var i = 0; i < blks.length; i++) {
        var x = blks[i], oa = A, ob = B, oc = C, od = D;
        A = ff(A, B, C, D, x[0], 7, -680876936); D = ff(D, A, B, C, x[1], 12, -389564586);
        C = ff(C, D, A, B, x[2], 17, 606105819); B = ff(B, C, D, A, x[3], 22, -1044525330);
        A = ff(A, B, C, D, x[4], 7, -176418897); D = ff(D, A, B, C, x[5], 12, 1200080426);
        C = ff(C, D, A, B, x[6], 17, -1473231341); B = ff(B, C, D, A, x[7], 22, -45705983);
        A = ff(A, B, C, D, x[8], 7, 1770035416); D = ff(D, A, B, C, x[9], 12, -1958414417);
        C = ff(C, D, A, B, x[10], 17, -42063); B = ff(B, C, D, A, x[11], 22, -1990404162);
        A = ff(A, B, C, D, x[12], 7, 1804603682); D = ff(D, A, B, C, x[13], 12, -40341101);
        C = ff(C, D, A, B, x[14], 17, -1502002290); B = ff(B, C, D, A, x[15], 22, 1236535329);
        A = gg(A, B, C, D, x[1], 5, -165796510); D = gg(D, A, B, C, x[6], 9, -1069501632);
        C = gg(C, D, A, B, x[11], 14, 643717713); B = gg(B, C, D, A, x[0], 20, -373897302);
        A = gg(A, B, C, D, x[5], 5, -701558691); D = gg(D, A, B, C, x[10], 9, 38016083);
        C = gg(C, D, A, B, x[15], 14, -660478335); B = gg(B, C, D, A, x[4], 20, -405537848);
        A = gg(A, B, C, D, x[9], 5, 568446438); D = gg(D, A, B, C, x[14], 9, -1019803690);
        C = gg(C, D, A, B, x[3], 14, -187363961); B = gg(B, C, D, A, x[8], 20, 1163531501);
        A = gg(A, B, C, D, x[13], 5, -1444681467); D = gg(D, A, B, C, x[2], 9, -51403784);
        C = gg(C, D, A, B, x[7], 14, 1735328473); B = gg(B, C, D, A, x[12], 20, -1926607734);
        A = hh(A, B, C, D, x[5], 4, -378558); D = hh(D, A, B, C, x[8], 11, -2022574463);
        C = hh(C, D, A, B, x[11], 16, 1839030562); B = hh(B, C, D, A, x[14], 23, -35309556);
        A = hh(A, B, C, D, x[1], 4, -1530992060); D = hh(D, A, B, C, x[4], 11, 1272893353);
        C = hh(C, D, A, B, x[7], 16, -155497632); B = hh(B, C, D, A, x[10], 23, -1094730640);
        A = hh(A, B, C, D, x[13], 4, 681279174); D = hh(D, A, B, C, x[0], 11, -358537222);
        C = hh(C, D, A, B, x[3], 16, -722521979); B = hh(B, C, D, A, x[6], 23, 76029189);
        A = hh(A, B, C, D, x[9], 4, -640364487); D = hh(D, A, B, C, x[12], 11, -421815835);
        C = hh(C, D, A, B, x[15], 16, 530742520); B = hh(B, C, D, A, x[2], 23, -995338651);
        A = ii(A, B, C, D, x[0], 6, -198630844); D = ii(D, A, B, C, x[7], 10, 1126891415);
        C = ii(C, D, A, B, x[14], 15, -1416354905); B = ii(B, C, D, A, x[5], 21, -57434055);
        A = ii(A, B, C, D, x[12], 6, 1700485571); D = ii(D, A, B, C, x[3], 10, -1894986606);
        C = ii(C, D, A, B, x[10], 15, -1051523); B = ii(B, C, D, A, x[1], 21, -2054922799);
        A = ii(A, B, C, D, x[8], 6, 1873313359); D = ii(D, A, B, C, x[15], 10, -30611744);
        C = ii(C, D, A, B, x[6], 15, -1560198380); B = ii(B, C, D, A, x[13], 21, 1309151649);
        A = ii(A, B, C, D, x[4], 6, -145523070); D = ii(D, A, B, C, x[11], 10, -1120210379);
        C = ii(C, D, A, B, x[2], 15, 718787259); B = ii(B, C, D, A, x[9], 21, -343485551);
        A = (A + oa) | 0; B = (B + ob) | 0; C = (C + oc) | 0; D = (D + od) | 0;
      }
      var hex = '';
      [A, B, C, D].forEach(function (v) {
        for (var j = 0; j < 4; j++) hex += ('0' + ((v >> (j * 8)) & 255).toString(16)).slice(-2);
      });
      return hex;
    };
  })();
  if (typeof module !== 'undefined' && module.exports) module.exports = md5;
  else global.ROARmd5 = md5;
})(typeof window !== 'undefined' ? window : this);
