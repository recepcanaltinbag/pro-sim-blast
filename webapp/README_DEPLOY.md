# ROAR-DB web — çalıştırma ve ücretsiz yayınlama

Uygulama `webapp/app.py` (FastAPI + Jinja2). Veri kaynağı `roar.sqlite`
(pipeline çıktısı, `build_operons.py` sonrası `ro.sequence`, `operon*`,
`neighbor_protein*` tabloları dolu olmalı).

## Yerelde

```bash
pip install -r webapp/requirements.txt       # fastapi, uvicorn, jinja2, python-multipart, httpx
cd webapp && uvicorn app:app --reload --port 8000
# http://localhost:8000     API dokümanı: /api/docs
```

Ortam değişkenleri: `ROAR_DB`, `ROAR_HMM`, `ROAR_MOTIF_HMM`, `ROAR_ECOLOGY`
(varsayılanlar üst dizine bakar).

## Seçenek A — GitHub Pages (ücretsiz, sunucusuz; önerilen)

Tüm gezinme, arama, küme/varyant/giriş sayfaları ve indirmeler statik olarak
üretilir. **Yalnızca "Classify a sequence" sayfası statikte yoktur** (hmmer
sunucu ister; bkz. Seçenek B).

```bash
cd webapp
python3 freeze.py --out site --base /RieskeDB   # --base: Pages'in yayınlanacağı repo adı
bash deploy_pages.sh                            # site/ → github.com/recepcanaltinbag/RieskeDB gh-pages dalı
```

Sonra GitHub → Settings → Pages → Branch: `gh-pages` / root. Adres:
`https://recepcanaltinbag.github.io/RieskeDB/`. (Kod `pro-sim-blast` reposunda kalır; Pages yalnızca üretilmiş `site/` içeriğini taşır.)

Özel alan adı ya da `<kullanıcı>.github.io` deposu kullanırsan `--base ""`.
Site boyutu ~1 GB sınırının altında kalır (bkz. `freeze_full.log`); büyürse
`--limit` ile giriş sayfası sayısı kısılabilir veya `download/` klasörü
Releases'e taşınabilir.

## Seçenek B — Hugging Face Spaces (ücretsiz, Docker; sınıflandırıcı dahil)

Tam uygulama (dizi sınıflandırma + JSON API) için. Ücretsiz CPU Space
(2 vCPU, 16 GB RAM) yeterli; 48 saat trafik olmazsa uyur, ilk istekte uyanır.

1. huggingface.co → New Space → SDK: **Docker**, görünürlük Public.
2. Space deposuna şunları koy (build bağlamı pipeline dizinidir):
   ```
   Dockerfile              ← webapp/Dockerfile kopyası (kökte olmalı)
   webapp/                 (app.py, templates/, static/, requirements.txt)
   ro_filter.py  ro_motif.py  cluster_ecology.csv
   RieskeDB71.hmm          (14 MB)
   ROs_71_Clean/ROmotif71.hmm
   roar.sqlite             (~200 MB → git lfs track "*.sqlite" "*.hmm")
   ```
   ```bash
   git lfs install && git lfs track "*.sqlite" "*.hmm"
   ```
3. Push. Space `PORT=7860`'ı otomatik kullanır. README.md başına
   `sdk: docker` ve `app_port: 7860` içeren YAML bloğu ekle.

Yerel test: `docker build -f webapp/Dockerfile -t roar-db . && docker run -p 7860:7860 roar-db`

## Seçenek C — Vercel / Netlify / Cloudflare Pages

Statik dışa aktarımı (`site/`) bunlara da yükleyebilirsin; GitHub Pages'e göre
artısı yok. Vercel'in Python serverless katmanı `hmmsearch` ikilisini
çalıştıramaz; sınıflandırıcı orada çalışmaz. Render/Fly.io gibi ücretsiz
konteyner katmanları Docker imajını çalıştırır ama 200 MB'lık DB'yi imaja
gömmek gerekir (yapılıyor) ve uyku/bant genişliği sınırları vardır.

## Veri güncellendiğinde

`roar.sqlite` değişince: statik site için `freeze.py` + `deploy_pages.sh`;
HF Space için yeni `roar.sqlite`'ı push et (LFS).
