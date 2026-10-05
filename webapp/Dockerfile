# ROAR-DB web application.
# Build context must be the pipeline directory (parent of webapp/):
#     docker build -f webapp/Dockerfile -t roar-db .
#     docker run -p 7860:7860 roar-db
# Port 7860 is what Hugging Face Spaces expects.
FROM python:3.11-slim

RUN apt-get update && apt-get install -y --no-install-recommends hmmer \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app
COPY webapp/requirements.txt /app/webapp/requirements.txt
RUN pip install --no-cache-dir -r /app/webapp/requirements.txt

# pipeline modules the classifier re-uses + models + data
COPY ro_filter.py ro_motif.py cluster_ecology.csv /app/
COPY RieskeDB71.hmm /app/RieskeDB71.hmm
COPY ROs_71_Clean/ROmotif71.hmm /app/ROs_71_Clean/ROmotif71.hmm
COPY roar.sqlite /app/roar.sqlite
COPY webapp /app/webapp

# hmmpress speeds up hmmsearch on the 71-model file
RUN hmmpress -f /app/RieskeDB71.hmm >/dev/null

ENV ROAR_DB=/app/roar.sqlite \
    ROAR_HMM=/app/RieskeDB71.hmm \
    ROAR_MOTIF_HMM=/app/ROs_71_Clean/ROmotif71.hmm \
    ROAR_ECOLOGY=/app/cluster_ecology.csv \
    PORT=7860

# non-root user (required by Hugging Face Spaces)
RUN useradd -m -u 1000 user && chown -R user /app
USER user
EXPOSE 7860
WORKDIR /app/webapp
CMD ["sh", "-c", "uvicorn app:app --host 0.0.0.0 --port ${PORT} --workers 2"]
