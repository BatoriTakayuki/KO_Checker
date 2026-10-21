FROM python:3.12-slim
WORKDIR /app
RUN apt-get update && apt-get install -y --no-install-recommends \
    ca-certificates fastp hisat2 samtools sra-toolkit \
    && rm -rf /var/lib/apt/lists/*
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt
COPY rnako.py pipeline.py webapp.py transfers.py reporting.py storage.py region_loss.py structure.py igv_view.py ./
COPY test_core.py test_integration.py test_web.py test_resume.py test_reporting.py test_workflow.py test_storage.py test_region_loss.py test_structure.py ./
COPY static ./static
COPY test_igv.py ./
EXPOSE 8080
CMD ["python", "-u", "/app/webapp.py", "--host", "0.0.0.0", "--data", "/workdata", "--references", "/references", "--jobs", "/results/web"]
