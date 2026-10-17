FROM python:3.12-slim
WORKDIR /app
RUN apt-get update && apt-get install -y --no-install-recommends \
    ca-certificates hisat2 samtools sra-toolkit \
    && rm -rf /var/lib/apt/lists/*
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt
COPY rnako.py pipeline.py webapp.py transfers.py ./
COPY test_core.py test_integration.py test_web.py test_resume.py ./
COPY static ./static
EXPOSE 8080
CMD ["python", "-u", "/app/webapp.py", "--host", "0.0.0.0", "--data", "/workdata", "--references", "/references", "--jobs", "/results/web"]
