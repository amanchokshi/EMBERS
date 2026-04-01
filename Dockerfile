FROM python:3.13-slim

LABEL org.opencontainers.image.title="embers" \
      org.opencontainers.image.version="1.0.1" \
      org.opencontainers.image.authors="Aman Chokshi" \
      org.opencontainers.image.source="https://github.com/amanchokshi/embers" \
      org.opencontainers.image.description="EMBERS 1.0.1 processing environment for OzSTAR"

ENV DEBIAN_FRONTEND=noninteractive \
    PIP_NO_CACHE_DIR=1 \
    PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    pkg-config \
    libhdf5-dev \
    git \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /opt

ARG EMBERS_REF=1.0.1

RUN git clone --depth 1 --branch ${EMBERS_REF} https://github.com/amanchokshi/embers.git \
    && cd embers \
    && python -m pip install --upgrade pip setuptools wheel \
    && python -m pip install .

WORKDIR /work

CMD ["bash"]
