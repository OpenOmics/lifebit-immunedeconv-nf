# syntax=docker/dockerfile:1.6
#
# Lightweight image for T-cell (CD4/CD8) deconvolution from bulk RNA-seq
# using the R immunedeconv package.
#
# Base: rocker/r-ver (Debian slim + pinned R) — much smaller than tidyverse/verse.
# Final image size: ~1.1-1.4 GB (vs ~3-4 GB for bioconductor/rstudio bases).
#
# Bundled methods (no external license required):
#   quanTIseq, MCP-counter, EPIC, xCell, ABIS, ConsensusTME, ESTIMATE
# NOT bundled (require user-supplied binary / license):
#   CIBERSORT, CIBERSORTx  — mount CIBERSORT.R + LM22.txt at runtime with -v
#
# Gene ID conversion:
#   Handles Ensembl IDs (ENSG.../ENSMUSG...) and Entrez IDs on input, mapping
#   to HGNC (human) or MGI (mouse) symbols via org.Hs.eg.db / org.Mm.eg.db.
#
# Build:
#   docker build -t immunedeconv:0.1 .
# Run (interactive):
#   docker run --rm -it -v $PWD:/data immunedeconv:0.1 R
# Run (script):
#   docker run --rm -v $PWD:/data immunedeconv:0.1 \
#       Rscript /opt/scripts/run_deconv.R --tpm /data/tpm.tsv --method quantiseq --out /data/out.tsv

ARG R_VERSION=4.4.1
FROM rocker/r-ver:${R_VERSION} AS builder

ENV DEBIAN_FRONTEND=noninteractive \
    R_LIBS_USER=/usr/local/lib/R/site-library \
    MAKEFLAGS="-j4"

# System deps needed to compile the R package graph (curl/xml/ssl/png/etc.)
RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential \
        gfortran \
        libcurl4-openssl-dev \
        libssl-dev \
        libxml2-dev \
        libgit2-dev \
        libfontconfig1-dev \
        libharfbuzz-dev \
        libfribidi-dev \
        libfreetype6-dev \
        libpng-dev \
        libtiff-dev \
        libjpeg-dev \
        zlib1g-dev \
        libbz2-dev \
        liblzma-dev \
        libicu-dev \
        libglpk-dev \
        pandoc \
        git \
        ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# Install BiocManager + remotes, then immunedeconv from Bioconductor.
# immunedeconv pulls in EPIC, MCPcounter, xCell, ConsensusTME, ESTIMATE,
# quanTIseq (bundled), ABIS, etc. as dependencies.
RUN R -e "install.packages(c('BiocManager','remotes','optparse','readr'), \
          repos='https://cloud.r-project.org', Ncpus=4)" \
 && R -e "BiocManager::install(version='3.19', ask=FALSE, update=FALSE)" \
 && R -e "BiocManager::install(c('immunedeconv','AnnotationDbi','org.Hs.eg.db','org.Mm.eg.db'), \
          ask=FALSE, update=FALSE, Ncpus=4)"

# -------- Final stage: strip build toolchain --------
FROM rocker/r-ver:${R_VERSION}

ENV DEBIAN_FRONTEND=noninteractive \
    R_LIBS_USER=/usr/local/lib/R/site-library

# Runtime-only libs (must match what the compiled R packages link against).
# rocker/r-ver:4.4.x is Ubuntu 22.04 (Jammy) based -> use jammy package names.
RUN apt-get update && apt-get install -y --no-install-recommends \
        libcurl4 \
        libssl3 \
        libxml2 \
        libgit2-1.1 \
        libfontconfig1 \
        libharfbuzz0b \
        libfribidi0 \
        libfreetype6 \
        libpng16-16 \
        libtiff5 \
        libjpeg-turbo8 \
        zlib1g \
        libbz2-1.0 \
        liblzma5 \
        libicu70 \
        libglpk40 \
        procps \
        ca-certificates \
    && rm -rf /var/lib/apt/lists/* \
    && rm -rf /tmp/* /var/tmp/*

# Copy the R library from the builder (contains immunedeconv + all deps)
COPY --from=builder /usr/local/lib/R/site-library /usr/local/lib/R/site-library

# Small driver script so the container is directly usable from Nextflow
RUN mkdir -p /opt/scripts /data
COPY run_deconv.R /opt/scripts/run_deconv.R
RUN chmod +x /opt/scripts/run_deconv.R

# Non-root user (nf-core / lifebit convention)
RUN groupadd -r deconv && useradd -r -g deconv -d /home/deconv -m -s /bin/bash deconv \
 && chown -R deconv:deconv /data
USER deconv
WORKDIR /data

LABEL org.opencontainers.image.title="immunedeconv" \
      org.opencontainers.image.description="Bulk RNA-seq immune cell deconvolution (CD4/CD8 T cells and more) via the R immunedeconv package" \
      org.opencontainers.image.source="https://github.com/omnideconv/immunedeconv" \
      org.opencontainers.image.licenses="BSD-3-Clause"

CMD ["R", "--no-save"]
