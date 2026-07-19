FROM rocker/r-ver:4.4.3

LABEL org.opencontainers.image.source="https://github.com/yharigaya/detectgxt" \
      org.opencontainers.image.description="Minimal R image with detectgxt and runtime dependencies"

# snpStats, a Bioconductor dependency, may compile C code that needs zlib headers
# when source packages are installed in minimal Linux environments.
RUN apt-get update -qq \
    && apt-get install -y --no-install-recommends zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

RUN Rscript -e 'install.packages(c("remotes", "BiocManager"), repos = "https://cloud.r-project.org")'

COPY . /tmp/detectgxt

RUN Rscript -e 'options(repos = BiocManager::repositories()); remotes::install_local("/tmp/detectgxt", dependencies = NA, upgrade = "never", build_vignettes = FALSE)' \
    && Rscript -e 'library(detectgxt)' \
    && rm -rf /tmp/detectgxt

CMD ["R"]
