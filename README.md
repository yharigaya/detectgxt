# detectgxt

## About

DetectGxT is a method for detecting gene-by-treatment (GxT) interactions on molecular count phenotypes, such as RNA-seq and ATAC-seq data.
The key feature of this method is the use of nonlinear regression, which captures the relationship between genotypes and transformed molecular count phenotypes more accurately than a commonly used linear regression model.
See [`Get started`](https://yharigaya.github.io/detectgxt/articles/detectgxt.html) and [`Reference`](https://yharigaya.github.io/detectgxt/reference/index.html) for usage instructions.

Portions of the code and documentation were generated with the assistance of Claude Code (Anthropic) and Codex (OpenAI).
We have reviewed and verified the materials produced by these tools.

## Publication

Yuriko Harigaya,
Michael I. Love\*,
William Valdar\*.
"DetectGxT: detecting gene-by-treatment interactions on molecular count phenotypes accounting for allelic additivity."
(\* These authors contributed equally to this work.)

## Data input

DetectGxT takes as input a list of feature-SNP pairs for which significant genotype effects have already been identified in at least one of the control and treated conditions.
Individual genotype and phenotype data for these pairs in both conditions are also required.
The input set of feature-SNP pairs can be identified by a standard method, such as [MatrixEQTL](https://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/runit.html), [TensorQTL](https://github.com/broadinstitute/tensorqtl), and [limix_qtl](https://github.com/single-cell-genetics/limix_qtl/wiki).
See [`Get started`](https://yharigaya.github.io/detectgxt/articles/detectgxt.html) for the format of input data.

## Installation

The *detectgxt* R package can be installed from GitHub using *devtools*.

```r
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}
devtools::install_github("yharigaya/detectgxt")
```

### System prerequisites

*detectgxt* imports [snpStats](https://bioconductor.org/packages/snpStats/), a Bioconductor package that includes compiled C code. 
Most users will not need to install additional system libraries.

If installation fails while compiling a dependency from source with an error such as `zlib.h: No such file or directory`, install the zlib development headers for your system and then re-run `devtools::install_github()`:

- Debian/Ubuntu: `sudo apt-get install zlib1g-dev`
- RHEL/Fedora/CentOS: `sudo dnf install zlib-devel`
  or `sudo yum install zlib-devel`
- macOS: install the Xcode Command Line Tools with
  `xcode-select --install`; if using Homebrew R, also try `brew install zlib`

On Windows, no separate zlib installation is usually needed. 
If R asks to build packages from source, install the version of [Rtools](https://cran.r-project.org/bin/windows/Rtools/) matching your R version.

If `snpStats` is not installed automatically, install it from Bioconductor first:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install("snpStats")
devtools::install_github("yharigaya/detectgxt")
```

### Using a container

A container image with *detectgxt* and its runtime dependencies is published to GitHub Container Registry.

With Docker:

```sh
docker pull ghcr.io/yharigaya/detectgxt:latest
docker run --rm ghcr.io/yharigaya/detectgxt:latest Rscript -e 'library(detectgxt)'
```

On HPC systems where Apptainer/Singularity is used instead of Docker, pull the same image directly:

```sh
apptainer pull detectgxt.sif docker://ghcr.io/yharigaya/detectgxt:latest
apptainer exec detectgxt.sif Rscript -e 'library(detectgxt)'
```

This image provides *detectgxt* and its runtime dependencies only; projects with additional package or data requirements should build their own image from `ghcr.io/yharigaya/detectgxt:latest`.
