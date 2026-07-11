# detectgxt

<!-- [ClassifyGxT](https://yharigaya.github.io/classifygxt) -->

## About

DetectGxT is a method for detecting gene-by-treatment (GxT) interactions on molecular count phenotypes, such as RNA-seq and ATAC-seq data.
<!--It takes as input a list of feature-SNP pairs (or SNPs) for which significant genotype effects have already been identified in at least one of the control and treated conditions.
Individual genotype and phenotype data for these pairs in both conditions are also required.-->
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

The *detectgxt* R package can be installed using *devtools*.

```
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
devtools::install_github("yharigaya/detectgxt")
```
