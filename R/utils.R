#' Create a genotype matrix
#'
#' This function takes genotype data in PLINK format as input and
#' returns a \code{geno} data frame, which can be passed to
#' \code{\link{map_qtl}} using the \code{geno} argument. See the
#' package vignette for details.
#'
#' @export
#' @importFrom magrittr "%>%"
#' @importFrom snpStats read.plink
#' @importFrom methods as
#' @importFrom dplyr pull
#' @inheritParams map_qtl
#' @param plink.prefix A character string specifying the prefix of
#' PLINK files.
#' @param plink.suffix A character string specifying the suffix of
#' PLINK files.
#' @param chromosomes An integer vector specifying chromosomes.
#'
#' @return A data frame formatted as input for MatrixEQTL. The first
#' column contains SNP IDs. The rest of the columns contain genotypes
#' coded as (0, 1, 2). Each column corresponds to each sample.
#'
#' @examples
#' \dontrun{
#' candidate <- data.frame(feat_id = "f1", snp_id = "g1")
#' geno <- get_geno(
#'     candidate = candidate,
#'     chromosomes = 1,
#'     plink.prefix = "path/to/plink",
#'     plink.suffix = "bed"
#' )
#' }
get_geno <- function(candidate, chromosomes,
                     plink.prefix, plink.suffix) {

    select.snps <- candidate %>% pull("snp_id")
    geno.list <- vector("list", length(chromosomes))

    for (i in seq_along(chromosomes)) {
        chr <- chromosomes[[i]]

        plink <- paste0(plink.prefix, ".", chr, ".", plink.suffix)
        tmp <- read.plink(bed=plink)
        map.snps <- rownames(tmp$map)
        selected.snps <- intersect(select.snps, map.snps)
        if (length(selected.snps) == 0) {
            next
        }
        tmp <- read.plink(
            bed=plink, select.snps=selected.snps)
        geno <- as(tmp$genotypes, "numeric") %>% t
        geno.list[[i]] <- geno

    }

    geno.list <- Filter(Negate(is.null), geno.list)
    if (length(geno.list) == 0) {
        return(data.frame(snp_id=character(0), stringsAsFactors=FALSE))
    }

    geno.matrix <- do.call("rbind", geno.list)

    n.sub <- ncol(geno.matrix)
    output <- geno.matrix %>%
        as.data.frame %>%
        mutate(snp_id=rownames(geno.matrix)) %>%
        `[`(, c(n.sub + 1, seq_len(n.sub))) %>%
        `rownames<-`(NULL)

    output

}

#' Pre-process molecular count data
#'
#' This function takes raw molecular count data as input and returns a
#' \code{pheno} data frame, which can be passed to
#' \code{\link{map_qtl}} using the \code{pheno} argument.
#'
#' @export
#' @importFrom magrittr "%>%"
#' @importFrom stats lm residuals prcomp
#' @importFrom dplyr left_join
#' @inheritParams map_qtl
#' @param count A matrix containing molecular counts (non-negative
#' integers). The row and column names correspond to feature IDs and
#' sample IDs, respectively.
#' @param num.pc An integer specifying the number of molecular
#' phenotype principal components to control for.
#' @param list A logical scalar. If \code{TRUE}, the PCA result,
#' i.e., the output from the \code{prcomp} function in the stats
#' package, will be returned as an element of a list.
#'
#' @return A data frame containing pre-processed phenotype data or a
#' list containing the following elements:
#' \itemize{
#' \item{\code{pheno} - A data frame containing pre-processed
#' phenotype data.}
#' \item{\code{pca} - A list object from the \code{prcomp} function in
#' the stats package.}
#' }
#'
#' @examples
#' anno <- data.frame(
#'     sample = paste0("s", 1:20),
#'     subject = rep(paste0("d", 1:10), each = 2),
#'     condition = rep(c(0, 1), 10)
#' )
#' count <- matrix(
#'     rpois(100, lambda = 100), nrow = 5,
#'     dimnames = list(paste0("f", 1:5), paste0("s", 1:20))
#' )
#' pheno <- preprocess_pheno(count = count, anno = anno, num.pc = 3)
preprocess_pheno <- function(count, anno, covar=NULL,
                             num.pc=10, list=FALSE) {

    get_residual <- function(y, d) {
        y <- data.frame(sample=names(y), y=as.numeric(y))
        d <- y %>%
            left_join(d, by="sample")
        d <- d[, colnames(d) != "sample", drop=FALSE]
        residuals(lm(y ~ ., data=d))
    }

    # scale and transform data according to Palowitch et al. (2018)
    n.feat <- nrow(count)
    n.sample <- ncol(count)
    if (is.null(colnames(count)) || any(colnames(count) == "")) {
        stop("`count` must have column names")
    }
    if (anyDuplicated(colnames(count))) {
        stop("column names of `count` must be unique")
    }
    if (!all(c("sample", "condition") %in% colnames(anno))) {
        stop("`anno` must contain 'sample' and 'condition'")
    }
    if (anyDuplicated(anno$sample)) {
        stop("sample IDs in `anno` must be unique")
    }
    if (!all(colnames(count) %in% anno$sample)) {
        stop("column names of `count` must match samples in `anno`")
    }
    if (!is.numeric(num.pc) || length(num.pc) != 1 || is.na(num.pc) ||
        num.pc < 0 || num.pc != as.integer(num.pc)) {
        stop("`num.pc` must be a non-negative integer scalar")
    }
    tot <- colSums(count)
    xbar <- sum(tot) / (n.feat * n.sample)
    scaled <- sweep(count, 2, tot, "/") * n.feat * xbar
    transformed <- log(scaled + 1)

    d <- anno[
        match(colnames(transformed), anno$sample),
        c("sample", "condition")]

    if (!is.null(covar)) {
        if (!("sample" %in% colnames(covar))) {
            colnames(covar)[1] <- "sample"
        }
        if (anyDuplicated(covar$sample)) {
            stop("sample IDs in `covar` must be unique")
        }
        if (!all(colnames(count) %in% covar$sample)) {
            stop("column names of `count` must match samples in `covar`")
        }
        d <- d %>% left_join(covar, by="sample")
    }

    residualized <- apply(
        transformed, 1, get_residual, d=d) %>%
        `rownames<-`(colnames(transformed))

    pca <- prcomp(residualized)
    num.pc.use <- min(num.pc, ncol(pca$x))
    if (num.pc.use == 0) {
        pc <- data.frame(row.names=colnames(transformed))
    } else {
        pc <- pca$x[colnames(transformed), seq_len(num.pc.use), drop=FALSE] %>%
            as.data.frame
    }

    d2 <- data.frame(sample=colnames(transformed), pc) %>%
        `rownames<-`(NULL)

    if (!is.null(covar)) {
        d2 <- d2 %>%
            left_join(covar, by="sample")
    }

    res <- apply(transformed, 1, get_residual, d=d2) %>%
        t %>%
        as.data.frame %>%
        `colnames<-`(colnames(transformed)) %>%
        mutate(feat_id=rownames(transformed)) %>%
        `[`(, c(n.sample + 1, seq_len(n.sample))) %>%
        `rownames<-`(NULL)

    if (list) {
        res <- list(pheno=res, pca=pca)
    }

    res

}

#' Simulate molecular count data
#'
#' This function generates a matrix of negative-binomial count data
#' using the nonlinear GxT parametrization consistent with the models
#' fit by \code{\link{map_qtl}} when \code{fn = "nonlinear"}.
#' Each feature is associated with its own SNP. When \code{geno} is
#' \code{NULL}, a separate genotype vector is simulated per feature.
#'
#' @export
#' @importFrom magrittr "%>%"
#' @importFrom stats rnorm rnbinom runif rbinom
#' @inheritParams make_sim_data
#' @param n.feat An integer specifying the number of features (e.g.,
#' genes and ATAC-seq peaks).
#' @param n.sample An integer specifying the number of samples.
#' @param geno An n.feat x n.sub matrix of genotype values (0, 1, 2).
#' Column names must match \code{unique(anno$subject)}. Row \code{i}
#' gives the subject-level genotypes for feature \code{i}. When
#' \code{NULL} (the default), a separate genotype vector is simulated
#' for each feature using \code{lb.maf} and \code{ub.maf}.
#' @param sd A numeric vector of length three specifying the standard
#' deviations of the genotype, treatment, and interaction effects,
#' respectively. Each effect is drawn independently per feature from
#' Normal(0, sd). Default \code{c(0, 0, 0)} gives no effects. For each
#' effect, exactly one of the corresponding elements of \code{sd} and
#' \code{coef} must be non-\code{NA}.
#' @param coef A numeric vector of length three specifying fixed values
#' for the genotype, treatment, and interaction effects, respectively.
#' For each effect, exactly one of the corresponding elements of
#' \code{sd} and \code{coef} must be non-\code{NA}.
#' @param lb.maf A scalar specifying the lower bound of the MAF used
#' when simulating genotypes (\code{geno = NULL}).
#' @param ub.maf A scalar specifying the upper bound of the MAF used
#' when simulating genotypes (\code{geno = NULL}).
#' @param filter.geno A logical scalar. When \code{TRUE} and
#' \code{geno = NULL}, genotype simulation for each feature is repeated
#' until all three levels (0, 1, 2) are present.
#' @param intercept.mean A scalar specifying the mean of the per-feature
#' intercept \code{b0} drawn from \code{Normal(intercept.mean,
#' intercept.sd)}.
#' @param intercept.sd A scalar specifying the standard deviation of the
#' per-feature intercept.
#' @param disp.mean.rel A function that maps \code{exp(b0)} (the
#' baseline mean count) to the negative-binomial dispersion
#' (\code{size = 1 / dispersion} in \code{rnbinom}).
#' @param size.factors A numeric vector of length \code{n.sample}
#' specifying multiplicative size factors applied to the count-scale
#' mean.
#'
#' @return A matrix of non-negative integers with \code{n.feat} rows
#' and \code{n.sample} columns. Row names are \code{"f1", "f2", ...}
#' and column names are the sample IDs in \code{anno$sample}.
#'
#' @examples
#' anno <- data.frame(
#'     sample = paste0("s", 1:20),
#'     subject = rep(paste0("d", 1:10), each = 2),
#'     condition = rep(c(0, 1), 10)
#' )
#' counts <- make_count_data(
#'     anno = anno, n.feat = 5, n.sample = 20,
#'     sd = c(0, 0, 0)
#' )
make_count_data <- function(anno, n.feat, n.sample,
                            geno=NULL,
                            sd=c(0, 0, 0),
                            coef=rep(NA, 3),
                            lb.maf=0.05,
                            ub.maf=0.5,
                            filter.geno=FALSE,
                            intercept.mean=4,
                            intercept.sd=2,
                            disp.mean.rel=function(x) 4/x + 0.1,
                            size.factors=rep(1, n.sample)) {

    if (nrow(anno) != n.sample) {
        stop("`n.sample` must match the number of rows in `anno`")
    }
    if (length(size.factors) != n.sample) {
        stop("`size.factors` must have length `n.sample`")
    }

    if (length(sd) != 3) stop("the length of `sd` must be three")
    if (length(coef) != 3) stop("the length of `coef` must be three")
    effect.names <- c("genotype", "treatment", "interaction")
    for (k in seq_len(3)) {
        if ((is.na(coef[k]) & is.na(sd[k])) |
            (!is.na(coef[k]) & !is.na(sd[k])))
            stop(paste(
                "only one of `coef` and `sd` must be specified for the",
                effect.names[k], "effect"))
    }
    if (!is.null(geno)) {
        if (!is.matrix(geno))
            stop("`geno` must be a matrix")
        if (nrow(geno) != n.feat)
            stop("`geno` must have n.feat rows")
        if (!all(colnames(geno) %in% unique(anno$subject)))
            stop("column names of `geno` must match the subjects in `anno`")
    }

    b0.vec <- rnorm(n.feat, intercept.mean, intercept.sd)
    dispersion <- disp.mean.rel(exp(b0.vec))

    get_bvec <- function(k) {
        if (!is.na(sd[k])) rnorm(n.feat, 0, sd[k]) else rep(coef[k], n.feat)
    }
    b1.vec <- get_bvec(1)
    b2.vec <- get_bvec(2)
    b3.vec <- get_bvec(3)

    subject <- anno$subject
    t.vec <- anno$condition
    n.sub <- length(unique(subject))
    sub.names <- unique(subject)
    max.try <- 10000L

    if (filter.geno && is.null(geno) && n.sub < 3) {
        stop("`filter.geno` requires at least three subjects")
    }

    count.data <- matrix(NA_integer_, nrow=n.feat, ncol=n.sample)

    for (i in seq_len(n.feat)) {

        if (is.null(geno)) {
            maf <- runif(1, lb.maf, ub.maf)
            if (!filter.geno) {
                geno.sub <- rbinom(n=n.sub, size=2, prob=maf)
            } else {
                geno.vec <- rep(0, 3)
                num.try <- 0L
                while (any(geno.vec == 0)) {
                    num.try <- num.try + 1L
                    geno.sub <- rbinom(n=n.sub, size=2, prob=maf)
                    geno.vec <- factor(geno.sub, levels=c(0, 1, 2)) %>%
                        table %>%
                        as.numeric
                    if (num.try >= max.try && any(geno.vec == 0)) {
                        stop("unable to simulate genotypes satisfying `filter.geno`")
                    }
                }
            }
            names(geno.sub) <- sub.names
            g.vec <- geno.sub[subject]
        } else {
            g.vec <- geno[i, subject]
        }

        if (any(is.na(g.vec)))
            stop("some subjects in `anno` have no matching entry in `geno`")

        param.i <- c(b0=b0.vec[i], b1=b1.vec[i], b2=b2.vec[i], b3=b3.vec[i])
        log.mu.i <- get_mean(
            g=g.vec, t=t.vec, param=param.i, m=c(1L, 1L, 1L), fn="nonlinear")
        mu.i <- exp(log.mu.i) * size.factors
        count.data[i, ] <- rnbinom(n.sample, mu=mu.i, size=1/dispersion[i])
    }

    mode(count.data) <- "integer"
    rownames(count.data) <- paste0("f", seq_len(n.feat))
    colnames(count.data) <- anno$sample

    count.data

}

#' Filter genotypes
#'
#' This function takes a \code{candidate} data frame (see
#' \code{\link{map_qtl}}) as input and removes SNPs for which either
#' of the homozygotes is missing in at least one of the control and
#' treated conditions.
#'
#' @export
#' @inheritParams map_qtl
#'
#' @return A list of data frames containing the following elements:
#' \itemize{
#' \item{\code{included} - A data frame containing filtered data.}
#' \item{\code{excluded} - A data frame containing excluded
#' feature-SNP pairs.}
#' }
#'
#' @examples
#' candidate <- data.frame(
#'     feat_id = c("f1", "f2"),
#'     snp_id = c("g1", "g2")
#' )
#' geno <- data.frame(
#'     snp_id = c("g1", "g2"),
#'     d1 = c(0, 1), d2 = c(1, 2), d3 = c(2, 0)
#' )
#' anno <- data.frame(
#'     sample = paste0("s", 1:6),
#'     subject = rep(c("d1", "d2", "d3"), each = 2),
#'     condition = rep(c(0, 1), 3)
#' )
#' filter_geno(candidate = candidate, geno = geno, anno = anno)
filter_geno <- function(candidate, geno, anno=NULL) {

    if (!is.null(anno)) {
        req <- c("sample", "subject", "condition")
        if (!all(req %in% colnames(anno))) {
            stop("`anno` must contain 'sample', 'subject', and 'condition'")
        }
    }

    logic.vec <- rep(TRUE, nrow(candidate))
    for (i in seq_len(nrow(candidate))) {

        snp.id <- candidate[i, 2]

        if (!(snp.id %in% geno[, 1])) {
            warning(paste(snp.id, "not found in the genotype matrix"))
            logic.vec[i] <- FALSE
            next
        }
        g <- geno[geno[, 1] == snp.id, -1] %>%
            as.numeric
        names(g) <- colnames(geno)[-1]

        if (is.null(anno)) {
            if (sum(g == 0) == 0 || sum(g == 2) == 0) {
                logic.vec[i] <- FALSE
            }
        } else {
            g.anno <- g[anno$subject]
            if (any(is.na(g.anno))) {
                warning(paste(
                    "subjects in `anno` not found in genotype columns for",
                    snp.id))
                logic.vec[i] <- FALSE
                next
            }
            d <- data.frame(condition=anno$condition, geno=g.anno)
            support.by.cond <- tapply(
                d$geno, d$condition,
                function(x) any(x == 0) & any(x == 2))
            if (!all(support.by.cond)) {
                logic.vec[i] <- FALSE
            }
        }
    }

    included <- candidate[logic.vec, ]
    excluded <- candidate[!logic.vec, ]

    list(included=included, excluded=excluded)

}

#' Nonlinear mean function for GxT model
#'
#' This function computes the expected log expression under the nonlinear
#' genotype-by-treatment model used by \code{\link{map_qtl}}.
#'
#' @export
#' @param g A vector of the genotype coded as (0, 1, 2) or the
#' imputation-based allelic dosage.
#' @param t A vector of treatment indicators (0 or 1).
#' @param b0 A scalar specifying the intercept.
#' @param b1 A scalar specifying the genotype effect.
#' @param b2 A scalar specifying the treatment effect.
#' @param b3 A scalar specifying the genotype-by-treatment interaction
#' effect.
#'
#' @return A numeric vector of expected log expression values.
get_mu <- function(g, t, b0, b1, b2, b3) {
    x1 <- (1 - g/2) * (1 - t)
    x2 <- (g/2) * (1 - t)
    x3 <- (1 - g/2) * t
    x4 <- (g/2) * t
    r.mu <- x1 +
        exp(2 * b1) * x2 +
        exp(b2) * x3 +
        exp(2 * b1 + b2 + 2 * b3) * x4
    b0 + log(r.mu)
}
