#' Generate data for simulation analysis
#'
#' This function generates simulated genotype and phenotype data for
#' all eight GxT model categories.
#'
#' @export
#' @import mvtnorm
#' @importFrom magrittr "%>%"
#' @importFrom stats runif rbinom rnorm
#' @importFrom purrr map_int map_dbl pluck
#' @importFrom dplyr mutate bind_rows bind_cols
#' @param num A named, ordered integer vector specifying the numbers of
#' simulations in the eight model categories.
#' @param fn A character string specifying the function. This must be
#' one of "nonlinear" and "linear", corresponding to nonlinear and
#' linear models, respectively.
#' @param lb.maf A scalar specifying the lower bound of MAF.
#' @param ub.maf A scalar specifying the upper bound of MAF.
#' @param filter.geno A logical scalar as to whether to ensure that
#' homozygotes for both alleles have at least one observation in both
#' control and treated conditions.
#' @param anno A data frame containing sample IDs (character strings
#' or integers), subject IDs (character strings or integers), and
#' treatment conditions (0 or 1). The columns must be named as "sample",
#' "subject", and "condition".
#' @param genotype A data frame containing the subject IDs (character
#' strings or integers) and the genotype coded as (0, 1, 2) or the
#' imputation-based allelic dosage in the first and second columns,
#' respectively. The columns must be named as "subject" and "g".
#' @param sd A numeric vector of length three specifying the standard
#' deviations of the genotype, treatment, and interaction effects,
#' respectively. For each effect, exactly one of the corresponding
#' elements of \code{sd} and \code{coef} must be non-\code{NA}.
#' @param coef A numeric vector of length three specifying fixed values
#' for the genotype, treatment, and interaction effects, respectively.
#' For each effect, exactly one of the corresponding elements of
#' \code{sd} and \code{coef} must be non-\code{NA}.
#' @param b0 A scalar specifying the intercept.
#' @param sigma A scalar specifying the residual error standard deviation.
#' @param ranef A logical scalar as to whether to include random
#' effect.
#' @param sigma.u A scalar specifying the random intercept standard
#' deviation. If \code{ranef} is \code{TRUE}, this is set to 1
#' by default.
#' @param kinship A matrix containing pairwise genetic relatedness
#' between subjects. The row and column names must match the set of
#' unique elements of the "subject" column in "anno" in the
#' corresponding order (i.e., \code{unique(anno$subject)}). This is
#' set to \code{NULL} by default, in which case the identity matrix is
#' used.
#' @param seed An integer specifying a seed for RNG.
#' @param output A character string specifying the output format. This
#' must be one of "list" and "data.frame".
#'
#' @return A list of lists containing:
#' \itemize{
#' \item{\code{y} - A vector of phenotypes.}
#' \item{\code{g} - A vector of genotypes.}
#' \item{\code{t} - A vector of treatment indicators.}
#' \item{\code{sample} - A vector of sequencing sample IDs.}
#' \item{\code{subject} - A vector of subject IDs.}
#' \item{\code{index} - An integer specifying one of the eight
#'     models.}
#' \item{\code{maf} - A scalar specifying the minor allele frequency
#'     used for generating the genotype data.}
#' \item{\code{beta} - A named numeric vector specifying the true
#'     coefficient values used for generating the phenotype data. The
#'     "b0" element represents the intercept. The "b1", "b2", and "b3"
#'     elements respectively represent the genotype, treatment, and
#'     interaction effect sizes. }
#' } or a list containing:
#' \itemize{
#' \item{\code{candidate} - A data frame containing feature IDs and
#' SNP IDs. See \code{\link{map_qtl}}. This also includes the integers
#' specifying models, the minor allele frequency, and the true
#' coefficient values.}
#' \item{\code{geno} - A data frame containing genotypes. See \code{\link{map_qtl}}.}
#' \item{\code{pheno} - A data frame containing phenotypes. See \code{\link{map_qtl}}.}
#' }
#'
#' @examples
#' anno <- data.frame(
#'     sample = paste0("s", 1:20),
#'     subject = rep(paste0("d", 1:10), each = 2),
#'     condition = rep(c(0, 1), 10)
#' )
#' dat <- make_sim_data(
#'     num = c(0, 0, 0, 0, 0, 0, 0, 1),
#'     anno = anno,
#'     fn = "nonlinear",
#'     coef = c(0.5, 0.3, 0.2)
#' )
make_sim_data <- function(num, anno, fn, genotype=NULL,
                          lb.maf=0.05, ub.maf=0.5,
                          filter.geno=FALSE,
                          sd=rep(NA, 3), coef=rep(NA, 3), b0=0, sigma=1,
                          ranef=FALSE, sigma.u=NULL, kinship=NULL,
                          seed=1, output="data.frame") {

    set.seed(seed)

    if (!(output %in% c("data.frame", "list"))) {
        stop("`output` must be one of `data.frame` and `list`")
    }

    if (!is.numeric(num) || length(num) != 8 || any(is.na(num)) ||
        any(num < 0) || any(num != as.integer(num))) {
        stop("`num` must be a non-negative integer vector of length 8")
    }

    if (!ranef) {
        if (!is.null(sigma.u)) {
            stop("`ranef` must be set to `TRUE` when `sigma.u` is supplied")
        }
        if (!is.null(kinship)) {
            stop("`ranef` must be set to `TRUE` when `kinship` is supplied")
        }
    }

    if (all(is.na(coef)) & all(is.na(sd))) {
        stop("at least one of the `coef` and `sd` needs to be specified.")
    }

    if (length(coef) != 3) {
        stop("the length of `coef` must be three")
    }

    if (length(sd) != 3) {
        stop("the length of `sd` must be three")
    }

    if ((is.na(coef[1]) & is.na(sd[1])) |
        ((!is.na(coef[1])) & (!is.na(sd[1])))) {
        stop("only one of the `coef` and `sd` needs to be specified for the genotype effect")
    }

    if ((is.na(coef[2]) & is.na(sd[2])) |
        ((!is.na(coef[2])) & (!is.na(sd[2])))) {
        stop("only one of the `coef` and `sd` needs to be specified for the treatment effect")
    }

    if ((is.na(coef[3]) & is.na(sd[3])) |
        ((!is.na(coef[3])) & (!is.na(sd[3])))) {
        stop("only one of the `coef` and `sd` needs to be specified for the interaction effect")
    }

    sample <- anno$sample
    subject <- anno$subject
    t <- anno$condition
    n.sample <- length(sample)
    n.sub <- length(unique(subject))

    if (ranef) {

        if (is.null(kinship)) {
            A <- diag(n.sub)
            rownames(A) <- colnames(A) <- unique(subject)
        } else {
            A <- kinship
            if (is.null(colnames(A)) | is.null(rownames(A))) {
                stop(paste(
                    "the kinship matrix must have",
                    "the row and column names"))
            }
            if (!all(colnames(A) == rownames(A))) {
                stop(paste(
                    "the row and column names of",
                    "the kinship matrix must be identical"))
            }
            if (!all(rownames(A) == unique(subject))) {
                stop(paste(
                    "the row names of",
                    "the kinship matrix must match",
                    "the subjects"))
            }
            asym <- max(abs(A - t(A)))
            if (asym > 1e-8) {
                stop(paste0(
                    "the kinship matrix is not symmetric ",
                    "(max |K - t(K)| = ", signif(asym, 3), "); ",
                    "please check the kinship matrix"),
                    call.=FALSE)
            }
        }

        if (is.null(sigma.u)) {
            sigma.u <- 1
        }

        Z <- matrix(0, nrow=n.sample, ncol=n.sub)
        for (i in seq_len(n.sample)) {
            for (j in seq_len(n.sub)) {
                if (subject[i] == unique(subject)[j]) Z[i, j] <- 1
            }
        }

        Sigma.u <- Z %*% A %*% t(Z)
        cov <-
            sigma.u^2 * Sigma.u +
            sigma^2 * diag(n.sample)
        if (!is.null(kinship)) {
            tryCatch(
                chol(cov),
                error=function(e) stop(
                    "the variance component sigma_u^2 * Z K Z' + sigma^2 * I ",
                    "is not positive definite; please check the kinship matrix",
                    call.=FALSE))
        }

    }

    model.name <- get_model_mat()

    data.list.list <- vector("list", length(num))
    max.try <- 10000L

    if (filter.geno && is.null(genotype)) {
        subject.by.condition <- split(subject, t)
        has.enough.subjects <- vapply(
            subject.by.condition,
            function(x) length(unique(x)) >= 2,
            logical(1))
        if (!all(has.enough.subjects)) {
            stop("`filter.geno` requires at least two subjects in each condition")
        }
    }

    for (j in seq_along(num)) {

        n.data <- num[j]
        if (n.data == 0) {
            data.list.list[[j]] <- NULL
            next
        }

        index <- j
        data.list <- vector("list", length(n.data))

        for (i in seq_len(n.data)) {

            m.vec <- as.numeric(model.name[index, ])

            if (is.null(genotype)) {
                maf <- runif(1, lb.maf, ub.maf)

                num.try <- 0L
                repeat {
                    num.try <- num.try + 1L
                    geno <- rbinom(n=n.sub, size=2, prob=maf)
                    names(geno) <- unique(subject)
                    if (!filter.geno) break
                    d.tmp <- merge(anno, geno, by.x="subject", by.y=0)
                    colnames(d.tmp)[ncol(d.tmp)] <- "geno"
                    if (all(tapply(d.tmp$geno, d.tmp$condition,
                                  function(g) any(g == 0) & any(g == 2)))) break
                    if (num.try >= max.try) {
                        stop("unable to simulate genotypes satisfying `filter.geno`")
                    }
                }

            } else if (is.data.frame(genotype)) {
                maf <- NA
                geno <- genotype$g
                names(geno) <- genotype$subject
            } else {
                stop("`genotype` has an incorrect format.")
            }

            d <- merge(anno, geno, by.x="subject", by.y=0) %>%
                `colnames<-`(
                    c("subject", "sample", "condition", "geno"))
            d <- d[match(anno$sample, d$sample), ]
            g <- d$geno

            if (m.vec[1] == 1) {
                b1 <- ifelse(
                    !is.na(sd[1]),
                    rnorm(n=1, mean=0, sd=sd[1]),
                    coef[1])
            } else {
                b1 <- 0
            }

            if (m.vec[2] == 1) {
                b2 <- ifelse(
                    !is.na(sd[2]),
                    rnorm(n=1, mean=0, sd=sd[2]),
                    coef[2])
            } else {
                b2 <- 0
            }

            if (m.vec[3] == 1) {
                b3 <- ifelse(
                    !is.na(sd[3]),
                    rnorm(n=1, mean=0, sd=sd[3]),
                    coef[3])
            } else {
                b3 <- 0
            }

            beta <- c(b0, b1, b2, b3) %>%
                `names<-`(c("b0", "b1", "b2", "b3"))

            y.mean <- get_mean(
                g=g, t=t, param=beta, m=m.vec, fn=fn)

            if (!ranef) {
                y <- y.mean +
                    rnorm(n=n.sample, mean=0, sd=sigma)
            } else {
                u.mean <- rep(0, n.sample)
                y <- y.mean +
                    rmvnorm(n=1, mean=u.mean, sigma=cov) %>%
                    as.numeric
            }

            data.list[[i]] <- list(
                y=y, g=g, t=t, sample=sample, subject=subject,
                index=index, maf=maf, beta=beta)

        }

        data.list.list[[j]] <- data.list
    }

    dat <- do.call("c", data.list.list)

    if (output == "data.frame") {
        dat <- dat %>%
            format_input(num=num, anno=anno)
    }

    dat

}

#' Get eigenvectors and eigenvalues of the covariance matrix
#'
#' This function computes eigenvectors and eigenvalues of the
#' covariance matrix, which are needed for running
#' \code{\link{map_qtl}} and \code{\link{map_qtl_for_each}} when
#' modeling a random effect. \code{kinship} must be set to the default
#' if a subject-specific random effect will be modeled. \code{kinship}
#' must be specified if a polygenic (kinship) random effect will be
#' modeled.
#' This is the same as the \code{get_tu_lambda} function from the classifygxt package.
#'
#' @export
#' @importFrom magrittr "%>%"
#' @inheritParams make_sim_data
#' @param kinship A matrix containing pairwise genetic relatedness
#' between subjects. The row and column names must match the set of
#' unique elements of the "subject" column in "anno" in the
#' corresponding order (i.e., \code{unique(anno$subject)}). This is
#' set to \code{NULL} by default, in which case the identity matrix is
#' used.
#'
#' @return A list object containing:
#' \itemize{
#' \item{\code{tU} - A matrix containing the transposed eigenvectors
#'     of the covariance matrix.}
#' \item{\code{lambda} - A vector containing the eigenvalues of the
#'     covariance matrix.}
#' }
#'
#' @examples
#' anno <- data.frame(
#'     sample = paste0("s", 1:20),
#'     subject = rep(paste0("d", 1:10), each = 2),
#'     condition = rep(c(0, 1), 10)
#' )
#' eig <- get_eigen(anno = anno)
get_eigen <- function(anno,
                      kinship=NULL) {

    t <- anno$condition
    n <- length(t)
    if (!is.numeric(t)) {
        stop("`anno` has an incorrect format")
    }

    subject <- anno$subject
    if (!(is.character(subject) | is.integer(subject))) {
        stop("`anno` has an incorrect format")
    }

    if (!(is.character(anno$sample) | is.integer(anno$sample))) {
        stop("`anno` has an incorrect format")
    }

    sub.name <- unique(subject)
    n.sub <- length(sub.name)
    sub.vec <- subject %>%
        factor(levels=sub.name) %>%
        as.numeric

    if (is.null(kinship)) {
        A <- diag(n.sub)
    } else if (is.matrix(kinship)) {
        A <- kinship
        if (is.null(colnames(A)) | is.null(rownames(A))) {
            stop(paste(
                "the kinship matrix must have",
                "the row and column names"))
        }
        if (!all(colnames(A) == rownames(A))) {
            stop(paste(
                "the row and column names of",
                "the kinship matrix must be identical"))
        }
        if (!all(rownames(A) == sub.name)) {
            stop(paste(
                 "the row names of",
                 "the kinship matrix must match",
                 "the subjects"))
        }
        asym <- max(abs(A - t(A)))
        if (asym > 1e-8) {
            stop(paste0(
                "the kinship matrix is not symmetric ",
                "(max |K - t(K)| = ", signif(asym, 3), "); ",
                "please check the kinship matrix"),
                call.=FALSE)
        }
    }

    Z <- matrix(0, nrow=n, ncol=n.sub)
    for (i in seq_len(n)) {
        for (j in seq_len(n.sub)) {
            if (sub.vec[i] == j)
                Z[i, j] <- 1
        }
    }

    Sigma.u <- Z %*% A %*% t(Z)

    eigen.sample.kernel <- eigen(Sigma.u, symmetric=TRUE)
    tU <- t(eigen.sample.kernel$vectors)
    lambda <- eigen.sample.kernel$values

    output <- list(tU=tU, lambda=lambda)
    output
}

#' Map marginal or response QTLs for feature-SNP pairs
#'
#' This function takes as input individual phenotype and genotype data
#' and tests for marginal and response QTLs.
#'
#' @export
#' @importFrom tibble as_tibble
#' @importFrom purrr map
#' @importFrom stats lm residuals
#' @inheritParams map_qtl_for_each
#' @param candidate A data frame containing feature IDs and SNP IDs.
#' The first column must contain feature IDs and the second column must
#' contain SNP IDs.
#' @param fdr A scalar specifying a false discovery rate (FDR)
#' threshold.
#' @param geno A data frame formatted as input for MatrixEQTL. The
#' first column must contain SNP IDs. The rest of the columns must
#' contain genotypes, which can be dosage or coded as (0, 1, 2). Each
#' column corresponds to each sample.
#' @param pheno A data frame formatted as input for MatrixEQTL. The
#' first column must contain feature IDs (e.g., genes and ATAC-peak
#' regions). The rest of the columns must contain pre-processed
#' molecular phenotypes. Each column corresponds to each sample.
#' @param anno A data frame containing sample IDs, subject IDs, and
#' treatment conditions. The columns must be named as "sample",
#' "subject", and "condition".
#' @param covar A data frame containing fixed effect covariates. The
#' first column must contain sample IDs. The rest of the columns must
#' contain covariate values. Each row corresponds to each sample. Each
#' column corresponds to each covariate.
#' @param ranef A logical scalar as to whether to include random
#' effect.
#' @param kinship A matrix containing pairwise genetic relatedness
#' between subjects. The row and column names must match the set of
#' unique elements of the "subject" column in "anno" in the
#' corresponding order (i.e., \code{unique(anno$subject)}). This is
#' set to \code{NULL} by default, in which case the identity matrix is
#' used.
#' @param filter.geno A logical scalar as to whether to ensure that
#' homozygotes for both alleles have at least one observation in both
#' control and treated conditions.
#' @param excluded A logical scalar as to whether to return a data
#' frame containing feature-SNP pairs that have been excluded from the
#' analysis. If \code{TRUE}, the output is a list of data frames.
#'
#' @return A data frame containing some of the following columns:
#' \itemize{
#' \item{\code{feat_id} - A character string specifying the feature ID.}
#' \item{\code{snp_id} - A character string specifying the SNP ID.}
#' \item{\code{beta0} - A scalar value representing the intercept
#' estimate.}
#' \item{\code{beta_g} - A scalar value representing the genotype
#' effect estimate based on the marginal model (when
#' \code{type = "marginal"}) or the interaction model (when
#' \code{type = "response"}).}
#' \item{\code{beta_t} - A scalar value representing the treatment
#' effect estimate based on the interaction model.}
#' \item{\code{beta_gxt} - A scalar value representing the interaction
#' effect estimate based on the interaction model.}
#' \item{\code{sigma} - A scalar value representing the estimate of
#' the residual error standard deviation.}
#' \item{\code{sigma_u} - A scalar value representing the estimate of
#' the random effect standard deviation.}
#' \item{\code{p_value} - A scalar value representing the nominal
#' p-value from LRT.}
#' \item{\code{return_code} - An integer: 0 if all optimizations
#' converged, \code{NA} otherwise.}
#' }
#'
#' @examples
#' anno <- data.frame(
#'     sample = paste0("s", 1:20),
#'     subject = rep(paste0("d", 1:10), each = 2),
#'     condition = rep(c(0, 1), 10)
#' )
#' dat <- make_sim_data(
#'     num = c(0, 0, 0, 0, 0, 0, 0, 1),
#'     anno = anno, fn = "nonlinear",
#'     coef = c(0.5, 0.3, 0.2)
#' )
#' result <- map_qtl(
#'     candidate = dat$candidate[1, ],
#'     fn = "nonlinear",
#'     geno = dat$geno,
#'     pheno = dat$pheno,
#'     anno = anno
#' )
map_qtl <- function(candidate,
                    fdr=NULL,
                    fn,
                    geno,
                    pheno,
                    anno,
                    covar=NULL,
                    scale=FALSE,
                    rint=FALSE,
                    ranef=FALSE,
                    kinship=NULL,
                    filter.geno=TRUE,
                    excluded=FALSE,
                    type="response",
                    seed=1,
                    epsilon=1e-2,
                    control=NULL,
                    gls=FALSE) {

    run_map_qtl <- function(i, candidate, fn, geno, pheno,
                            anno, covar, tu.lambda=NULL,
                            type, seed) {

        feat.id <- candidate[i, 1]
        snp.id <- candidate[i, 2]

        if (!(snp.id %in% geno[, 1])) {
            warning(paste(snp.id, "not found in the genotype matrix"))
            return(NULL)
        }
        g <- geno[geno[, 1] == snp.id, -1] %>%
            as.numeric
        names(g) <- colnames(geno)[-1]

        if (!(feat.id %in% pheno[, 1])) {
            warning(paste(feat.id, "not found in the phenotype matrix"))
            return(NULL)
        }
        y <- pheno[pheno[, 1] == feat.id, -1] %>%
            as.numeric
        names(y) <- colnames(pheno)[-1]

        d <- anno %>%
            `[`(, c("sample", "subject", "condition")) %>%
            merge(g, by.x="subject", by.y=0) %>%
            `colnames<-`(
                c("subject", "sample", "condition", "geno")) %>%
            merge(y, by.x="sample", by.y=0) %>%
            `colnames<-`(
                c("sample", "subject", "condition",
                  "geno", "pheno"))
        d <- d[match(anno$sample, d$sample), ]

        if (!all(d$subject == anno$subject)) {
            stop("subjects must match in `geno` and `anno`")
        }

        if (!all(d$sample == anno$sample)) {
            stop("samples must match in `pheno` and `anno`")
        }

        if (is.null(covar)) {
            input <- list(feat.id=feat.id, snp.id=snp.id,
                          y=d$pheno, g=d$geno,
                          t=d$condition, subject=d$subject)
        } else {
            d <- d %>%
                merge(covar, by.x="sample", by.y=1)
            d <- d[match(anno$sample, d$sample), ]

            regress <- d[, !colnames(d) %in%
                             c("sample", "subject", "condition", "geno")]
            d$residual <- residuals(lm(pheno ~ ., data=regress))
            input <- list(feat.id=feat.id, snp.id=snp.id,
                          y=d$residual, g=d$geno,
                          t=d$condition, subject=d$subject)
        }

        res <- map_qtl_for_each(
            input=input, fn=fn,
            scale=scale, rint=rint,
            tu.lambda=tu.lambda,
            type=type, seed=seed,
            epsilon=epsilon,
            control=control,
            gls=gls) %>%
            as_tibble %>%
            as.data.frame

        res

    }

    if (!is.null(fdr)) {
        if (!is.numeric(fdr) || length(fdr) != 1 || is.na(fdr) ||
            fdr < 0 || fdr > 1) {
            stop("`fdr` must be a numeric scalar between 0 and 1")
        }
    }

    if (filter.geno) {
        filtered <- filter_geno(candidate=candidate, geno=geno, anno=anno)
        candidate <- filtered$included
    }

    if (ranef) {

        if (!is.null(kinship)) {
            if (!all(!is.na(
                      match(rownames(kinship),
                            unique(anno$subject))))) {
                stop("subjects must match in `kinship` and `anno`")
            }

            kinship <- kinship[
                unique(anno$subject), unique(anno$subject)]
        }

        tu.lambda <- get_eigen(anno, kinship)
        output <- nrow(candidate) %>%
            seq_len %>%
            map(run_map_qtl, candidate=candidate, fn=fn,
                geno=geno, pheno=pheno, anno=anno, covar=covar,
                tu.lambda=tu.lambda,
                type=type, seed=seed) %>%
            bind_rows %>%
            `rownames<-`(NULL)
    } else {
        if (!is.null(kinship)) {
            warning("ignoring `kinship`, as `ranef` is 'FALSE'")
        }
        output <- nrow(candidate) %>%
            seq_len %>%
            map(run_map_qtl, candidate=candidate, fn=fn,
                geno=geno, pheno=pheno, anno=anno, covar=covar,
                type=type, seed=seed) %>%
            bind_rows %>%
            `rownames<-`(NULL)
    }

    apply_fdr <- function(d, fdr) {
        if (is.null(fdr)) return(d)
        if (!"p_value" %in% colnames(d)) return(d)
        d$q_value <- stats::p.adjust(d$p_value, method="fdr")
        d <- d[!is.na(d$q_value) & d$q_value <= fdr, , drop=FALSE]
        rownames(d) <- NULL
        d
    }

    output <- apply_fdr(output, fdr)

    if (excluded) {
        if (!filter.geno) {
            stop("`excluded` can only be `TRUE` when `filter.geno` is `TRUE`")
        }
        excluded <- filtered$excluded
        output <- list(fit=output, excluded=excluded)
    }

    output

}

#' Map marginal or response QTLs for a feature-SNP pair
#'
#' This function takes as input individual phenotype and genotype data
#' for a given feature-SNP pair and tests for marginal and response QTLs.
#'
#' @export
#' @importFrom stats qnorm pchisq sd
#' @param feat.id A character string specifying the feature name.
#' @param snp.id A character string specifying the SNP name.
#' @param input A list containing phenotype, genotype, treatment, and
#'     subject. The elements must be named "y", "g", "t", and
#'     "subject".
#' @param fn A character string specifying the function. This must be
#'     one of "nonlinear" and "linear", corresponding to nonlinear and
#'     linear models, respectively.
#' @param scale A logical scalar as to whether to standardize the
#'     phenotype, \code{y}. This can be set to \code{TRUE} only if
#'     \code{fn} is set to \code{"linear"}.
#' @param rint A logical scalar as to whether the phenotypes need to be
#'     RINT-transformed. This cannot be set to \code{TRUE} if \code{fn}
#'     is set to \code{"nonlinear"}.
#' @param tu.lambda A list obtained from
#'     \code{\link{get_eigen}}. It must contain the transposed
#'     eigenvector matrix and the eigenvalues of the covariance
#'     matrix. The element names must be "tU" and "lambda".
#' @param type A character string specifying the type of QTLs. This
#'     must be one of \code{"marginal"} and \code{"response"}.
#' @param seed An integer specifying a seed for RNG.
#' @param epsilon A small positive scalar. The upper bound of the
#'     \eqn{h^2} search interval is set to \eqn{1 - \epsilon}. Default
#'     is \code{1e-2}. Ignored when \code{ranef} is \code{FALSE}.
#' @param control A list of control parameters passed to
#'     \code{\link[stats]{optim}} for the inner BFGS step. If
#'     \code{NULL}, the defaults are used. Ignored when
#'     \code{gls = TRUE}.
#' @param gls A logical scalar. If \code{TRUE}, use the closed-form
#'     generalized least squares (GLS) estimator instead of BFGS.
#'     Only valid when \code{fn} is \code{"linear"}. Default is
#'     \code{FALSE}.
#'
#' @return A list containing some of the following elements:
#' \itemize{
#' \item{\code{feat_id} - A character string specifying the feature ID.}
#' \item{\code{snp_id} - A character string specifying the SNP ID.}
#' \item{\code{beta0} - A scalar value representing the intercept
#' estimate.}
#' \item{\code{beta_g} - A scalar value representing the genotype
#' effect estimate based on the marginal model (when
#' \code{type = "marginal"}) or the interaction model (when
#' \code{type = "response"}).}
#' \item{\code{beta_t} - A scalar value representing the treatment
#' effect estimate based on the interaction model.}
#' \item{\code{beta_gxt} - A scalar value representing the interaction
#' effect estimate based on the interaction model.}
#' \item{\code{sigma} - A scalar value representing the estimate of
#' the residual error standard deviation.}
#' \item{\code{sigma_u} - A scalar value representing the estimate of
#' the random effect standard deviation.}
#' \item{\code{p_value} - A scalar value representing the nominal
#' p-value from LRT.}
#' \item{\code{return_code} - An integer: 0 if all optimizations
#' converged, \code{NA} otherwise.}
#' }
#'
#' @examples
#' input <- list(
#'     feat.id = "f1",
#'     snp.id = "g1",
#'     y = rnorm(20),
#'     g = rep(c(0, 1, 2), length.out = 20),
#'     t = rep(c(0, 1), 10),
#'     subject = rep(paste0("d", 1:10), each = 2)
#' )
#' result <- map_qtl_for_each(input = input, fn = "linear")
map_qtl_for_each <- function(input,
                             feat.id=NULL,
                             snp.id=NULL,
                             fn,
                             scale=FALSE,
                             rint=FALSE,
                             tu.lambda=NULL,
                             type="response",
                             seed=1,
                             epsilon=1e-2,
                             control=NULL,
                             gls=FALSE) {

    if (is.null(feat.id)) {
        feat.id <- input$feat.id
    }

    if (is.null(snp.id)) {
        snp.id <- input$snp.id
    }

    g <- input$g
    t <- input$t
    y <- input$y
    n <- length(y)
    if (!(is.numeric(g) & is.numeric(t) & is.numeric(y))) {
        stop("`input` has an incorrect format")
    }

    if (!(fn %in% c("linear", "nonlinear"))) {
        stop("`fn` must be one of 'linear' and 'nonlinear'")
    }

    if (!(type %in% c("marginal", "response"))) {
        stop("`type` must be one of 'marginal' and 'response'")
    }

    if (scale) {
        if (fn != "linear") {
            stop("`scale` can be `TRUE` only if `fn` is 'linear'")
        }
        y <- (y - mean(y)) / sd(y)
        input$y <- y
    }

    if (rint) {
        if (fn == "nonlinear") {
            stop("`rint` cannot be `TRUE` if `fn` is 'nonlinear'")
        }
        y <- qnorm((rank(y) - (1/2)) / length(y))
        input$y <- y
    }

    if (is.null(tu.lambda)) {
        ranef <- FALSE
    } else if (is.list(tu.lambda)) {
        ranef <- TRUE
    } else {
        stop("`tu.lambda` has an incorrect format")
    }

    model0 <- model1 <- model2 <- 0
    if (type == "marginal") {
        model0 <- model1 <- 1
    } else if (type == "response") {
        model1 <- model2 <- 1
    }

    if (model0 == 1) {
        m.vec <- c(0, 1, 0)
        fit0 <- get_mle(input=input, fn.gp=fn, m=m.vec,
                        tu.lambda=tu.lambda, epsilon=epsilon,
                        control=control, gls=gls)
    }

    if (model1 == 1) {
        m.vec <- c(1, 1, 0)
        fit1 <- get_mle(input=input, fn.gp=fn, m=m.vec,
                        tu.lambda=tu.lambda, epsilon=epsilon,
                        control=control, gls=gls)
    }

    if (model2 == 1) {
        m.vec <- c(1, 1, 1)
        fit2 <- get_mle(input=input, fn.gp=fn, m=m.vec,
                        tu.lambda=tu.lambda, epsilon=epsilon,
                        control=control, gls=gls)
    }

    get_pval <- function(fit.null, fit.full) {
        stat <- 2 * (fit.null$value - fit.full$value)
        pval <- pchisq(q=stat, df=1, lower.tail=FALSE)
        pval
    }

    get_sigma_estimates <- function(fit, ranef) {
        s2 <- exp(fit$par["ln_s2"])
        h2 <- fit$h2
        sigma <- if (ranef) sqrt((1 - h2) * s2) else sqrt(s2)
        sigma.u <- if (ranef) sqrt(h2 * s2) else NA
        list(sigma=unname(sigma), sigma_u=unname(sigma.u))
    }

    is_valid_fit <- function(fit) {
        is.list(fit) &&
            !inherits(fit, "try-error") &&
            !is.null(fit$par) &&
            !is.null(fit$value) &&
            !is.null(fit$convergence)
    }

    if (type == "marginal") {

        output <- list(
            feat_id=feat.id,
            snp_id=snp.id,
            beta0=NA,
            beta_g=NA,
            sigma=NA,
            sigma_u=NA,
            p_value=NA,
            return_code=NA
        )

        if (is_valid_fit(fit0) & is_valid_fit(fit1)) {
            est <- get_sigma_estimates(fit1, ranef)
            output$beta0 <- unname(fit1$par["b0"])
            output$beta_g <- unname(fit1$par["b1"])
            output$sigma <- est$sigma
            output$sigma_u <- est$sigma_u
            output$p_value <- get_pval(fit0, fit1)
            if (fit0$convergence == 0 & fit1$convergence == 0) {
                output$return_code <- 0L
            }
        }

    } else if (type == "response") {

        output <- list(
            feat_id=feat.id,
            snp_id=snp.id,
            beta0=NA,
            beta_g=NA,
            beta_t=NA,
            beta_gxt=NA,
            sigma=NA,
            sigma_u=NA,
            p_value=NA,
            return_code=NA
        )

        if (is_valid_fit(fit1) & is_valid_fit(fit2)) {
            est <- get_sigma_estimates(fit2, ranef)
            output$beta0 <- unname(fit2$par["b0"])
            output$beta_g <- unname(fit2$par["b1"])
            output$beta_t <- unname(fit2$par["b2"])
            output$beta_gxt <- unname(fit2$par["b3"])
            output$sigma <- est$sigma
            output$sigma_u <- est$sigma_u
            output$p_value <- get_pval(fit1, fit2)
            if (fit1$convergence == 0 & fit2$convergence == 0) {
                output$return_code <- 0L
            }
        }

    }

    output

}

#' Get MLE estimates
#'
#' This function optimizes the log-likelihood over the fixed-effect
#' coefficients and the total variance \eqn{s^2 = \sigma_u^2 + \sigma^2}.
#' When \code{tu.lambda} is supplied (random-effects model), heritability
#' \eqn{h^2 = \sigma_u^2 / s^2} is profiled out via either a grid
#' search or Brent's method.
#'
#' @export
#' @importFrom magrittr "%>%"
#' @importFrom stats lm coef optim dnorm dgamma optimize setNames sigma
#' @param input A list containing phenotype, genotype, treatment, and
#'     subject. The elements must be named "y", "g", "t", and
#'     "subject".
#' @param fn.gp A character string specifying the function to model
#'     the relationship between the genotype and phenotype. This must
#'     be one of "nonlinear" and "linear".
#' @param m A vector of binary indicator variables specifying the
#'     exclusion (0) and inclusion (1) of the genotype, treatment, and
#'     interaction terms.
#' @param tu.lambda A list containing the transposed eigenvector
#'     matrix and the eigenvalues of the covariance matrix (from
#'     \code{\link{get_eigen}}). \code{NULL} means no random effect.
#' @param h2.method A character string. \code{"brent"} (default) uses
#'     Brent's method via \code{optimize} to find the MLE of \eqn{h^2};
#'     \code{"grid"} evaluates the profile likelihood on a fixed grid.
#'     Ignored when \code{tu.lambda} is \code{NULL}.
#' @param n.grid An integer specifying the number of grid points in
#'     \eqn{[0, 1-\epsilon]} when \code{h2.method = "grid"}. Default is 100.
#' @param return.profile A logical scalar. If \code{TRUE} \emph{and}
#'     \code{h2.method = "grid"}, attach the grid data frame
#'     (\code{h2}, \code{log_lik}) to the returned object as
#'     \code{$profile}. Default \code{FALSE}.
#' @param epsilon A small positive scalar. The upper bound of the
#'     \eqn{h^2} search interval is set to \eqn{1 - \epsilon}. Default
#'     is \code{1e-2}. Ignored when \code{tu.lambda} is \code{NULL}.
#' @param control A list of control parameters passed to
#'     \code{\link[stats]{optim}} for the inner BFGS step. If
#'     \code{NULL}, the defaults are used. Ignored when
#'     \code{gls = TRUE}.
#' @param gls A logical scalar. If \code{TRUE}, use the closed-form
#'     generalized least squares (GLS) estimator instead of BFGS to
#'     solve for the fixed effects and variance. Only valid when
#'     \code{fn.gp = "linear"}. Default is \code{FALSE}.
#'
#' @return A list with the elements returned by \code{optim} plus:
#' \itemize{
#'   \item{\code{h2} - MLE of heritability (only when \code{ranef=TRUE}).}
#'   \item{\code{profile} - Data frame of \code{h2} and \code{log_lik}
#'     values (only when \code{h2.method="grid"} and
#'     \code{return.profile=TRUE}).}
#' }
get_mle <- function(input, fn.gp, m,
                    tu.lambda=NULL,
                    h2.method=c("brent", "grid"),
                    n.grid=100,
                    return.profile=FALSE,
                    epsilon=1e-2,
                    control=NULL,
                    gls=FALSE) {

    h2.method <- match.arg(h2.method)

    if (gls && fn.gp != "linear") {
        stop("`gls` can only be `TRUE` when `fn.gp` is 'linear'")
    }

    ranef <- !is.null(tu.lambda)
    if (!is.null(tu.lambda) && !is.list(tu.lambda)) {
        stop("`tu.lambda` has an incorrect format")
    }

    ini <- get_ini(input=input, m=m, ranef=ranef)

    # optimize (betas, ln_s2) via BFGS for fixed h2
    inner_opt <- function(h2) {
        warn0 <- getOption("warn")
        options(warn=-1)
        on.exit(options(warn=warn0))
        args <- list(
            par=ini,
            fn=get_obj,
            input=input,
            m=m,
            fn.gp=fn.gp,
            ranef=ranef,
            tu.lambda=tu.lambda,
            h2=h2,
            method="BFGS",
            hessian=FALSE)
        if (!is.null(control)) args$control <- control
        res <- try(do.call(optim, args), silent=TRUE)
        return(res)
    }

    # solve analytically for (betas, s2) for fixed h2 (GLS/WLS)
    inner_gls <- function(h2) {
        y <- input$y
        g <- input$g
        t <- input$t
        n <- length(y)

        cols <- list(b0=rep(1, n))
        if (m[1] == 1) cols$b1 <- g
        if (m[2] == 1) cols$b2 <- t
        if (m[3] == 1) cols$b3 <- g * t
        X <- do.call(cbind, cols)

        if (!ranef) {
            beta.hat <- tryCatch(
                solve(crossprod(X), crossprod(X, y)),
                error=function(e) NULL)
            if (is.null(beta.hat)) {
                return(structure(list(), class="try-error"))
            }
            resid <- y - X %*% beta.hat
            s2.hat <- sum(resid^2) / n
        } else {
            tU <- tu.lambda$tU
            lambda <- tu.lambda$lambda
            denom <- h2 * lambda + (1 - h2)
            if (any(denom <= 0)) {
                return(structure(list(), class="try-error"))
            }
            w.vec <- 1 / denom
            Ry <- as.vector(tU %*% y)
            RX <- tU %*% X
            sRX <- sqrt(w.vec) * RX
            qr.fit <- qr(sRX)
            if (qr.fit$rank < ncol(RX)) {
                return(structure(list(), class="try-error"))
            }
            beta.hat <- qr.coef(qr.fit, sqrt(w.vec) * Ry)
            resid.R <- Ry - RX %*% beta.hat
            s2.hat <- sum(w.vec * resid.R^2) / n
        }

        par <- c(as.vector(beta.hat), log(s2.hat))
        names(par) <- c(names(cols), "ln_s2")

        nll <- get_obj(param=par, input=input, m=m, fn.gp=fn.gp,
                       ranef=ranef, tu.lambda=tu.lambda, h2=h2)

        list(par=par, value=nll, convergence=0L, message=NULL)
    }

    solve_inner <- if (gls) inner_gls else inner_opt

    if (!ranef) {
        res <- solve_inner(h2=NULL)
        if (!inherits(res, "try-error")) res$h2 <- NA_real_
        return(res)
    }

    # compute profile log-likelihood as a scalar function of h2
    h2.upper <- 1 - epsilon
    profile_ll <- function(h2) {
        r <- solve_inner(h2)
        if (inherits(r, "try-error")) return(-Inf)
        return(-r$value)  # minimize nll; negate to get ll
    }

    if (h2.method == "grid") {

        h2.grid <- seq(0, h2.upper, length.out=n.grid)
        ll.grid <- vapply(h2.grid, profile_ll, numeric(1))
        h2.hat <- h2.grid[which.max(ll.grid)]

    } else {

        opt.h2 <- optimize(
            f=profile_ll,
            interval=c(0, h2.upper),
            maximum=TRUE)
        h2.hat <- opt.h2$maximum

    }

    res <- solve_inner(h2.hat)
    if (inherits(res, "try-error")) {
        return(NA)
    }
    res$h2 <- h2.hat

    if (h2.method == "grid" && return.profile) {
        res$profile <- data.frame(h2=h2.grid, log_lik=ll.grid)
    }

    return(res)
}

#' Plot the profile log-likelihood over h2
#'
#' This function takes the output of \code{\link{get_mle}} called with
#' \code{h2.method = "grid"} and \code{return.profile = TRUE} and
#' draws a line plot of the profile log-likelihood against \eqn{h^2}.
#'
#' @export
#' @importFrom graphics abline
#' @param fit A list returned by \code{\link{get_mle}} that contains a
#'     \code{$profile} element (i.e., fitted with
#'     \code{h2.method = "grid"} and \code{return.profile = TRUE}).
#' @param xlab A character string specifying the x-axis label.
#'     Default \code{"h2"}.
#' @param ylab A character string specifying the y-axis label.
#'     Default \code{"Profile log-likelihood"}.
#' @param main A character string specifying the plot title.
#'     Default \code{""}.
#' @param mark.max A logical scalar. If \code{TRUE} (default), draw a
#'     vertical dashed line at the \eqn{h^2} value that maximizes the
#'     profile log-likelihood.
#' @param ... Further graphical parameters passed to \code{plot}.
#'
#' @return The profile data frame (invisibly).
plot_profile_h2 <- function(fit,
                             xlab="h2",
                             ylab="Profile log-likelihood",
                             main="",
                             mark.max=TRUE,
                             ...) {

    if (is.null(fit$profile)) {
        stop(paste(
            "`fit` does not contain a `$profile` element.",
            "Re-run get_mle() with h2.method='grid' and",
            "return.profile=TRUE."))
    }

    prof <- fit$profile
    plot(prof$h2, prof$log_lik,
         type="l",
         xlab=xlab, ylab=ylab, main=main, ...)

    if (mark.max) {
        h2.max <- prof$h2[which.max(prof$log_lik)]
        abline(v=h2.max, lty=2, col="grey40")
    }

    invisible(prof)
}
