tolerance <- 1e-2

n.sample <- 160
n.sub <- 80

anno <- data.frame(
    sample=paste0("s", seq_len(n.sample)),
    subject=paste0("d", rep(seq_len(n.sub), each=2)),
    condition=rep(c(0, 1), times=n.sub))

# check make_sim_data() input

num.min <- c(0, 0, 0, 0, 0, 0, 0, 3)
sd.valid <- c(1, 1, 1)

test_that("make_sim_data errors on invalid output argument", {
    expect_error(
        make_sim_data(num=num.min, anno=anno, fn="nonlinear",
                      sd=sd.valid, output="invalid"),
        "`output` must be one of `data.frame` and `list`"
    )
})

test_that("make_sim_data errors when sigma.u is supplied without ranef=TRUE", {
    expect_error(
        make_sim_data(num=num.min, anno=anno, fn="nonlinear",
                      sd=sd.valid, sigma.u=1),
        "`ranef` must be set to `TRUE` when `sigma.u` is supplied"
    )
})

test_that("make_sim_data errors when kinship is supplied without ranef=TRUE", {
    kin <- diag(n.sub)
    rownames(kin) <- colnames(kin) <- paste0("d", seq_len(n.sub))
    expect_error(
        make_sim_data(num=num.min, anno=anno, fn="nonlinear",
                      sd=sd.valid, kinship=kin),
        "`ranef` must be set to `TRUE` when `kinship` is supplied"
    )
})

test_that("make_sim_data errors when neither coef nor sd is specified", {
    expect_error(
        make_sim_data(num=num.min, anno=anno, fn="nonlinear"),
        "at least one of the `coef` and `sd` needs to be specified"
    )
})

test_that("make_sim_data errors when coef has wrong length", {
    expect_error(
        make_sim_data(num=num.min, anno=anno, fn="nonlinear",
                      coef=c(1, 2)),
        "the length of `coef` must be three"
    )
})

test_that("make_sim_data errors when sd has wrong length", {
    expect_error(
        make_sim_data(num=num.min, anno=anno, fn="nonlinear",
                      sd=c(1, 2)),
        "the length of `sd` must be three"
    )
})

test_that("make_sim_data errors when num is not length 8", {
    expect_error(
        make_sim_data(num=c(1, 2), anno=anno, fn="nonlinear", sd=sd.valid),
        "`num` must be a non-negative integer vector of length 8"
    )
})

test_that("make_sim_data errors when both coef and sd given for genotype effect", {
    expect_error(
        make_sim_data(num=num.min, anno=anno, fn="nonlinear",
                      sd=sd.valid, coef=c(1, NA, NA)),
        "genotype effect"
    )
})

test_that("make_sim_data errors when both coef and sd given for treatment effect", {
    expect_error(
        make_sim_data(num=num.min, anno=anno, fn="nonlinear",
                      sd=sd.valid, coef=c(NA, 1, NA)),
        "treatment effect"
    )
})

test_that("make_sim_data errors when both coef and sd given for interaction effect", {
    expect_error(
        make_sim_data(num=num.min, anno=anno, fn="nonlinear",
                      sd=sd.valid, coef=c(NA, NA, 1)),
        "interaction effect"
    )
})

# test filter.geno

num.filter <- c(0, 0, 0, 0, 0, 0, 0, 20)
coef.filter <- c(0.5, 0.3, 0.2)

anno.paired <- data.frame(
    sample=paste0("s", 1:8),
    subject=rep(paste0("d", 1:4), each=2),
    condition=rep(c(0, 1), 4)
)

anno.unpaired <- data.frame(
    sample=paste0("s", 1:8),
    subject=paste0("d", 1:8),
    condition=c(0, 0, 0, 0, 1, 1, 1, 1)
)

anno.partial <- data.frame(
    sample=paste0("s", 1:6),
    subject=c("d1", "d1", "d2", "d2", "d3", "d4"),
    condition=c(0, 1, 0, 1, 0, 1)
)

anno.repeated <- data.frame(
    sample=paste0("s", 1:8),
    subject=c("d1", "d1", "d1", "d2", "d2", "d3", "d3", "d4"),
    condition=c(0, 0, 1, 0, 1, 0, 1, 1)
)

check_filter_geno <- function(dat, anno) {
    geno.mat <- dat$geno
    subjects <- colnames(geno.mat)[-1]
    all(apply(geno.mat[, -1, drop=FALSE], 1, function(g.row) {
        names(g.row) <- subjects
        d.tmp <- merge(anno, g.row, by.x="subject", by.y=0)
        colnames(d.tmp)[ncol(d.tmp)] <- "geno"
        all(tapply(d.tmp$geno, d.tmp$condition,
                   function(g) any(g == 0) & any(g == 2)))
    }))
}

test_that("filter.geno guarantees both homozygotes in each condition (paired)", {
    dat <- make_sim_data(
        num=num.filter, anno=anno.paired, fn="nonlinear",
        coef=coef.filter, filter.geno=TRUE)
    expect_true(check_filter_geno(dat, anno.paired))
})

test_that("filter.geno guarantees both homozygotes in each condition (unpaired)", {
    dat <- make_sim_data(
        num=num.filter, anno=anno.unpaired, fn="nonlinear",
        coef=coef.filter, filter.geno=TRUE)
    expect_true(check_filter_geno(dat, anno.unpaired))
})

test_that("filter.geno guarantees both homozygotes in each condition (partially paired)", {
    dat <- make_sim_data(
        num=num.filter, anno=anno.partial, fn="nonlinear",
        coef=coef.filter, filter.geno=TRUE)
    expect_true(check_filter_geno(dat, anno.partial))
})

test_that("filter.geno guarantees both homozygotes in each condition (repeated measurements)", {
    dat <- make_sim_data(
        num=num.filter, anno=anno.repeated, fn="nonlinear",
        coef=coef.filter, filter.geno=TRUE)
    expect_true(check_filter_geno(dat, anno.repeated))
})

test_that("make_sim_data errors when filter.geno cannot be satisfied", {
    anno.small <- data.frame(
        sample=paste0("s", 1:2),
        subject=paste0("d", 1:2),
        condition=c(0, 1)
    )
    expect_error(
        make_sim_data(
            num=c(0, 0, 0, 0, 0, 0, 0, 1),
            anno=anno.small, fn="nonlinear",
            coef=coef.filter, filter.geno=TRUE),
        "`filter.geno` requires at least two subjects in each condition"
    )
})

# check map_qtl_for_each() input

input.valid <- list(
    feat.id="f1",
    snp.id="g1",
    y=rnorm(n.sample),
    g=rep(c(0, 1, 2), length.out=n.sample),
    t=rep(c(0, 1), n.sample / 2),
    subject=rep(paste0("d", seq_len(n.sub)), each=2)
)

test_that("map_qtl_for_each errors when y is non-numeric", {
    bad <- input.valid
    bad$y <- as.character(bad$y)
    expect_error(
        map_qtl_for_each(input=bad, fn="linear"),
        "`input` has an incorrect format"
    )
})

test_that("map_qtl_for_each errors when g is non-numeric", {
    bad <- input.valid
    bad$g <- as.character(bad$g)
    expect_error(
        map_qtl_for_each(input=bad, fn="linear"),
        "`input` has an incorrect format"
    )
})

test_that("map_qtl_for_each errors when t is non-numeric", {
    bad <- input.valid
    bad$t <- as.character(bad$t)
    expect_error(
        map_qtl_for_each(input=bad, fn="linear"),
        "`input` has an incorrect format"
    )
})

test_that("map_qtl_for_each errors when scale=TRUE with fn='nonlinear'", {
    expect_error(
        map_qtl_for_each(input=input.valid, fn="nonlinear", scale=TRUE),
        "`scale` can be `TRUE` only if `fn` is 'linear'"
    )
})

test_that("map_qtl_for_each errors when rint=TRUE with fn='nonlinear'", {
    expect_error(
        map_qtl_for_each(input=input.valid, fn="nonlinear", rint=TRUE),
        "`rint` cannot be `TRUE` if `fn` is 'nonlinear'"
    )
})

test_that("map_qtl_for_each errors when tu.lambda is not NULL or a list", {
    expect_error(
        map_qtl_for_each(input=input.valid, fn="linear", tu.lambda="bad"),
        "`tu.lambda` has an incorrect format"
    )
})

test_that("map_qtl_for_each errors when fn is invalid", {
    expect_error(
        map_qtl_for_each(input=input.valid, fn="invalid"),
        "`fn` must be one of 'linear' and 'nonlinear'"
    )
})

test_that("map_qtl_for_each errors when type is invalid", {
    expect_error(
        map_qtl_for_each(input=input.valid, fn="linear", type="invalid"),
        "`type` must be one of 'marginal' and 'response'"
    )
})

test_that("map_qtl_for_each with rank-deficient GLS returns failed fit output", {
    input.rank.deficient <- input.valid
    input.rank.deficient$t <- rep(0, n.sample)

    result <- map_qtl_for_each(
        input=input.rank.deficient, fn="linear",
        type="response", gls=TRUE)

    expect_true(is.na(result$p_value))
    expect_true(is.na(result$return_code))
})

# get map_qtl() comparison fixtures

sd <- c(1.5, 2.0, 1.0)
num <- c(0, 0, 0, 0, 0, 0, 0, 100)
sigma.u <- 1

d1 <- make_sim_data(
    num=num, anno=anno, fn="nonlinear", sd=sd, filter.geno=TRUE)

d2 <- make_sim_data(
    num=num, anno=anno, fn="nonlinear", sd=sd,
    ranef=TRUE, sigma.u=sigma.u, filter.geno=TRUE)

d3 <- make_sim_data(
    num=num, anno=anno, fn="linear", sd=sd, filter.geno=TRUE)

d4 <- make_sim_data(
    num=num, anno=anno, fn="linear", sd=sd,
    ranef=TRUE, sigma.u=sigma.u, filter.geno=TRUE)

# test map_qtl() comparisons

test_that("map_qtl works for data generated without random effect", {
    skip_if_not_installed("nlme")

    nlme1 <- seq_len(nrow(d1$candidate)) %>%
        map(get_input, dat=d1, anno=anno) %>%
        map(call_nlme) %>%
        bind_rows %>%
        as.data.frame

    res1 <- map_qtl(
        candidate=d1$candidate, fn="nonlinear",
        geno=d1$geno, pheno=d1$pheno, anno=anno,
        ranef=TRUE, type="response")

    expect_true(
        max(abs(res1$sigma_u - nlme1$sigma_u)) < tolerance)

})

test_that("map_qtl works for data generated with random effect", {
    skip_if_not_installed("nlme")

    nlme2 <- seq_len(nrow(d2$candidate)) %>%
        map(get_input, dat=d2, anno=anno) %>%
        map(call_nlme) %>%
        bind_rows %>%
        as.data.frame

    res2 <- map_qtl(
        candidate=d2$candidate, fn="nonlinear",
        geno=d2$geno, pheno=d2$pheno, anno=anno,
        ranef=TRUE, type="response")

    expect_true(
        max(abs(res2$sigma_u - nlme2$sigma_u)) < tolerance)

})

test_that("map_qtl works for marginal type with linear model", {

    lm3.add <- seq_len(nrow(d3$candidate)) %>%
        map(get_input, dat=d3, anno=anno) %>%
        map(call_lm_additive) %>%
        bind_rows %>%
        as.data.frame

    res3.marginal <- map_qtl(
        candidate=d3$candidate, fn="linear",
        geno=d3$geno, pheno=d3$pheno, anno=anno,
        type="marginal")

    expect_true(max(abs(res3.marginal$beta_g - lm3.add$b1)) < tolerance)

})

test_that("map_qtl works for marginal type with nonlinear model", {

    res1.marginal <- map_qtl(
        candidate=d1$candidate, fn="nonlinear",
        geno=d1$geno, pheno=d1$pheno, anno=anno,
        type="marginal")

    expect_true(all(!is.na(res1.marginal$beta_g)))
    expect_true(all(res1.marginal$return_code == 0))

})

test_that("map_qtl works for linear model without random effect", {

    lm3 <- seq_len(nrow(d3$candidate)) %>%
        map(get_input, dat=d3, anno=anno) %>%
        map(call_lm) %>%
        bind_rows %>%
        as.data.frame

    res3 <- map_qtl(
        candidate=d3$candidate, fn="linear",
        geno=d3$geno, pheno=d3$pheno, anno=anno,
        type="response")

    expect_true(max(abs(res3$beta_g - lm3$b1)) < tolerance)
    expect_true(max(abs(res3$beta_t - lm3$b2)) < tolerance)
    expect_true(max(abs(res3$beta_gxt - lm3$b3)) < tolerance)

})

test_that("map_qtl works for linear model with random effect", {
    skip_if_not_installed("nlme")

    lme4 <- seq_len(nrow(d4$candidate)) %>%
        map(get_input, dat=d4, anno=anno) %>%
        map(call_lme) %>%
        bind_rows %>%
        as.data.frame

    res4 <- map_qtl(
        candidate=d4$candidate, fn="linear",
        geno=d4$geno, pheno=d4$pheno, anno=anno,
        ranef=TRUE, type="response")

    expect_true(max(abs(res4$sigma_u - lme4$sigma_u)) < tolerance)

})

# get p-value fixtures

sd.pval <- c(1.5, 2.0, 1.0)

d.nl.signal <- make_sim_data(
    num=c(0, 0, 0, 0, 0, 0, 0, 20), anno=anno,
    fn="nonlinear", sd=sd.pval, filter.geno=TRUE)

d.nl.null <- make_sim_data(
    num=c(0, 0, 0, 20, 0, 0, 0, 0), anno=anno,
    fn="nonlinear", sd=sd.pval, filter.geno=TRUE)

d.lm.signal <- make_sim_data(
    num=c(0, 0, 0, 0, 0, 0, 0, 20), anno=anno,
    fn="linear", sd=sd.pval, filter.geno=TRUE)

res.nl.signal <- map_qtl(
    candidate=d.nl.signal$candidate, fn="nonlinear",
    geno=d.nl.signal$geno, pheno=d.nl.signal$pheno, anno=anno,
    type="response")

res.nl.null <- map_qtl(
    candidate=d.nl.null$candidate, fn="nonlinear",
    geno=d.nl.null$geno, pheno=d.nl.null$pheno, anno=anno,
    type="response")

res.lm.signal <- map_qtl(
    candidate=d.lm.signal$candidate, fn="linear",
    geno=d.lm.signal$geno, pheno=d.lm.signal$pheno, anno=anno,
    type="response")

res.lm.signal.fdr <- map_qtl(
    candidate=d.lm.signal$candidate, fn="linear",
    geno=d.lm.signal$geno, pheno=d.lm.signal$pheno, anno=anno,
    type="response", fdr=0.2)

res.nl.marginal <- map_qtl(
    candidate=d.nl.signal$candidate, fn="nonlinear",
    geno=d.nl.signal$geno, pheno=d.nl.signal$pheno, anno=anno,
    type="marginal")

# test p-values

test_that("p-values are in [0, 1] for nonlinear response type", {
    expect_true(all(res.nl.signal$p_value >= 0 & res.nl.signal$p_value <= 1))
})

test_that("return_code is 0 for all pairs in nonlinear response type", {
    expect_true(all(res.nl.signal$return_code == 0))
})

test_that("p-values are small when GxT interaction is present (nonlinear)", {
    expect_true(mean(res.nl.signal$p_value < 0.05) > 0.5)
})

test_that("p-values are not systematically small under the null (nonlinear)", {
    expect_true(mean(res.nl.null$p_value) > 0.1)
})

test_that("p-values are in [0, 1] for linear response type", {
    expect_true(all(res.lm.signal$p_value >= 0 & res.lm.signal$p_value <= 1))
})

test_that("return_code is 0 for all pairs in linear response type", {
    expect_true(all(res.lm.signal$return_code == 0))
})

test_that("fdr filtering adds q_value and keeps only rows at threshold", {
    expect_true("q_value" %in% colnames(res.lm.signal.fdr))
    expect_true(all(!is.na(res.lm.signal.fdr$q_value)))
    expect_true(all(res.lm.signal.fdr$q_value <= 0.2))
    expect_true(nrow(res.lm.signal.fdr) <= nrow(res.lm.signal))
})

test_that("p-values are small when GxT interaction is present (linear)", {
    expect_true(mean(res.lm.signal$p_value < 0.05) > 0.5)
})

test_that("map_qtl errors when fdr is outside [0, 1]", {
    expect_error(
        map_qtl(
            candidate=d.lm.signal$candidate[1, , drop=FALSE], fn="linear",
            geno=d.lm.signal$geno, pheno=d.lm.signal$pheno, anno=anno,
            type="response", fdr=1.1),
        "`fdr` must be a numeric scalar between 0 and 1"
    )
})

test_that("p-values are in [0, 1] for marginal type", {
    expect_true(all(res.nl.marginal$p_value >= 0 & res.nl.marginal$p_value <= 1))
})

test_that("return_code is 0 for all pairs in marginal type", {
    expect_true(all(res.nl.marginal$return_code == 0))
})

# get map_qtl() fixtures with excluded / covar

n.sample.ut <- 40
n.sub.ut <- 20

anno.ut <- data.frame(
    sample=paste0("s", seq_len(n.sample.ut)),
    subject=paste0("d", rep(seq_len(n.sub.ut), each=2)),
    condition=rep(c(0, 1), times=n.sub.ut)
)

geno.ut <- data.frame(
    snp_id=c("g_ok", "g_bad"),
    stringsAsFactors=FALSE
)
geno.ut[1, paste0("d", seq_len(n.sub.ut))] <- rep(c(0, 1, 2), length.out=n.sub.ut)
geno.ut[2, paste0("d", seq_len(n.sub.ut))] <- rep(1, n.sub.ut)

set.seed(42)
pheno.ut <- data.frame(
    feat_id=c("f1", "f2"),
    matrix(rnorm(2 * n.sample.ut), nrow=2),
    stringsAsFactors=FALSE
)
colnames(pheno.ut)[-1] <- paste0("s", seq_len(n.sample.ut))

candidate.ut <- data.frame(
    feat_id=c("f1", "f2"),
    snp_id=c("g_ok", "g_bad"),
    stringsAsFactors=FALSE
)

# test map_qtl() with excluded=TRUE

test_that("map_qtl with excluded=TRUE returns a named list", {
    result <- map_qtl(
        candidate=candidate.ut, fn="linear",
        geno=geno.ut, pheno=pheno.ut, anno=anno.ut,
        filter.geno=TRUE, excluded=TRUE, type="response")
    expect_type(result, "list")
    expect_named(result, c("fit", "excluded"))
    expect_s3_class(result$fit, "data.frame")
    expect_s3_class(result$excluded, "data.frame")
})

test_that("map_qtl with excluded=TRUE correctly partitions included and excluded pairs", {
    result <- map_qtl(
        candidate=candidate.ut, fn="linear",
        geno=geno.ut, pheno=pheno.ut, anno=anno.ut,
        filter.geno=TRUE, excluded=TRUE, type="response")
    expect_equal(nrow(result$fit), 1)
    expect_equal(nrow(result$excluded), 1)
    expect_equal(result$fit$snp_id, "g_ok")
    expect_equal(result$excluded$snp_id, "g_bad")
})

test_that("map_qtl with excluded=TRUE reports missing SNPs as excluded", {
    candidate.missing <- data.frame(
        feat_id="f1", snp_id="g_missing",
        stringsAsFactors=FALSE)

    expect_warning(
        result <- map_qtl(
            candidate=candidate.missing, fn="linear",
            geno=geno.ut, pheno=pheno.ut, anno=anno.ut,
            filter.geno=TRUE, excluded=TRUE, type="response"),
        "g_missing not found"
    )

    expect_equal(nrow(result$fit), 0)
    expect_equal(nrow(result$excluded), 1)
    expect_equal(result$excluded$snp_id, "g_missing")
})

# test map_qtl() with covar

covar.ut <- data.frame(
    sample=anno.ut$sample,
    batch=rep(c(0, 1), n.sample.ut / 2)
)

candidate.one <- data.frame(feat_id="f1", snp_id="g_ok", stringsAsFactors=FALSE)

test_that("map_qtl with covar returns a data frame with p_value column", {
    result <- map_qtl(
        candidate=candidate.one, fn="linear",
        geno=geno.ut, pheno=pheno.ut, anno=anno.ut,
        covar=covar.ut, type="response")
    expect_s3_class(result, "data.frame")
    expect_true("p_value" %in% colnames(result))
})

test_that("map_qtl with covar produces a p-value in [0, 1]", {
    result <- map_qtl(
        candidate=candidate.one, fn="linear",
        geno=geno.ut, pheno=pheno.ut, anno=anno.ut,
        covar=covar.ut, type="response")
    expect_true(result$p_value >= 0 & result$p_value <= 1)
})

# covariate correction fixtures

set.seed(42)
n.sample.cv <- 40
n.sub.cv <- 20

anno.cv <- data.frame(
    sample=paste0("s", seq_len(n.sample.cv)),
    subject=paste0("d", rep(seq_len(n.sub.cv), each=2)),
    condition=rep(c(0, 1), times=n.sub.cv)
)

count.cv <- make_count_data(
    anno=anno.cv, n.feat=40, n.sample=n.sample.cv,
    sd=c(1, 0, 0))

pp.cv <- preprocess_pheno(count=count.cv, anno=anno.cv, list=TRUE)
pc1.cv <- pp.cv$pca$x[anno.cv$sample, 1]
pc1.cv <- pc1.cv / sd(pc1.cv)

d.cv <- make_sim_data(
    num=c(0, 0, 0, 0, 0, 0, 0, 20), anno=anno.cv,
    fn="linear", sd=c(1.5, 2.0, 1.0), filter.geno=TRUE)

res.cv.clean <- map_qtl(
    candidate=d.cv$candidate, fn="linear",
    geno=d.cv$geno, pheno=d.cv$pheno, anno=anno.cv,
    type="response")

multiplier <- 5 * mean(res.cv.clean$sigma)
pheno.cv.corrupt <- d.cv$pheno
pheno.cv.corrupt[, -1] <- sweep(pheno.cv.corrupt[, -1], 2, pc1.cv * multiplier, "+")

covar.cv <- data.frame(sample=anno.cv$sample, pc1=pc1.cv)

res.cv.no.covar <- map_qtl(
    candidate=d.cv$candidate, fn="linear",
    geno=d.cv$geno, pheno=pheno.cv.corrupt, anno=anno.cv,
    type="response")

res.cv.with.covar <- map_qtl(
    candidate=d.cv$candidate, fn="linear",
    geno=d.cv$geno, pheno=pheno.cv.corrupt, anno=anno.cv,
    covar=covar.cv, type="response")

# test covariate correction

test_that("covar correction changes results relative to no correction", {
    expect_false(isTRUE(all.equal(res.cv.no.covar$sigma, res.cv.with.covar$sigma)))
})

test_that("covar correction recovers sigma closer to clean data than no correction", {
    err.no.covar <- mean(abs(res.cv.no.covar$sigma - res.cv.clean$sigma))
    err.with.covar <- mean(abs(res.cv.with.covar$sigma - res.cv.clean$sigma))
    expect_true(err.with.covar < err.no.covar)
})

test_that("p-values with covar correction are in [0, 1] and NA-free", {
    expect_true(all(!is.na(res.cv.with.covar$p_value)))
    expect_true(all(res.cv.with.covar$p_value >= 0 &
                    res.cv.with.covar$p_value <= 1))
})

# t=0-only fixtures

anno.t0 <- anno[anno$condition == 0, ]

d.t0.nl <- make_sim_data(
    num=c(0, 0, 0, 0, 0, 0, 0, 20), anno=anno.t0,
    fn="nonlinear", sd=sd.pval, filter.geno=TRUE)

d.t0.lm <- make_sim_data(
    num=c(0, 0, 0, 0, 0, 0, 0, 20), anno=anno.t0,
    fn="linear", sd=sd.pval, filter.geno=TRUE)

res.t0.nl <- map_qtl(
    candidate=d.t0.nl$candidate, fn="nonlinear",
    geno=d.t0.nl$geno, pheno=d.t0.nl$pheno, anno=anno.t0,
    type="marginal")

res.t0.lm <- map_qtl(
    candidate=d.t0.lm$candidate, fn="linear",
    geno=d.t0.lm$geno, pheno=d.t0.lm$pheno, anno=anno.t0,
    type="marginal")

# test t=0-only data

test_that("map_qtl with t=0-only data returns a data frame (nonlinear)", {
    expect_s3_class(res.t0.nl, "data.frame")
    expect_equal(nrow(res.t0.nl), 20)
})

test_that("p-values are in [0, 1] for t=0-only nonlinear marginal", {
    expect_true(all(!is.na(res.t0.nl$p_value)))
    expect_true(all(res.t0.nl$p_value >= 0 & res.t0.nl$p_value <= 1))
})

test_that("return_code is 0 for all pairs in t=0-only nonlinear marginal", {
    expect_true(all(res.t0.nl$return_code == 0))
})

test_that("map_qtl with t=0-only data returns a data frame (linear)", {
    expect_s3_class(res.t0.lm, "data.frame")
    expect_equal(nrow(res.t0.lm), 20)
})

test_that("p-values are in [0, 1] for t=0-only linear marginal", {
    expect_true(all(!is.na(res.t0.lm$p_value)))
    expect_true(all(res.t0.lm$p_value >= 0 & res.t0.lm$p_value <= 1))
})

test_that("return_code is 0 for all pairs in t=0-only linear marginal", {
    expect_true(all(res.t0.lm$return_code == 0))
})

# t=1-only fixtures

anno.t1 <- anno.t0
anno.t1$condition <- 1L

res.t1.nl <- map_qtl(
    candidate=d.t0.nl$candidate, fn="nonlinear",
    geno=d.t0.nl$geno, pheno=d.t0.nl$pheno, anno=anno.t1,
    type="marginal")

res.t1.lm <- map_qtl(
    candidate=d.t0.lm$candidate, fn="linear",
    geno=d.t0.lm$geno, pheno=d.t0.lm$pheno, anno=anno.t1,
    type="marginal")

# test t=1-only data

test_that("p-values are in [0, 1] for t=1-only nonlinear marginal", {
    expect_true(all(!is.na(res.t1.nl$p_value)))
    expect_true(all(res.t1.nl$p_value >= 0 & res.t1.nl$p_value <= 1))
})

test_that("return_code is 0 for all pairs in t=1-only nonlinear marginal", {
    expect_true(all(res.t1.nl$return_code == 0))
})

test_that("p-values are in [0, 1] for t=1-only linear marginal", {
    expect_true(all(!is.na(res.t1.lm$p_value)))
    expect_true(all(res.t1.lm$p_value >= 0 & res.t1.lm$p_value <= 1))
})

test_that("return_code is 0 for all pairs in t=1-only linear marginal", {
    expect_true(all(res.t1.lm$return_code == 0))
})

test_that("recoding t=0 to t=1 does not change p-values (nonlinear)", {
    expect_equal(res.t1.nl$p_value, res.t0.nl$p_value, tolerance=tolerance)
})

test_that("recoding t=0 to t=1 does not change beta_g (nonlinear)", {
    expect_equal(res.t1.nl$beta_g, res.t0.nl$beta_g, tolerance=tolerance)
})

test_that("recoding t=0 to t=1 does not change p-values (linear)", {
    expect_equal(res.t1.lm$p_value, res.t0.lm$p_value, tolerance=tolerance)
})

test_that("recoding t=0 to t=1 does not change beta_g (linear)", {
    expect_equal(res.t1.lm$beta_g, res.t0.lm$beta_g, tolerance=tolerance)
})

# kinship validation fixtures

n.sub.kn <- 5
sub.names.kn <- paste0("d", seq_len(n.sub.kn))
anno.kn <- data.frame(
    sample=paste0("s", seq_len(2 * n.sub.kn)),
    subject=rep(sub.names.kn, each=2),
    condition=rep(c(0, 1), n.sub.kn)
)

K.valid <- diag(n.sub.kn)
rownames(K.valid) <- colnames(K.valid) <- sub.names.kn

K.asym <- K.valid
K.asym[1, 2] <- 0.5  # K[1,2] != K[2,1]

K.npsd <- K.valid
K.npsd[1, 2] <- K.npsd[2, 1] <- -2  # leading 2x2 eigenvalue = -1

# test get_eigen() kinship validation

test_that("get_eigen errors on non-symmetric kinship", {
    expect_error(
        get_eigen(anno=anno.kn, kinship=K.asym),
        "not symmetric"
    )
})

test_that("get_eigen succeeds without warning on non-PSD kinship", {
    expect_no_condition(
        get_eigen(anno=anno.kn, kinship=K.npsd)
    )
})

test_that("get_eigen succeeds and returns tU and lambda for valid kinship", {
    result <- get_eigen(anno=anno.kn, kinship=K.valid)
    expect_type(result, "list")
    expect_named(result, c("tU", "lambda"))
})

# test make_sim_data() kinship validation

test_that("make_sim_data errors on non-symmetric kinship", {
    expect_error(
        make_sim_data(
            num=c(0, 0, 0, 0, 0, 0, 0, 1), anno=anno.kn, fn="nonlinear",
            coef=c(0.5, 0.3, 0.2), ranef=TRUE, kinship=K.asym),
        "not symmetric"
    )
})

test_that("make_sim_data errors on non-PD variance component", {
    expect_error(
        make_sim_data(
            num=c(0, 0, 0, 0, 0, 0, 0, 1), anno=anno.kn, fn="nonlinear",
            coef=c(0.5, 0.3, 0.2), ranef=TRUE, kinship=K.npsd),
        "not positive definite"
    )
})

test_that("make_sim_data succeeds with valid PSD kinship", {
    result <- make_sim_data(
        num=c(0, 0, 0, 0, 0, 0, 0, 1), anno=anno.kn, fn="nonlinear",
        coef=c(0.5, 0.3, 0.2), ranef=TRUE, kinship=K.valid)
    expect_type(result, "list")
    expect_named(result, c("candidate", "geno", "pheno"))
})
