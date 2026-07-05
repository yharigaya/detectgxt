make_candidate <- function(feat.ids, snp.ids) {
    data.frame(feat_id=feat.ids, snp_id=snp.ids, stringsAsFactors=FALSE)
}

make_geno_df <- function(snp.ids, geno.list) {
    n.sub <- length(geno.list[[1]])
    mat <- do.call(rbind, geno.list)
    df <- data.frame(snp_id=snp.ids, mat, row.names=NULL, stringsAsFactors=FALSE)
    colnames(df) <- c("snp_id", paste0("d", seq_len(n.sub)))
    df
}

# test filter_geno()

test_that("filter_geno returns a named list of two data frames", {
    geno <- make_geno_df("g1", list(c(0, 1, 2)))
    candidate <- make_candidate("f1", "g1")
    result <- filter_geno(candidate, geno)
    expect_type(result, "list")
    expect_named(result, c("included", "excluded"))
    expect_s3_class(result$included, "data.frame")
    expect_s3_class(result$excluded, "data.frame")
})

test_that("filter_geno includes SNP when all three genotype levels present", {
    geno <- make_geno_df("g1", list(c(0, 0, 1, 2, 2)))
    candidate <- make_candidate("f1", "g1")
    result <- filter_geno(candidate, geno)
    expect_equal(nrow(result$included), 1)
    expect_equal(nrow(result$excluded), 0)
})

test_that("filter_geno includes SNP when only 0s and 2s present (no heterozygotes)", {
    geno <- make_geno_df("g1", list(c(0, 0, 2, 2)))
    candidate <- make_candidate("f1", "g1")
    result <- filter_geno(candidate, geno)
    expect_equal(nrow(result$included), 1)
    expect_equal(nrow(result$excluded), 0)
})

test_that("filter_geno excludes SNP missing homozygous reference (no 0s)", {
    geno <- make_geno_df("g1", list(c(1, 1, 2, 2, 2)))
    candidate <- make_candidate("f1", "g1")
    result <- filter_geno(candidate, geno)
    expect_equal(nrow(result$included), 0)
    expect_equal(nrow(result$excluded), 1)
})

test_that("filter_geno excludes SNP missing homozygous alternate (no 2s)", {
    geno <- make_geno_df("g1", list(c(0, 0, 0, 1, 1)))
    candidate <- make_candidate("f1", "g1")
    result <- filter_geno(candidate, geno)
    expect_equal(nrow(result$included), 0)
    expect_equal(nrow(result$excluded), 1)
})

test_that("filter_geno warns and treats pair as excluded when SNP absent from geno", {
    geno <- make_geno_df("g_other", list(c(0, 1, 2)))
    candidate <- make_candidate("f1", "g_missing")
    expect_warning(
        result <- filter_geno(candidate, geno),
        "g_missing not found"
    )
    expect_equal(nrow(result$included), 0)
    expect_equal(nrow(result$excluded), 1)
})

test_that("filter_geno correctly separates a mix of passing and failing SNPs", {
    geno <- make_geno_df(
        c("g1", "g2", "g3"),
        list(
            c(0, 0, 1, 2, 2),
            c(1, 1, 2, 2, 2),
            c(0, 0, 0, 1, 1)
        )
    )
    candidate <- make_candidate(
        c("f1", "f2", "f3"),
        c("g1", "g2", "g3")
    )
    result <- filter_geno(candidate, geno)
    expect_equal(nrow(result$included), 1)
    expect_equal(nrow(result$excluded), 2)
    expect_equal(result$included$snp_id, "g1")
    expect_true(all(c("g2", "g3") %in% result$excluded$snp_id))
})

test_that("filter_geno with anno enforces homozygote support in each condition", {
    anno.cond <- data.frame(
        sample=c("s1", "s2", "s3", "s4"),
        subject=c("d1", "d2", "d3", "d4"),
        condition=c(0, 0, 1, 1)
    )
    candidate <- make_candidate("f1", "g1")
    geno <- data.frame(
        snp_id="g1", d1=0, d2=0, d3=2, d4=2,
        stringsAsFactors=FALSE
    )

    result.overall <- filter_geno(candidate, geno)
    expect_equal(nrow(result.overall$included), 1)

    result.by.cond <- filter_geno(candidate, geno, anno=anno.cond)
    expect_equal(nrow(result.by.cond$included), 0)
    expect_equal(nrow(result.by.cond$excluded), 1)
})

test_that("filter_geno with anno errors when required columns are missing", {
    anno.bad <- data.frame(sample="s1", subject="d1")
    candidate <- make_candidate("f1", "g1")
    geno <- make_geno_df("g1", list(c(0, 1, 2)))
    expect_error(
        filter_geno(candidate, geno, anno=anno.bad),
        "`anno` must contain 'sample', 'subject', and 'condition'"
    )
})

# shared fixtures

n.sample <- 40
n.sub <- 20

anno.ut <- data.frame(
    sample=paste0("s", seq_len(n.sample)),
    subject=paste0("d", rep(seq_len(n.sub), each=2)),
    condition=rep(c(0, 1), times=n.sub)
)

n.feat <- 10
count.ut <- make_count_data(anno=anno.ut, n.feat=n.feat, n.sample=n.sample)

# test preprocess_pheno()

test_that("preprocess_pheno returns a data frame with correct dimensions", {
    result <- preprocess_pheno(count=count.ut, anno=anno.ut)
    expect_s3_class(result, "data.frame")
    expect_equal(nrow(result), n.feat)
    expect_equal(ncol(result), n.sample + 1)
})

test_that("preprocess_pheno first column is named feat_id", {
    result <- preprocess_pheno(count=count.ut, anno=anno.ut)
    expect_equal(colnames(result)[1], "feat_id")
})

test_that("preprocess_pheno output contains no NAs", {
    result <- preprocess_pheno(count=count.ut, anno=anno.ut)
    expect_false(anyNA(result))
})

test_that("preprocess_pheno with list=TRUE returns pheno data frame and prcomp object", {
    result <- preprocess_pheno(count=count.ut, anno=anno.ut, list=TRUE)
    expect_type(result, "list")
    expect_named(result, c("pheno", "pca"))
    expect_s3_class(result$pheno, "data.frame")
    expect_s3_class(result$pca, "prcomp")
})

test_that("preprocess_pheno with GxT count data has correct dims, no NAs, and pca$x has n_sample rows", {
    set.seed(8)
    count.gxt <- make_count_data(
        anno=anno.ut, n.feat=n.feat, n.sample=n.sample,
        sd=c(0.5, 0.5, 0.3))
    result <- preprocess_pheno(count=count.gxt, anno=anno.ut, list=TRUE)
    expect_s3_class(result$pheno, "data.frame")
    expect_equal(nrow(result$pheno), n.feat)
    expect_equal(ncol(result$pheno), n.sample + 1)
    expect_false(anyNA(result$pheno))
    expect_equal(nrow(result$pca$x), n.sample)
})

test_that("preprocess_pheno aligns samples by count column names", {
    set.seed(9)
    count.perm <- count.ut[, rev(seq_len(ncol(count.ut)))]

    result <- preprocess_pheno(count=count.ut, anno=anno.ut, num.pc=1)
    result.perm <- preprocess_pheno(count=count.perm, anno=anno.ut, num.pc=1)
    result.perm <- result.perm[, c("feat_id", colnames(count.ut))]

    expect_equal(result, result.perm, tolerance=1e-8)
})

test_that("preprocess_pheno works with custom sample IDs from make_count_data", {
    anno.custom <- anno.ut
    anno.custom$sample <- paste0("x", seq_len(n.sample))
    count.custom <- make_count_data(
        anno=anno.custom, n.feat=5, n.sample=n.sample)

    result <- preprocess_pheno(
        count=count.custom, anno=anno.custom, num.pc=2)

    expect_equal(colnames(count.custom), anno.custom$sample)
    expect_equal(colnames(result)[-1], anno.custom$sample)
})

test_that("preprocess_pheno supports num.pc=0", {
    result <- preprocess_pheno(count=count.ut, anno=anno.ut, num.pc=0)
    expect_s3_class(result, "data.frame")
    expect_equal(nrow(result), n.feat)
    expect_equal(ncol(result), n.sample + 1)
    expect_false(anyNA(result))
})

# test make_count_data()

test_that("make_count_data returns an integer matrix of correct dimensions", {
    set.seed(1)
    result <- make_count_data(anno=anno.ut, n.feat=5, n.sample=n.sample)
    expect_true(is.matrix(result))
    expect_equal(typeof(result), "integer")
    expect_equal(nrow(result), 5)
    expect_equal(ncol(result), n.sample)
})

test_that("make_count_data values are non-negative", {
    set.seed(2)
    result <- make_count_data(anno=anno.ut, n.feat=5, n.sample=n.sample)
    expect_true(all(result >= 0))
})

test_that("make_count_data row and column names use feature and sample IDs", {
    set.seed(3)
    n.f <- 4
    result <- make_count_data(anno=anno.ut, n.feat=n.f, n.sample=n.sample)
    expect_equal(rownames(result), paste0("f", seq_len(n.f)))
    expect_equal(colnames(result), anno.ut$sample)
})

test_that("make_count_data with sd=c(0.5, 0.5, 0.3) returns non-negative integers", {
    set.seed(4)
    result <- make_count_data(
        anno=anno.ut, n.feat=5, n.sample=n.sample,
        sd=c(0.5, 0.5, 0.3))
    expect_true(is.matrix(result))
    expect_equal(typeof(result), "integer")
    expect_true(all(result >= 0))
    expect_equal(nrow(result), 5)
    expect_equal(ncol(result), n.sample)
})

test_that("make_count_data with user-supplied geno matrix returns correct output", {
    set.seed(5)
    n.f <- 3
    sub.names <- unique(anno.ut$subject)
    n.s <- length(sub.names)
    geno.mat <- matrix(
        rep(c(0L, 1L, 2L), length.out=n.f * n.s),
        nrow=n.f, ncol=n.s,
        dimnames=list(NULL, sub.names))
    result <- make_count_data(
        anno=anno.ut, n.feat=n.f, n.sample=n.sample,
        geno=geno.mat, sd=c(0.5, 0.5, 0.3))
    expect_true(is.matrix(result))
    expect_equal(typeof(result), "integer")
    expect_equal(nrow(result), n.f)
    expect_equal(ncol(result), n.sample)
    expect_true(all(result >= 0))
})

test_that("make_count_data with filter.geno=TRUE returns successfully", {
    set.seed(6)
    result <- make_count_data(
        anno=anno.ut, n.feat=2, n.sample=n.sample,
        filter.geno=TRUE)
    expect_true(is.matrix(result))
    expect_equal(nrow(result), 2)
})

test_that("make_count_data errors when filter.geno cannot be satisfied", {
    anno.small <- data.frame(
        sample=paste0("s", 1:2),
        subject=paste0("d", 1:2),
        condition=c(0, 1)
    )
    expect_error(
        make_count_data(
            anno=anno.small, n.feat=1, n.sample=2, filter.geno=TRUE),
        "`filter.geno` requires at least three subjects"
    )
})

test_that("non-zero sd produces different counts than null model", {
    set.seed(7)
    res.null <- make_count_data(anno=anno.ut, n.feat=5, n.sample=n.sample)
    set.seed(7)
    res.effect <- make_count_data(
        anno=anno.ut, n.feat=5, n.sample=n.sample,
        sd=c(2, 2, 2))
    expect_false(isTRUE(all.equal(res.null, res.effect)))
})

test_that("make_count_data errors when geno is not a matrix", {
    expect_error(
        make_count_data(anno=anno.ut, n.feat=2, n.sample=n.sample,
                        geno=c(0, 1, 2)),
        regexp="`geno` must be a matrix")
})

test_that("make_count_data errors when geno has wrong number of rows", {
    sub.names <- unique(anno.ut$subject)
    geno.bad <- matrix(0L, nrow=5, ncol=length(sub.names),
                       dimnames=list(NULL, sub.names))
    expect_error(
        make_count_data(anno=anno.ut, n.feat=3, n.sample=n.sample,
                        geno=geno.bad),
        regexp="`geno` must have n.feat rows")
})

test_that("make_count_data errors when geno has non-matching column names", {
    geno.bad <- matrix(0L, nrow=2, ncol=n.sub,
                       dimnames=list(NULL, paste0("x", seq_len(n.sub))))
    expect_error(
        make_count_data(anno=anno.ut, n.feat=2, n.sample=n.sample,
                        geno=geno.bad),
        regexp="column names of `geno` must match")
})

test_that("make_count_data errors when both sd and coef are NA for an effect", {
    expect_error(
        make_count_data(anno=anno.ut, n.feat=2, n.sample=n.sample,
                        sd=c(NA, 0, 0), coef=c(NA, NA, NA)),
        regexp="only one of")
})

test_that("make_count_data errors when both sd and coef are specified for an effect", {
    expect_error(
        make_count_data(anno=anno.ut, n.feat=2, n.sample=n.sample,
                        sd=c(0.5, 0, 0), coef=c(0.5, NA, NA)),
        regexp="only one of")
})

# test get_geno()

plink.prefix <- file.path(
    system.file("extdata", "plink", package="detecther"), "geno")

test_that("get_geno returns a data frame", {
    candidate <- make_candidate("f1", "g1")
    result <- get_geno(candidate, chromosomes=1,
                       plink.prefix=plink.prefix, plink.suffix="bed")
    expect_s3_class(result, "data.frame")
})

test_that("get_geno first column is named snp_id", {
    candidate <- make_candidate("f1", "g1")
    result <- get_geno(candidate, chromosomes=1,
                       plink.prefix=plink.prefix, plink.suffix="bed")
    expect_equal(colnames(result)[1], "snp_id")
})

test_that("get_geno remaining columns are sample IDs", {
    candidate <- make_candidate("f1", "g1")
    result <- get_geno(candidate, chromosomes=1,
                       plink.prefix=plink.prefix, plink.suffix="bed")
    expect_equal(colnames(result)[-1], paste0("d", 1:6))
})

test_that("get_geno returns one row per requested SNP", {
    candidate <- make_candidate(c("f1", "f2"), c("g1", "g3"))
    result <- get_geno(candidate, chromosomes=1,
                       plink.prefix=plink.prefix, plink.suffix="bed")
    expect_equal(nrow(result), 2)
    expect_setequal(result$snp_id, c("g1", "g3"))
})

test_that("get_geno returns correct dosage values (0, 1, 2)", {
    candidate <- make_candidate("f1", "g1")
    result <- get_geno(candidate, chromosomes=1,
                       plink.prefix=plink.prefix, plink.suffix="bed")
    expect_equal(as.numeric(result[1, -1]), c(0, 0, 1, 1, 2, 2))
})

test_that("get_geno stacks SNPs across chromosomes", {
    candidate <- make_candidate(paste0("f", 1:5), c("g1", "g2", "g3", "g4", "g5"))
    result <- get_geno(candidate, chromosomes=1:2,
                       plink.prefix=plink.prefix, plink.suffix="bed")
    expect_equal(nrow(result), 5)
    expect_setequal(result$snp_id, c("g1", "g2", "g3", "g4", "g5"))
})
