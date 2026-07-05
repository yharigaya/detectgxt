get_input <- function(i, dat, anno) {
    candidate <- dat$candidate
    geno <- dat$geno
    pheno <- dat$pheno
    feat.id <- candidate[i, 1]
    snp.id <- candidate[i, 2]
    if (!(snp.id %in% geno[, 1])) {
        warning(paste(snp.id, "not found in the genotype matrix"))
        return(NULL)
    }
    g <- geno[geno[, 1] == snp.id, -1] %>% as.numeric
    names(g) <- colnames(geno)[-1]
    if (!(feat.id %in% pheno[, 1])) {
        warning(paste(feat.id, "not found in the phenotype matrix"))
        return(NULL)
    }
    y <- pheno[pheno[, 1] == feat.id, -1] %>% as.numeric
    names(y) <- colnames(pheno)[-1]
    d <- anno %>%
        `[`(, c("sample", "subject", "condition")) %>%
        merge(g, by.x="subject", by.y=0) %>%
        `colnames<-`(
            c("subject", "sample", "condition", "geno")) %>%
        merge(y, by.x="sample", by.y=0) %>%
        `colnames<-`(
            c("sample", "subject", "condition", "geno", "pheno"))
    d <- d[match(anno$sample, d$sample), ]
    d
}

call_lm <- function(input) {

    colnames(input) <- c("sample", "subject", "t", "g", "y")
    fit <- lm(y ~ g * t, data=input)
    beta <- stats::coef(fit)
    res <- beta %>%
        `names<-`(c("b0", "b1", "b2", "b3"))
    res

}

call_lm_additive <- function(input) {

    colnames(input) <- c("sample", "subject", "t", "g", "y")
    fit <- lm(y ~ g + t, data=input)
    beta <- stats::coef(fit)
    res <- beta %>%
        `names<-`(c("b0", "b1", "b2"))
    res

}

call_lme <- function(input) {

    colnames(input) <- c("sample", "subject", "t", "g", "y")
    fit <- nlme::lme(y ~ g * t, random=~1|subject,
                     data=input, method="ML")
    sigma.u <- nlme::VarCorr(fit)[1, 2] %>% as.numeric
    res <- c(sigma.u) %>%
        `names<-`("sigma_u")
    res

}

call_nlme <- function(input) {

    colnames(input) <- c("sample", "subject", "t", "g", "y")
    fit.lm <- lm(y ~ g * t, data=input)
    b.hat <- stats::coef(fit.lm)
    fit <- nlme::nlme(
        y ~ get_mu(g, t, b0, b1, b2, b3),
        fixed=b0 + b1 + b2 + b3 ~ 1,
        random=list(subject=b0 ~ 1),
        data=input,
        start=c(b0=b.hat[1], b1=b.hat[2], b2=b.hat[3], b3=b.hat[4]))

    beta <- summary(fit)$coefficients$fixed
    sigma <- sigma(fit)
    sigma.u <- nlme::VarCorr(fit)[1, 2] %>% as.numeric
    pval <- summary(fit)$tTable[, 5]
    loglik <- fit$logLik

    res <- c(beta, sigma, sigma.u, pval, loglik) %>%
        `names<-`(
            c(paste0("b", 0:3), "sigma", "sigma_u",
              paste0("p", 0:3), "loglik"))
    res

}
