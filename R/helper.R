
# compute negative log-likelihood
#
# parameterization:
#   ranef=FALSE  param contains ln_s2  (s2 = sigma^2, total = residual variance)
#   ranef=TRUE   param contains ln_s2  (s2 = sigma_u^2 + sigma^2)
#                h2 is passed separately and is fixed during inner optimization;
#                sigma_u^2 = h2 * s2,  sigma^2 = (1 - h2) * s2
get_obj <- function(param, input, m,
                    fn.gp, ranef, tu.lambda, h2=NULL) {

    s2 <- exp(param["ln_s2"])
    y <- input$y
    g <- input$g
    t <- input$t
    n <- length(y)
    mu <- get_mean(g=g, t=t, param=param, m=m, fn=fn.gp)

    if (ranef) {
        tU <- tu.lambda$tU
        lambda <- tu.lambda$lambda
        # diagonal of rotated covariance: s2 * (h2*lambda + (1-h2))
        d.vec <- s2 * (h2 * lambda + (1 - h2))
        if (any(d.vec <= 0)) return(Inf)
        R <- tU %*% (y - mu)
        nxh <- (1/2) * sum(R * R / d.vec) +
               (1/2) * sum(log(d.vec))    +
               (1/2) * n * log(2 * pi)
    } else {
        nxh <- -sum(dnorm(x=y, mean=mu, sd=sqrt(s2), log=TRUE))
    }

    unname(nxh)
}

# compute starting values for inner optimization over (betas, ln_s2)
#
# for both ranef=TRUE and ranef=FALSE the free parameters are the same:
# betas determined by m, plus ln_s2.  h2 is not a free parameter here.
get_ini <- function(input, m, ranef) {
    d <- data.frame(y=input$y, g=input$g, t=input$t)
    all.terms <- c("g", "t", "g:t")
    all.bnames <- c("b1", "b2", "b3")

    sel.terms <- all.terms[m == 1]
    sel.bnames <- all.bnames[m == 1]

    # only include terms whose component variables are not constant;
    # a constant predictor (e.g. t=0 throughout) makes the design matrix
    # rank-deficient and yields NA coefficients
    is_variable <- function(term) {
        vars <- strsplit(term, ":")[[1]]
        all(vapply(vars, function(v) stats::var(d[[v]]) > 0, logical(1)))
    }
    active <- vapply(sel.terms, is_variable, logical(1))
    fit.terms <- sel.terms[active]
    fit.bnames <- sel.bnames[active]

    frml <- if (length(fit.terms) == 0) {
        "y ~ 1"
    } else {
        paste(c("y ~ 1", paste("+", fit.terms)), collapse=" ")
    }
    lm.fit <- lm(frml, data=d)

    # ln_s2 = log(sigma^2): same starting value whether ranef or not;
    # s2 here approximates the total variance via the residual sd.
    s2.ini <- sigma(lm.fit)^2

    # initialize all selected betas to 0, then fill in lm estimates
    ini <- setNames(
        c(rep(0, length(sel.bnames) + 1), log(s2.ini)),
        c("b0", sel.bnames, "ln_s2"))
    ini["b0"] <- unname(coef(lm.fit)[1])
    if (length(fit.bnames) > 0)
        ini[fit.bnames] <- unname(coef(lm.fit)[-1])
    ini
}

get_mean <- function(g, t, param, m, fn) {
    beta <- get_param(param, m)
    b0 <- beta[1]; b1 <- beta[2]
    b2 <- beta[3]; b3 <- beta[4]

    if (!(fn %in% c("linear", "nonlinear"))) {
        stop("fn must be either 'linear' or 'nonlinear'")
    }

    if (fn == "linear") {
        output <- b0 + b1 * g + b2 * t + b3 * g * t
    } else if (fn == "nonlinear") {
        tmp <-  exp(b0) * (1 - (g/2)) * (1 - t) +
            exp(b0 + 2 * b1) * (g/2) * (1 - t) +
            exp(b0 + b2) * (1 - (g/2)) * t +
            exp(b0 + 2 * b1 + b2 + 2 * b3) * (g/2) * t
        output <- log(tmp)
    }
    output
}

get_param <- function(param, m) {
    beta <- rep(0, 4)
    beta[1] <- param["b0"]
    if (m[1] == 1) beta[2] <- param["b1"]
    if (m[2] == 1) beta[3] <- param["b2"]
    if (m[3] == 1) beta[4] <- param["b3"]
    beta
}

get_model_mat <- function() {
    model.mat <- cbind(
        m1=rep(rep(c(0, 1), each=2^0), times=2^2),
        m2=rep(rep(c(0, 1), each=2^1), times=2^1),
        m3=rep(rep(c(0, 1), each=2^2), times=2^0))
    model.mat
}

format_input <- function(d, num, anno) {

    n.sample <- length(anno$sample)
    n.sub <- length(unique(anno$subject))

    feat.id <- paste0("f", seq_len(sum(num)))
    snp.id <- paste0("g", seq_len(sum(num)))
    candidate <- data.frame(feat_id=feat.id, snp_id=snp.id)

    pheno <- d %>%
        map(pluck, "y") %>%
        map(`names<-`, anno$sample) %>%
        bind_rows %>%
        as.data.frame %>%
        mutate(feat_id=feat.id) %>%
        `[`(, c(n.sample + 1, seq_len(n.sample)))

    geno <- d %>%
        map(pluck, "g") %>%
        map(`[`, !duplicated(anno$subject)) %>%
        map(`names<-`, anno$subject[!duplicated(anno$subject)]) %>%
        bind_rows %>%
        as.data.frame %>%
        mutate(snp_id=snp.id) %>%
        `[`(, c(n.sub + 1, seq_len(n.sub)))

    model <- d %>%
        map_int(pluck, "index")

    maf <- d %>%
        map_dbl(pluck, "maf")

    beta <- d %>%
        map(pluck, "beta") %>%
        bind_rows

    candidate <- candidate %>%
        mutate(model=model) %>%
        mutate(maf=maf) %>%
        bind_cols(beta)

    output <- list(
        candidate=candidate, geno=geno, pheno=pheno)

    output
}
