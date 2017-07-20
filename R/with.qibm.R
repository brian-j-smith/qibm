#' Evaluate an Expression in Quantitative Imaging Biomarker Model MCMC Samples
#' 
#' Evaluate an \R expression in an environment constructed from MCMC
#' samples, possibly modifying (a copy of) the original samples.
#'
#' @param data \code{\link{qibm}} object of MCMC sample model parameter values.
#' @param expr expression to evaluate.
#' @param select vector of quantification methods to select for the evaluation
#'   (default: all).
#' @param ref identifies the reference quantification method as \code{select[ref]}.
#' @param ... further arguments to include in the environment.
#' 
#' @return An \code{\link[coda]{mcmc.list}} object of MCMC sampled values of the
#' evaluated \code{expr}.
#' 
#' @seealso \code{\link{qibm}}.
#' @aliases with
#' 
with.qibm <- function(data, expr, select, ref = 1, ...) {
  if(missing(select)) select <- 1:data@M

  parms <- MCMCParmList$new(data, select, data@M)
  
  x <- eval(substitute(expr), c(parms$slice(1, 1), list(ref = ref, ...)))
  chains <- as.mcmc.list(
    replicate(nchain(data),
              mcmc(matrix(NA, niter(data), length(x)),
                   start = start(data), thin = thin(data)),
              simplify = FALSE)
  )
  
  for (i in 1:nchain(chains)) {
    for (j in 1:niter(chains)) {
      chains[[i]][j, ] <-
        eval(substitute(expr), c(parms$slice(i, j), list(ref = ref, ...)))
    }
  }
  
  new("qibmTransform", chains, M = data@M, select = select, ref = ref)
}


setnames <- function(chains, name) {
  varnames(chains) <- paste0(name, "[", chains@select, "]")
  chains
}


#' @rdname Bias
#' 
setMethod("Bias", signature(x = "qibm"),
function(x, ...) {
  with(x, mu - mu[ref], ...) %>% setnames("Bias")
})


#' @rdname CIndex
#' 
setMethod("CIndex", signature(x = "qibm"),
function(x, ...) {
  with(x, {
    M <- ncol(gamma.img)
    N <- nrow(gamma.img)
    retval <- rep(NA, M)
    sigma.within <- sqrt(sigma.opr^2 + sigma.imgopr^2 + sigma.err^2)
    b.ref <- rnorm(N, gamma.img[, ref], sigma.within[ref])
    retval[ref] <- 1
    for(i in setdiff(1:M, ref)) {
      b <- rnorm(N, gamma.img[, i], sigma.within[i])
      retval[i] <- rcorr.cens(b, b.ref, outx = TRUE)["C Index"]
    }
    retval
  }, ...) %>% setnames("CIndex")
})


#' @rdname Cor
#' 
#' @param diag whether to include correlation matrix diagonals.
#' 
setMethod("Cor", signature(x = "qibm"),
function(x, diag = FALSE, ...) {
  chains <- with(x, {
    sigma2.img <- diag(Sigma.img)
    r <- Sigma.img / sqrt(sigma2.img %o% sigma2.img)
    r[lower.tri(r, diag)]
  }, diag = diag, ...)
  N <- length(chains@select)
  ind <- matrix(chains@select, N, N)
  varnames(chains) <- paste0("Cor[", t(ind)[lower.tri(ind, diag = diag)], ",",
                             ind[lower.tri(ind, diag = diag)], "]")
  chains
})


#' @rdname GOF
#' 
setMethod("GOF", signature(x = "qibm"),
function(x, ...) {
  chains <- with(x , c(gof.rep, gof.obs), ...)
  tags <- expand.grid(chains@select, c("gof.rep", "gof.obs"))
  varnames(chains) <- paste0(tags[, 2], "[", tags[, 1], "]")
  new("qibmGOF", chains)
})

#' @rdname GOF
#' 
#' @param alpha tail probabilities for credible intervals.
#' 
setMethod("describe", signature(x = "qibmGOF"),
function(x, alpha = 0.05) {
  chains <- with(x, as.numeric(gof.rep >= gof.obs), select = seq_along(x@select))
  varnames(chains) <- paste0("pval[", x@select, "]")
  describe(chains, alpha = alpha)
})

#' @describeIn GOF Plot observed versus replicated goodness of fit quantities.
#' 
setMethod("plot", signature(x = "qibmGOF"),
function(x) {
  data <- as.data.frame(as.matrix(x))
  
  vars <- names(data)
  n <- length(vars) / 2
  datalong <- reshape(data, varying = list(head(vars, n), tail(vars, n)),
                      v.names = c("gof.rep", "gof.obs"), times = x@select,
                      direction = "long")
  
  pval <- tapply(datalong$gof.rep >= datalong$gof.obs, datalong$time, mean)
  ann_text <- data.frame(
    gof.obs = -Inf,
    gof.rep = Inf,
    time = x@select,
    lab = paste0("p = ", round(pval, 3))
  )
  
  ggplot(datalong, aes(gof.obs, gof.rep)) +
    geom_abline(slope = 1, intercept = 0) +
    geom_point() +
    labs(x = expression(T(y*"|"*theta)),
         y = expression(T(y^{rep}*"|"*theta))) +
    geom_text(data = ann_text, aes(label = lab, hjust = 0, vjust = 1)) +
    facet_wrap(~ time, ncol = ceiling(sqrt(n)))
})


#' @rdname ICC
#' 
setMethod("ICC", signature(x = "qibm"),
function(x, ...) {
  with(x, {
    sigma2.img <- diag(Sigma.img) 
    sigma2.tot <- sigma2.img + sigma.opr^2 + sigma.imgopr^2 + sigma.err^2
    sigma2.img / sigma2.tot
  }, ...) %>% setnames("ICC")
})


#' @rdname RC
#' 
setMethod("RC", signature(x = "qibm"),
function(x, ...) {
  with(x, 2.77 * sigma.err, ...) %>% setnames("RC")
})


#' @rdname RDC
#' 
setMethod("RDC", signature(x = "qibm"),
function(x, ...) {
  with(x, 2.77 * sqrt(sigma.opr^2 + sigma.imgopr^2 + sigma.err^2), ...) %>%
    setnames("RDC")
})


#' @rdname wCV
#' 
setMethod("wCV", signature(x = "qibm"),
function(x, ...) {
  with(x, {
    sigma2.within <- sigma.opr^2 + sigma.imgopr^2 + sigma.err^2
    sqrt(exp(sigma2.within) - 1)
  }, ...) %>% setnames("wCV")
})


#' @rdname LRM
#' 
#' @param coef two-element vector specifying the logistic regression
#'   intercept and slope with which to simulate outcomes.
#' @param N scalar or vector of simulated sample sizes.
#' @param seed random number seed for the simulations.
#' 
#' @seealso \code{\link{LRMCoef}}.
#' 
setMethod("LRM", signature(x = "qibm"),
function(x, coef, N, seed = 123, ...) {
  set.seed(seed)
  chains <- with(x, {
    stats <- matrix(NA, 2 * length(mu), length(N))
    for (j in seq_along(N)) {
      n <- N[j]
      gamma.img <- rmvnorm(n, mu, Sigma.img, method = "chol")
      y <- rbinom(n, 1, invlogit(coef[1] + coef[2] * gamma.img[, ref]))
      sigma.tot <- sqrt(sigma.opr^2 + sigma.imgopr^2 + sigma.err^2)
      gamma.tot <- gamma.img +
        matrix(rnorm(length(gamma.img), 0.0, sigma.tot), n, byrow = TRUE)
      stats[, j] <- apply(gamma.tot, 2, .stats_LRM, y = y)
    }
    stats
  }, coef = coef, N = N, ...)
  tags <- expand.grid(c("Coef", "Var"), chains@select, N)
  varnames(chains) <- paste0(tags[, 1], "[", tags[, 2], ",", tags[, 3], "]")
  new("qibmLRM", chains, coef = coef, N = N)
})

.stats_LRM <- function(x, y) {
  fit <- glm.fit(cbind(1, x), y, family = binomial())
  c(fit$coef[2], chol2inv(fit$qr$qr)[2, 2])
}

#' @rdname LRM
#' 
#' @param alpha significance for statistical testing of the odds ratio.
#' @param alternative direction of the alternative hypothesis.
#' @param scale value at which to compute the odds ratio
#' 
setMethod("describe", signature(x = "qibmLRM"),
function(x, alpha = 0.05, alternative = c("two.sided", "less", "greater"),
         scale = 1) {
  alternative <- match.arg(alternative)
  alpha2 <- if (alternative == "two.sided") alpha / 2 else alpha
  
  coef.ref <- x@coef[2]
  
  f <- function(x) {
    inds <- seq(1, ncol(x), by = 2)
    coef <- x[, inds, drop = FALSE]
    se <- sqrt(x[, -inds, drop = FALSE])
    err <- qnorm(1 - alpha2) * se
    innull <- 1
    incoverage <- 1
    if (alternative != "less") {
      lcl <- coef - err
      innull <- innull * (lcl < 0.0)
      incoverage <- incoverage * (lcl < coef.ref)
    }
    if (alternative != "greater") {
      ucl <- coef + err
      innull <- innull * (0.0 < ucl)
      incoverage <- incoverage * (coef.ref < ucl)
    }
    cbind(1.0 - apply(innull, 2, mean),
          coda:::safespec0(innull),
          apply(incoverage, 2, mean),
          coda:::safespec0(incoverage))
  }
  
  id <- expand.grid(Method = x@select, N = x@N)
  
  coef <- as.matrix(x)[, seq(1, nvar(x), by = 2), drop = FALSE]
  or <- exp(scale * coef)
  or.mean <- apply(or, 2, mean)
  or.hpd <- apply(or, 2, hpd, alpha = alpha, alternative = alternative)
  rmse <- sqrt((or.mean - exp(scale * coef.ref))^2 + apply(or, 2, var))
  
  stats <- lapply(x, f) %>% simplify2array %>% apply(c(1, 2), mean)
  
  data.frame(
    id,
    OR = or.mean,
    Lower = or.hpd[1, ],
    Upper = or.hpd[2, ],
    RMSE = rmse,
    Power = stats[, 1],
    Power.MCSE = sqrt(stats[, 2] / (niter(x) * nchain(x))),
    Coverage = stats[, 3],
    Coverage.MCSE = sqrt(stats[, 4] / (niter(x) * nchain(x))),
    row.names = 1:nrow(id)
  )
})
