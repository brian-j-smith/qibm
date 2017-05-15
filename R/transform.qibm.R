#' Transform Parameters from a Quantitative Imaging Biomarker Model Fit
#' 
#' \code{transform} is used to transform model parameters from a \link{qibm} 
#' fit.
#'
#' @param _data \code{\link{qibm}} object of MCMC sample model parameter
#'   values.
#' @param f transform function to be applied to the model parameters.
#' @param select vector of quantification methods to select for the transform
#'   (default: all).
#' @param ref identifies the reference quantification method as \code{select[ref]}.
#' @param seed numeric random number generator seed.
#' @param ... further arguments to be passed to \code{f}.
#' 
#' @return An \code{\link[coda]{mcmc.list}} object containing transformed MCMC
#'   sampled model parameter values.

setMethod("transform", signature(`_data` = "qibm"),
function(`_data`, f, select, ref = 1, seed = 123, ...) {
  set.seed(seed)
  
  M <- length(parminds(chains, "mu"))
  if(missing(select)) select <- 1:M
  
  parms <- list()
  parms$mu <- MCMCParmVec$new(chains, "mu", select, M)
  parms$gamma.img <- MCMCParmColMat$new(chains, "gamma.img", select, M)
  parms$Sigma.img <- MCMCParmVar$new(chains, "Sigma.img", select, M)
  parms$sigma.opr <- MCMCParmVec$new(chains, "sigma.opr", select, M)
  parms$sigma.imgopr <- MCMCParmVec$new(chains, "sigma.imgopr", select, M)
  parms$sigma.err <- MCMCParmVec$new(chains, "sigma.err", select, M)
  parms$gof.obs <- MCMCParmVec$new(chains, "gof.obs", select, M)
  parms$gof.rep <- MCMCParmVec$new(chains, "gof.rep", select, M)
  
  getSample <- function(chain, iter) {
    theta <- new("qibmSample", ref = ref)
    for (parmname in names(parms)) {
      slot(theta, parmname, check = FALSE) <-
        parms[[parmname]]$slice(chain, iter)
    }
    theta
  }

  x <- do.call(f, c(getSample(1, 1), list(...)))
  newchains <- as.mcmc.list(
    replicate(nchain(chains),
              mcmc(matrix(NA, niter(chains), length(x))),
              simplify = FALSE)
  )
  
  for (k in 1:nchain(chains)) {
    for (i in 1:niter(chains)) {
      newchains[[k]][i, ] <- do.call(f, c(getSample(k, i), list(...)))
    }
  }

  new("qibmTransform", newchains, select = select, ref = ref)
})


.transform_qibm <- function(x, f, ...) {
  chains <- transform(x, f, ...)
  prefix <- eval(as.character(substitute(f)))
  varnames(chains) <- paste0(prefix, "[", chains@select, "]")
  chains
}  


setMethod("Bias", signature(x = "qibm"),
function(x, ...) {
  .transform_qibm(x, Bias, ...)
})

setMethod("Bias", signature(x = "qibmSample"),
function(x) {
  x@mu - x@mu[x@ref]
})


setMethod("CIndex", signature(x = "qibm"),
function(x, ...) {
  .transform_qibm(x, CIndex, ...)
})

setMethod("CIndex", signature(x = "qibmSample"),
function(x) {
  M <- ncol(x@gamma.img)
  N <- nrow(x@gamma.img)
  retval <- rep(NA, M)
  sigma.within <- sqrt(x@sigma.opr^2 + x@sigma.imgopr^2 + x@sigma.err^2)
  b.ref <- rnorm(N, x@gamma.img[, x@ref], sigma.within[x@ref])
  retval[x@ref] <- 1
  for(i in setdiff(1:M, x@ref)) {
    b <- rnorm(N, x@gamma.img[, i], sigma.within[i])
    retval[i] <- rcorr.cens(b, b.ref, outx = TRUE)["C Index"]
  }
  retval
})


setMethod("Cor", signature(x = "qibm"),
function(x, diag = FALSE, ...) {
  chains <- transform(x, Cor, diag = diag, ...)
  N <- length(chains@select)
  ind <- matrix(chains@select, N, N)
  varnames(chains) <- paste0("Cor[", t(ind)[lower.tri(ind, diag = diag)], ",",
                             ind[lower.tri(ind, diag = diag)], "]")
  chains
})

setMethod("Cor", signature(x = "qibmSample"),
function(x, diag = FALSE, ...) {
  sigma2.img <- diag(x@Sigma.img)
  r <- x@Sigma.img / sqrt(sigma2.img %o% sigma2.img)
  r[lower.tri(r, diag)]
})


setMethod("GOF", signature(x = "qibm"),
function(x, ...) {
  chains <- transform(x, GOF, ...)
  tags <- expand.grid(chains@select, c("gof.rep", "gof.obs"))
  varnames(chains) <- paste0(tags[, 2], "[", tags[, 1], "]")
  new("qibmGOF", chains)
})

setMethod("GOF", signature(x = "qibmSample"),
function(x) {
  c(x@gof.rep,  x@gof.obs)
})

setMethod("describe", signature(x = "qibmGOF"),
function(x, alpha = 0.05, digits = options()$digits, scientific = FALSE) {
  chains <- transform(x, function(x) as.numeric(x@gof.rep >= x@gof.obs))
  describe(chains, alpha = alpha, digits = digits, scientific = scientific)
})

setMethod("plot", signature(x = "qibmGOF"),
function(x, y, ...) {
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


setMethod("ICC", signature(x = "qibm"),
function(x, ...) {
  .transform_qibm(x, ICC, ...)
})

setMethod("ICC", signature(x = "qibmSample"),
function(x) {
  sigma2.img <- diag(x@Sigma.img) 
  sigma2.tot <- sigma2.img + x@sigma.opr^2 + x@sigma.imgopr^2 + x@sigma.err^2
  sigma2.img / sigma2.tot
})


setMethod("RC", signature(x = "qibm"),
function(x, ...) {
  .transform_qibm(x, RC, ...)
})

setMethod("RC", signature(x = "qibmSample"),
function(x) {
  2.77 * sqrt(x@sigma.err^2)
})


setMethod("RDC", signature(x = "qibm"),
function(x, ...) {
  .transform_qibm(x, RDC, ...)
})

setMethod("RDC", signature(x = "qibmSample"),
function(x) {
  2.77 * sqrt(x@sigma.opr^2 + x@sigma.imgopr^2 + x@sigma.err^2)
})


setMethod("wCV", signature(x = "qibm"),
function(x, ...) {
  .transform_qibm(x, wCV, ...)
})

setMethod("wCV", signature(x = "qibmSample"),
function(x) {
  sigma2.within <- x@sigma.opr^2 + x@sigma.imgopr^2 + x@sigma.err^2
  sqrt(exp(sigma2.within) - 1)
})


setMethod("LRM", signature(x = "qibm"),
function(x, coef, N, ...) {
  chains <- transform(x, LRM, coef = coef, N = N, ...)
  tags <- expand.grid(c("Coef", "Var"), chains@select, N)
  varnames(chains) <- paste0(tags[, 1], "[", tags[, 2], ",", tags[, 3], "]")
  new("qibmLRM", chains, coef = coef, N = N)
})

.stats_LRM <- function(x, y) {
  fit <- glm.fit(cbind(1, x), y, family = binomial())
  c(fit$coef[2], chol2inv(fit$qr$qr)[2, 2])
}

setMethod("LRM", signature(x = "qibmSample"),
function(x, coef, N) {
  stats <- matrix(NA, 2 * length(x@mu), length(N))
  for (j in seq_along(N)) {
    n <- N[j]
    gamma.img <- rmvnorm(n, x@mu, x@Sigma.img, method = "chol")
    y <- rbinom(n, 1, invlogit(coef[1] + coef[2] * gamma.img[, x@ref]))
    sigma.tot <- sqrt(x@sigma.opr^2 + x@sigma.imgopr^2 + x@sigma.err^2)
    gamma.tot <- gamma.img +
      matrix(rnorm(length(gamma.img), 0, sigma.tot), n, byrow = TRUE)
    stats[, j] <- apply(gamma.tot, 2, .stats_LRM, y = y)
  }
  stats
})
