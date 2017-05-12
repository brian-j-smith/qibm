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
function(`_data`, f, select = NULL, ref = 1, seed = 123, ...) {
  set.seed(seed)
  
  chain <- as.matrix(chains)
  
  parmnames <- c("mu", "gamma.img", "Sigma.img", "sigma.opr", "sigma.imgopr",
                 "sigma.err", "gof.obs", "gof.rep")
  
  M <- NULL
  parms <- list()
  for (parmname in parmnames) {
    idx <- grep(paste0("^", parmname), colnames(chain))
    if (is.null(M)) M <- length(idx)
    if (is.null(select)) select <- 1:M
    parms[[parmname]] <-
      if (parmname == "gamma.img") {
        select2 <- matrix(FALSE, length(idx) / M, M)
        select2[, select] <- TRUE
        chain[, idx[select2], drop = FALSE]
      } else if (parmname == "Sigma.img") {
        select2 <- matrix(FALSE, M, M)
        select2[select, select] <- TRUE
        chain[, idx[select2], drop = FALSE]
      } else if (parmname %in% c("sigma.opr", "sigma.imgopr")) {
        parm <- matrix(0, nrow(chain), M)
        parm[, 1:length(idx)] <- chain[, idx]
        parm[, select, drop = FALSE]
      } else {
        chain[, idx[select], drop = FALSE]
      }
  }
  
  getSample <- function(i) {
    theta <- new("qibmSample", ref = ref)
    for(parmname in names(parms)) {
      slot(theta, parmname, check = FALSE) <-
        if (parmname %in% c("gamma.img", "Sigma.img")) {
          matrix(parms[[parmname]][i,], ncol = length(select))
        } else {
          parms[[parmname]][i,]
        }
    }
    theta
  }

  n <- nrow(chain)
  x <- do.call(f, c(getSample(1), list(...)))
  newchain <- matrix(NA, n, length(x))
  newchain[1,] <- x
  for(i in 2:n) {
    newchain[i,] <- do.call(f, c(getSample(i), list(...)))
  }
  
  newchains <- as.data.frame(newchain) %>%
    split(rep(1:nchain(chains), each=niter(chains))) %>%
    lapply(as.mcmc) %>%
    as.mcmc.list
  new("qibmTransform", newchains,
      select = select, ref = ref)
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
  .transform_qibm(x, GOF, ...)
})

setMethod("GOF", signature(x = "qibmSample"),
function(x) {
  as.numeric(x@gof.rep >= x@gof.obs)
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
  varnames(chains) <- paste0(c("Coef", "Var"),
                             "[", rep(chains@select, each = 2), "]")
  new("qibmLRM", chains, coef = coef, N = N)
})

.stats_LRM <- function(x, y) {
  fit <- glm.fit(cbind(1, x), y, family = binomial())
  c(fit$coef[2], chol2inv(fit$qr$qr)[2, 2])
}

setMethod("LRM", signature(x = "qibmSample"),
function(x, coef, N) {
  gamma.img <- rmvnorm(N, x@mu, x@Sigma.img, method = "chol")
  y <- rbinom(N, 1, invlogit(coef[1] + coef[2] * gamma.img[, x@ref]))
  sigma.tot <- sqrt(x@sigma.opr^2 + x@sigma.imgopr^2 + x@sigma.err^2)
  gamma.tot <- gamma.img +
    matrix(rnorm(length(gamma.img), 0, sigma.tot), N, byrow = TRUE)
  apply(gamma.tot, 2, .stats_LRM, y = y)
})
