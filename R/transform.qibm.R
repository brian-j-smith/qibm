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
  
  as.data.frame(newchain) %>%
    split(rep(1:nchain(chains), each=niter(chains))) %>%
    lapply(as.mcmc) %>%
    as.mcmc.list
})


setMethod("Bias", signature(object = "qibmSample"),
function(object) {
  object@mu - object@mu[object@ref]
})


setMethod("CIndex", signature(object = "qibmSample"),
function(object) {
  M <- ncol(object@gamma.img)
  N <- nrow(object@gamma.img)
  retval <- rep(NA, M)
  sigma.within <- sqrt(object@sigma.opr^2 + object@sigma.imgopr^2 +
                         object@sigma.err^2)
  b.ref <- rnorm(N, object@gamma.img[, object@ref],
                 sigma.within[object@ref])
  retval[object@ref] <- 1
  for(i in setdiff(1:M, object@ref)) {
    b <- rnorm(N, object@gamma.img[, i], sigma.within[i])
    retval[i] <- rcorr.cens(b, b.ref, outx = TRUE)["C Index"]
  }
  retval
})


setMethod("Cor", signature(object = "qibmSample"),
function(object, diag = FALSE, ...) {
  sigma2.img <- diag(object@Sigma.img)
  r <- object@Sigma.img / sqrt(sigma2.img %o% sigma2.img)
  r[lower.tri(r, diag)]
})


setMethod("GOF", signature(object = "qibmSample"),
function(object) {
  as.numeric(object@gof.rep >= object@gof.obs)
})


setMethod("ICC", signature(object = "qibmSample"),
function(object) {
  sigma2.img <- diag(object@Sigma.img) 
  sigma2.tot <- sigma2.img + object@sigma.opr^2 + object@sigma.imgopr^2 +
    object@sigma.err^2
  sigma2.img / sigma2.tot
})


setMethod("RC", signature(object = "qibmSample"),
function(object) {
  2.77 * sqrt(object@sigma.err^2)
})


setMethod("RDC", signature(object = "qibmSample"),
function(object) {
  2.77 * sqrt(object@sigma.opr^2 + object@sigma.imgopr^2 + object@sigma.err^2)
})


setMethod("wCV", signature(object = "qibmSample"),
function(object) {
  sigma2.within <- object@sigma.opr^2 + object@sigma.imgopr^2 +
    object@sigma.err^2
  sqrt(exp(sigma2.within) - 1)
})


setMethod("LRM", signature(object = "qibmSample"),
function(object, beta, N) {
  X <- rmvnorm(N, object@mu, object@Sigma.img, method = "chol")
  y <- rbinom(N, 1, invlogit(beta[1] + beta[2] * X[, object@ref]))
  sigma.tot <- sqrt(object@sigma.opr^2 + object@sigma.imgopr^2 +
                      object@sigma.err^2)
  W <- X + matrix(rnorm(length(X), 0, sigma.tot), N, byrow = TRUE)
  apply(W, 2, LRMstats, y = y)
})

LRMstats <- function(x, y) {
  fit <- glm(y ~ x, family = binomial)
  c(coef(fit)[2], sqrt(vcov(fit)[2, 2]))
}
