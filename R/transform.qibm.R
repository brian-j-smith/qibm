transform.qibm <- function(chains, f, select = NULL, ref = 1, seed = 123, ...) {
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
  
  getparms <- function(i) {
    retval <- list(ref = ref)
    for(parmname in names(parms)) {
      retval[[parmname]] <- 
        if (parmname %in% c("gamma.img", "Sigma.img")) {
          matrix(parms[[parmname]][i,], ncol = length(select))
        } else {
          parms[[parmname]][i,]
        }
    }
    retval
  }

  n <- nrow(chain)
  x <- do.call(f, c(getparms(1), list(...)))
  newchain <- matrix(NA, n, length(x))
  newchain[1,] <- x
  for(i in 2:n) {
    newchain[i,] <- do.call(f, c(getparms(i), list(...)))
  }
  
  as.data.frame(newchain) %>%
    split(rep(1:nchain(chains), each=niter(chains))) %>%
    lapply(as.mcmc) %>%
    as.mcmc.list
}


bias.qibm <- function(mu, ref, ...) {
  mu - mu[ref]
}


C.qibm <- function(gamma.img, sigma.opr, sigma.imgopr, sigma.err, ref, ...) {
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
}


cor.qibm <- function(Sigma.img, diag = FALSE, ...) {
  sigma2.img <- diag(Sigma.img)
  r <- Sigma.img / sqrt(sigma2.img %o% sigma2.img)
  r[lower.tri(r, diag)]
}


gof.qibm <- function(gof.obs, gof.rep, ...) {
  as.numeric(gof.rep >= gof.obs)
}


ICC.qibm <- function(Sigma.img, sigma.opr, sigma.imgopr, sigma.err, ...) {
  sigma2.img <- diag(Sigma.img) 
  sigma2.tot <- sigma2.img + sigma.opr^2 + sigma.imgopr^2 + sigma.err^2
  sigma2.img / sigma2.tot
}


RC.qibm <- function(sigma.err, ...) {
  2.77 * sqrt(sigma.err^2)
}


RDC.qibm <- function(sigma.opr, sigma.imgopr, sigma.err, ...) {
  2.77 * sqrt(sigma.opr^2 + sigma.imgopr^2 + sigma.err^2)
}


wCV.qibm <- function(sigma.opr, sigma.imgopr, sigma.err, ...) {
  sigma2.within <- sigma.opr^2 + sigma.imgopr^2 + sigma.err^2
  sqrt(exp(sigma2.within) - 1)
}


lrm.stats <- function(x, y) {
  fit <- glm(y ~ x, family=binomial)
  c(coef(fit)[2], sqrt(vcov(fit)[2, 2]))
}

lrm.qibm <- function(mu, Sigma.img, sigma.opr, sigma.imgopr, sigma.err,
                     beta, N, ref, ...) {
  X <- rmvnorm(N, mu, Sigma.img, method = "chol")
  y <- rbinom(N, 1, invlogit(beta[1] + beta[2] * X[, ref]))
  sigma.tot <- sqrt(sigma.opr^2 + sigma.imgopr^2 + sigma.err^2)
  W <- X + matrix(rnorm(length(X), 0, sigma.tot), N, byrow = TRUE)
  apply(W, 2, lrm.stats, y = y)
}
