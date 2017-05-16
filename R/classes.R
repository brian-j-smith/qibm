setOldClass("mcmc")
setOldClass("mcmc.list")


setClass("MCMCDescribe",
         slots = c(alpha = "numeric"),
         contains = "data.frame"
)


setClass("qibm", contains = "mcmc.list")

setClass("qibmTransform",
         slots = c(select = "numeric",
                   ref = "numeric"),
         contains = "qibm"
)

setClass("qibmGOF",
         contains = "qibmTransform"
)

setClass("qibmLRM",
         slots = c(coef = "numeric",
                   N = "numeric"),
         contains = "qibmTransform"
)

setClass("qibmSample",
  representation(
    ref = "numeric",
    mu = "numeric",
    gamma.img = "matrix",
    Sigma.img = "matrix",
    sigma.opr = "numeric",
    sigma.imgopr = "numeric",
    sigma.err = "numeric",
    gof.obs = "numeric",
    gof.rep = "numeric"
  )
)


MCMCParmVec <- R6Class("MCMCParmVec",
  public = list(
    chains = NULL,
    n = NULL,
    src = NULL,
    dest = NULL,
    initialize = function(chains, suffix, select, n) {
      inds <- parminds(chains, suffix)
      self$chains <- chains
      self$n <- length(select)
      self$dest <- which(select <= length(inds))
      self$src <- inds[select[self$dest]]
    },
    slice = function(chain, iter) {
      x <- numeric(self$n)
      x[self$dest] <- self$chains[[chain]][iter, self$src]
      x
    }
  )
)

MCMCParmColMat <- R6Class("MCMCParmMat",
  public = list(
    chains = NULL,
    n = NULL,
    inds = NULL,
    initialize = function(chains, suffix, select, n) {
      inds <- parminds(chains, suffix)
      select2 <- matrix(FALSE, length(inds) / n, n)
      select2[, select] <- TRUE
      self$chains <- chains
      self$n <- length(select)
      self$inds <- inds[select2]
    },
    slice = function(chain, iter) {
      matrix(self$chains[[chain]][iter, self$inds], ncol = self$n)
    }
  )
)

MCMCParmVar <- R6Class("MCMCParmVar",
  public = list(
    chains = NULL,
    n = NULL,
    inds = NULL,
    initialize = function(chains, suffix, select, n) {
      inds <- parminds(chains, suffix)
      select2 <- matrix(FALSE, n, n)
      select2[select, select] <- TRUE
      self$chains <- chains
      self$n <- length(select)
      self$inds <- inds[select2]
    },
    slice = function(chain, iter) {
      matrix(self$chains[[chain]][iter, self$inds], self$n)
    }
  )
)


setGeneric("describe", function(x, ...) {
  standardGeneric("describe")
})


setGeneric("Bias", function(x, ...) {
  standardGeneric("Bias")
})

setGeneric("CIndex", function(x, ...) {
  standardGeneric("CIndex")
})

setGeneric("Cor", function(x, ...) {
  standardGeneric("Cor")
})

setGeneric("GOF", function(x, ...) {
  standardGeneric("GOF")
})

setGeneric("ICC", function(x, ...) {
  standardGeneric("ICC")
})

setGeneric("RC", function(x, ...) {
  standardGeneric("RC")
})

setGeneric("RDC", function(x, ...) {
  standardGeneric("RDC")
})

setGeneric("wCV", function(x, ...) {
  standardGeneric("wCV")
})

setGeneric("LRM", function(x, ...) {
  standardGeneric("LRM")
})
