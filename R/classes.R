setOldClass("mcmc")
setOldClass("mcmc.list")


setClass("MCMCDescribe",
         slots = c(alpha = "numeric"),
         contains = "data.frame"
)


setClass("qibm", 
         slots = c(M = "numeric"),
         contains = "mcmc.list"
)

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


MCMCParmVec <- R6Class("MCMCParmVec",
  public = list(
    chains = NULL,
    src = NULL,
    dest = NULL,
    n = NULL,
    initialize = function(chains, src, select, n) {
      self$chains <- chains
      self$dest <- which(select <= length(src))
      self$src <- src[select[self$dest]]
      self$n <- length(select)
    },
    slice = function(chain, iter) {
      x <- numeric(self$n)
      x[self$dest] <- self$chains[[chain]][iter, self$src]
      x
    }
  )
)

MCMCParmColMat <- R6Class("MCMCParmColMat",
  public = list(
    chains = NULL,
    src = NULL,
    dim = NULL,
    initialize = function(chains, src, select, n) {
      m <- length(src) / n
      select2 <- matrix(FALSE, m, n)
      select2[, select] <- TRUE
      self$chains <- chains
      self$src <- src[select2]
      self$dim <- c(m, length(select))
    },
    slice = function(chain, iter) {
      matrix(self$chains[[chain]][iter, self$src], self$dim[1], self$dim[2])
    }
  )
)

MCMCParmVarMat <- R6Class("MCMCParmVarMat",
  public = list(
    chains = NULL,
    src = NULL,
    n = NULL,
    initialize = function(chains, src, select, n) {
      select2 <- matrix(FALSE, n, n)
      select2[select, select] <- TRUE
      self$chains <- chains
      self$src <- src[select2]
      self$n <- length(select)
    },
    slice = function(chain, iter) {
      matrix(self$chains[[chain]][iter, self$src], self$n, self$n)
    }
  )
)

MCMCParmList <- R6Class("MCMCParmList",
  public = list(
    parms = NULL,
    initialize = function(chains, select, n) {
      tokens <- strsplit(varnames(chains), "[][,]")
      prefix <- sapply(tokens, function(x) x[1])
      parms <- list()
      for (parmname in unique(prefix)) {
        src <- which(parmname == prefix)
        parms[[parmname]] <- switch(parmname,
          "mu" = MCMCParmVec$new(chains, src, select, n),
          "gamma.img" = MCMCParmColMat$new(chains, src, select, n),
          "Sigma.img" = MCMCParmVarMat$new(chains, src, select, n),
          "sigma.opr" = MCMCParmVec$new(chains, src, select, n),
          "sigma.imgopr" = MCMCParmVec$new(chains, src, select, n),
          "sigma.err" = MCMCParmVec$new(chains, src, select, n),
          "gof.obs" = MCMCParmVec$new(chains, src, select, n),
          "gof.rep" = MCMCParmVec$new(chains, src, select, n),
          NULL
        )
      }
      self$parms <- parms
    },
    slice = function(chain, iter) {
      lapply(self$parms, function(x) x$slice(chain, iter))
    }
  )
)


#' Descriptive Summaries
#' 
#' @param x an \R object.
#' @param ... further arguments passed to or from other methods.
#' 
setGeneric("describe", function(x, ...) {
  standardGeneric("describe")
})


#' Relative Mean Bias
#' 
#' @param x an \R object.
#' @param ... further arguments passed to \code{\link{with}}.
#' 
setGeneric("Bias", function(x, ...) {
  standardGeneric("Bias")
})

#' Concordance Index
#' 
#' @param x an \R object.
#' @param ... further arguments passed to \code{\link{with}}.
#' 
setGeneric("CIndex", function(x, ...) {
  standardGeneric("CIndex")
})

#' Between-Method Correlation
#' 
#' @param x an \R object.
#' @param ... further arguments passed to \code{\link{with}}.
#' 
setGeneric("Cor", function(x, ...) {
  standardGeneric("Cor")
})

#' Model Goodness of Fit
#' 
#' @param x an \R object.
#' @param ... further arguments passed to \code{\link{with}}.
#' 
setGeneric("GOF", function(x, ...) {
  standardGeneric("GOF")
})

#' Intra-Class Correlation Coefficient
#' 
#' @param x an \R object.
#' @param ... further arguments passed to \code{\link{with}}.
#' 
setGeneric("ICC", function(x, ...) {
  standardGeneric("ICC")
})

#' Repeatability Coefficient
#' 
#' @param x an \R object.
#' @param ... further arguments passed to \code{\link{with}}.
#' 
setGeneric("RC", function(x, ...) {
  standardGeneric("RC")
})

#' Reproducibility Coefficient
#' 
#' @param x an \R object.
#' @param ... further arguments passed to \code{\link{with}}.
#' 
setGeneric("RDC", function(x, ...) {
  standardGeneric("RDC")
})

#' Within-Subject Coefficient of Variation
#' 
#' @param x an \R object.
#' @param ... further arguments passed to \code{\link{with}}.
#' 
setGeneric("wCV", function(x, ...) {
  standardGeneric("wCV")
})

#' Logistic Regression Model Estimates
#' 
#' @param x an \R object.
#' @param ... further arguments passed to \code{\link{with}}.
#' 
setGeneric("LRM", function(x, ...) {
  standardGeneric("LRM")
})
