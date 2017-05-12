setOldClass("mcmc")
setOldClass("mcmc.list")

setClass("qibm", contains="mcmc.list")

setClass("qibmTransform",
         slots = c(select = "numeric",
                   ref = "numeric"),
         contains="qibm"
)

setClass("qibmLRM",
         slots = c(coef = "numeric",
                   N = "numeric"),
         contains="qibmTransform"
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


setGeneric("summarize")


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
