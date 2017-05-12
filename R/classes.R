setOldClass("mcmc")
setOldClass("mcmc.list")

setClass("qibm", contains="mcmc.list")

setClass("qibmTransform",
         slots = c(select = "numeric",
                   ref = "numeric"),
         contains="qibm"
)

setClass("qibmLRM",
         slots = c(beta = "numeric",
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


setGeneric("describe", function(object, ...) {
  standardGeneric("describe")
})


setGeneric("Bias", function(object, ...) {
  standardGeneric("Bias")
})

setGeneric("CIndex", function(object, ...) {
  standardGeneric("CIndex")
})

setGeneric("Cor", function(object, ...) {
  standardGeneric("Cor")
})

setGeneric("GOF", function(object, ...) {
  standardGeneric("GOF")
})

setGeneric("ICC", function(object, ...) {
  standardGeneric("ICC")
})

setGeneric("RC", function(object, ...) {
  standardGeneric("RC")
})

setGeneric("RDC", function(object, ...) {
  standardGeneric("RDC")
})

setGeneric("wCV", function(object, ...) {
  standardGeneric("wCV")
})

setGeneric("LRM", function(object, ...) {
  standardGeneric("LRM")
})
