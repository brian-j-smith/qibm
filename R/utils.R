hpd <- function(x, alpha = 0.05, alternative = c("two.sided", "less", "greater")) {
  alternative <- match.arg(alternative)
  
  n <- length(x)
  m <- max(1, ceiling(alpha * n))
  y <- sort(x)
  
  switch(alternative,
    "two.sided" = {
      a <- y[1:m]
      b <- y[(n - m + 1):n]
      i <- which.min(b - a)
      c(a[i], b[i])
    },
    "less" = c(-Inf, y[n - m + 1]),
    "greater" = c(y[m], Inf)
  ) %>% structure(names = c("Lower", "Upper"))
}


logit <- function(x) log(x / (1 - x))

invlogit <- function(x) 1 / (1 + exp(-x))


#' @rdname describe
#' 
#' @param alpha tail probabilities for credible intervals.
#' 
setMethod("describe", signature(x = "mcmc"),
function(x, alpha = 0.05) {
  describe(as.mcmc.list(x), alpha = alpha)
})

#' @rdname describe
#' 
setMethod("describe", signature(x = "mcmc.list"),
function(x, alpha = 0.05) {
  mcvar <- sapply(x, coda:::safespec0) %>% apply(1, mean)
  f <- function(x) c(mean(x), sd(x), hpd(x, alpha = alpha))
  stats <- as.matrix(x) %>%
    apply(2, f) %>%
    t %>%
    cbind(sqrt(mcvar / (niter(x) * nchain(x)))) %>%
    as.data.frame %>%
    structure(names = c("Mean", "SD", "Lower", "Upper", "MCSE"))
  new("MCMCDescribe", stats, alpha = alpha)
})


#' Pint Values
#' 
#' @param x an object returned by \code{\link{describe}}.
#' @param digits integer indicated the number of decimal places to display.
#' @param scientific whether to encode values in scientific format.
#' 
#' @aliases print
#' 
setMethod("print", signature(x = "MCMCDescribe"),
function(x, digits = options()$digits, scientific = FALSE) {
  stats <- round(x, digits = digits) %>%
    format(trim = TRUE, nsmall = digits, scientific = scientific)
  stats$Lower <- paste0("(", stats$Lower, ", ", stats$Upper, ")")
  stats$Upper <- NULL
  names(stats)[3] <- paste0(100 * (1 - x@alpha), "% HPD")
  print(stats)
})
