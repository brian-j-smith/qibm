logit <- function(x) log(x / (1 - x))

invlogit <- function(x) 1 / (1 + exp(-x))


parminds <- function(chains, suffix) {
  grep(paste0("^", suffix), varnames(chains))
}


setMethod("describe", signature(x = "mcmc"),
function(x, alpha = 0.05) {
  describe(as.mcmc.list(x), alpha = alpha)
})

setMethod("describe", signature(x = "mcmc.list"),
function(x, alpha = 0.05) {
  mcse <- sapply(x, coda:::safespec0) %>% apply(1, mean)
  f <- function(x) c(mean(x), sd(x), HPDinterval(as.mcmc(x), prob = 1 - alpha))
  stats <- as.matrix(x) %>%
    apply(2, f) %>%
    t %>%
    cbind(mcse) %>%
    as.data.frame %>%
    structure(names = c("Mean", "SD", "Lower", "Upper", "MCSE"))
  new("MCMCDescribe", stats, alpha = alpha)
})


setMethod("print", signature(x = "MCMCDescribe"),
function(x, digits = options()$digits, scientific = FALSE) {
  stats <- round(x, digits = digits) %>%
    format(trim = TRUE, nsmall = digits, scientific = scientific)
  stats$Lower <- paste0("(", stats$Lower, ", ", stats$Upper, ")")
  stats$Upper <- NULL
  names(stats)[3] <- paste0(100 * (1 - x@alpha), "% HPD")
  print(stats)
})
