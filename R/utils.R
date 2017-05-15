logit <- function(x) log(x / (1 - x))

invlogit <- function(x) 1 / (1 + exp(-x))


parminds <- function(chains, suffix) {
  grep(paste0("^", suffix), varnames(chains))
}


.describe_mcmc <- function(x, alpha, digits, scientific) {
  f <- function(x) c(mean(x), sd(x), HPDinterval(as.mcmc(x)))
  stats <- apply(as.matrix(x), 2, f) %>%
    round(digits = digits) %>%
    format(nsmall = digits, scientific = scientific) %>%
    trimws
  data.frame(stats[1,], stats[2,], paste0("(", stats[3,], ", ", stats[4,], ")"),
             stringsAsFactors = FALSE) %>%
    structure(names = c("Mean", "SD", paste0(100 * (1 - alpha), "% HPD")))
}
 
setMethod("describe", signature(x = "mcmc"),
function(x, alpha = 0.05, digits = options()$digits, scientific = FALSE) {
  .describe_mcmc(x, alpha, digits, scientific)
})

setMethod("describe", signature(x = "mcmc.list"),
function(x, alpha = 0.05, digits = options()$digits, scientific = FALSE) {
  .describe_mcmc(x, alpha, digits, scientific)
})
