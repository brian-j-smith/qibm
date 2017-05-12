logit <- function(x) log(x / (1 - x))

invlogit <- function(x) 1 / (1 + exp(-x))


.describe_mcmc <- function(object, alpha, digits, scientific) {
  f <- function(x) c(mean(x), sd(x), HPDinterval(as.mcmc(x)))
  stats <- apply(as.matrix(object), 2, f) %>%
    round(digits = digits) %>%
    format(nsmall = digits, scientific = scientific) %>%
    trimws
  data.frame(stats[1,], stats[2,], paste0("(", stats[3,], ", ", stats[4,], ")"),
             stringsAsFactors = FALSE) %>%
    structure(names = c("Mean", "SD", paste0(100 * (1 - alpha), "% HPD")))
}
 
setMethod("describe", signature(object = "mcmc"),
function(object, alpha = 0.05, digits = options()$digits, scientific = FALSE) {
  .describe_mcmc(object, alpha, digits, scientific)
})

setMethod("describe", signature(object = "mcmc.list"),
function(object, alpha = 0.05, digits = options()$digits, scientific = FALSE) {
  .describe_mcmc(object, alpha, digits, scientific)
})
