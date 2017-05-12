#' Quantitative Imaging Biomarker Model
#' 
#' \code{qibm} is used to fit a Bayesian measurement error model for the 
#' comparison of biomarker measurements derived from different quantitative
#' imaging methods.
#'
#' @param fixed fomula of the form \code{y ~ method} where \code{y} is a
#'   vector of biomarker measurements, or a transformation thereof, and
#'   \code{method} is a vector of identifiers for the imaging methods.
#' @param image quoted or unquoted name of a vector of image identifiers.
#' @param operator quoted or unquoted name of a vector of operator identifiers.
#' @param data data.frame containing the data vectors.
#' @param priors list of prior parameter values.
#' @param parameters vector naming the model parameters to be returned.
#' @param n.burnin number of MCMC samples to discard as a burn-in sequence.
#' @param n.iter total number of samples to generate.
#' @param n.thin period at which to save samples.
#' @param n.chains number of parallel MCMC chains to generate.
#' @param seed numeric random number generator seed.
#' 
#' @return A \code{qibm} object that inherits from \code{\link[coda]{mcmc.list}} and
#'   contains the MCMC sampled model parameter values.

qibm <- function(fixed, image, operator, data, priors = list(),
                 parameters = c("mu", "sigma.opr", "sigma.imgopr", "sigma.err"),
                 n.burnin = 2500, n.iter = 10000, n.thin = 5, n.chains = 3,
                 seed = 123) {
  ### Set random number seed
  set.seed(seed)
  
  ### Data structures
  mf <- model.frame(fixed, data=data, na.action=NULL)
  
  methodvar <- attr(terms(fixed), "term.labels")
  imgvar <- evalq(as.character(substitute(image)))
  oprvar <- evalq(as.character(substitute(operator)))
  
  data <- data.frame(
    y = model.response(mf),
    method = as.factor(data[[methodvar]]),
    img = as.factor(data[[imgvar]]),
    opr = as.factor(data[[oprvar]])
  ) %>% na.omit

  idx <- tapply(data$opr, data$method, function(x) length(unique(x)) > 1)
  model2method <- levels(data$method)[idx]
  inmodel2 <- data$method %in% model2method
  
  opr <- factor(ifelse(inmodel2, data$opr, NA))
  imgopr <- factor(ifelse(inmodel2, paste0(data$img, ":", data$opr), NA))
  
  ### Function to extract prior parameter values
  getprior <- function(param, default) {
    user <- priors[[param]]
    if(is.null(user)) default else user
  }

  ### Model data
  model.data <- within(list(), {
    y <- data$y
    N <- nrow(data)
    
    model1IDX <- which(!inmodel2)
    model1N <- length(model1IDX)
    model2IDX <- which(inmodel2)
    model2N <- length(model2IDX)
    
    methodIDX <- as.numeric(data$method)
    methodN <- nlevels(data$method)
    imgIDX <- as.numeric(data$img)
    imgN <- nlevels(data$img)
    
    oprIDX <- as.numeric(opr)
    oprN <- nlevels(opr)
    opr2method <- methodIDX[match(levels(opr), opr)]
    
    imgoprIDX <- as.numeric(imgopr)
    imgoprN <- nlevels(imgopr)
    imgopr2method <- methodIDX[match(levels(imgopr), imgopr)]
  
    methodvarIDX <- which(levels(data$method) %in% model2method)
    methodvarN <- length(methodvarIDX)
    
    mu.mean <- getprior("mu.mean", 0.0)
    mu.var <- getprior("mu.var", 1.0e6)
    Sigma.img.scale <- getprior("Sigma.img.scale", 1.0) *
      (diag(methodN) + getprior("Sigma.img.cor", 0.0) * (1 - diag(methodN)))
    Sigma.img.df <- getprior("Sigma.img.df", methodN)
    sigma2.opr.shape <- getprior("sigma2.opr.shape", 1.0e-3)
    sigma2.opr.rate <- getprior("sigma2.opr.rate", 1.0e-3)
    sigma2.imgopr.shape <- getprior("sigma2.imgopr.shape", 1.0e-3)
    sigma2.imgopr.rate <- getprior("sigma2.imgopr.rate", 1.0e-3)
    sigma2.err.shape <- getprior("sigma2.err.shape", 1.0e-3)
    sigma2.err.rate <- getprior("sigma2.err.rate", 1.0e-3)
    
    gof <- array(0, dim=c(N, methodN, 2))
    Map(function(i) gof[i, methodIDX[i], ] <<- NA, 1:N)
  })
  
  ### Initial values generator
  inits <- function(chain) {
    stats <- data.frame(y = model.data$y,
                        method = model.data$methodIDX,
                        img = model.data$imgIDX) %>%
             ddply(.(method, img), summarize,
                   y = y,
                   z = y - mean(y)) %>%
             ddply(.(method), summarize,
                   mean = mean(y),
                   se = sd(y) / sqrt(length(y)),
                   prec = pmin(1 / var(z), 1e10),
                   hw = sqrt(sqrt(prec * 12) / 2))
    stats$prec1 <- pmax(stats$prec - stats$hw, 0)
    stats$prec2 <- stats$prec1 + 2 * stats$hw
    within(list(), {
      mu <- rnorm(model.data$methodN, stats$mean, stats$se)
      tau.imgopr <- rep(NA, max(model.data$methodvarIDX))
      tau.imgopr[model.data$methodvarIDX] <-
        runif(model.data$methodvarN, stats$prec1, stats$prec2)
      tau.err <- runif(model.data$methodN, stats$prec1, stats$prec2)
      .RNG.name <- "base::Wichmann-Hill"
      .RNG.seed <- chain
    })
  }
  
  ### JAGS model specification
  qibm.jags <- "
    model{
      for (i in 1:model1N) {
        y[model1IDX[i]] ~ dnorm(gamma[model1IDX[i]], tau.err[methodIDX[model1IDX[i]]])
        gamma[model1IDX[i]] <- gamma.img[imgIDX[model1IDX[i]], methodIDX[model1IDX[i]]]
      }
      
      for (i in 1:model2N) {
        y[model2IDX[i]] ~ dnorm(gamma[model2IDX[i]], tau.err[methodIDX[model2IDX[i]]])
        gamma[model2IDX[i]] <- gamma.img[imgIDX[model2IDX[i]], methodIDX[model2IDX[i]]] +
                               gamma.opr[oprIDX[model2IDX[i]]] +
                               gamma.imgopr[imgoprIDX[model2IDX[i]]]
      }
      
      for (i in 1:N) {
        y.rep[i] ~ dnorm(gamma[i], tau.err[methodIDX[i]])
        gof[i, methodIDX[i], 1] <- pow(y.rep[i] - gamma[i], 2) * tau.err[methodIDX[i]]
        gof[i, methodIDX[i], 2] <- pow(y[i] - gamma[i], 2) * tau.err[methodIDX[i]]
      }  

      for (i in 1:imgN) {
        gamma.img[i, 1:methodN] ~ dmnorm(mu[], Omega.img[,])
      }
      
      for (i in 1:oprN) {
        gamma.opr[i] ~ dnorm(0.0, tau.opr[opr2method[i]])
      }
      
      for (i in 1:imgoprN) {
        gamma.imgopr[i] ~ dnorm(0.0, tau.imgopr[imgopr2method[i]])
      }
      
      for (i in 1:methodvarN) {
        tau.opr[methodvarIDX[i]] ~ dgamma(sigma2.opr.shape, sigma2.opr.rate)
        sigma.opr[methodvarIDX[i]] <- 1 / sqrt(tau.opr[methodvarIDX[i]])
        tau.imgopr[methodvarIDX[i]] ~ dgamma(sigma2.imgopr.shape, sigma2.imgopr.rate)
        sigma.imgopr[methodvarIDX[i]] <- 1 / sqrt(tau.imgopr[methodvarIDX[i]])
      }
      
      for (j in 1:methodN) {
        mu[j] ~ dnorm(mu.mean, 1 / mu.var)
        tau.err[j] ~ dgamma(sigma2.err.shape, sigma2.err.rate)
        sigma.err[j] <- 1 / sqrt(tau.err[j])
        
        gof.rep[j] <- sum(gof[, j, 1])
        gof.obs[j] <- sum(gof[, j, 2])
      }
      
      Omega.img[1:methodN, 1:methodN] ~ dwish(Sigma.img.scale[,], Sigma.img.df)
      Sigma.img[1:methodN, 1:methodN] <- inverse(Omega.img[,])
    }
  "
  
  ### Model fit
  chains <- jags.model(textConnection(qibm.jags), model.data, inits, n.chains) %>%
    coda.samples(parameters, n.iter, n.thin) %>%
    window(n.burnin + 1)
  new("qibm", chains)
}
