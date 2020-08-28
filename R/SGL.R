#' Function to compute the SGL solution
#'
#' Originally written by Simon for R package SGL.

SGLinner <- function(data, group.length, type = "linear", maxit = 1000, thresh = 0.001, min.frac = 0.1,
                gamma = 0.8, verbose = FALSE, step = 1, reset = 10,
                lambda1 = 0.2, lambda2 = 0.01, groupW = NULL, betas.old = rep(0, ncol(data$x)),
                weights = rep(1, NROW(data$x)))
{
  intercept = 0
  if (type == "linear")
  {
    Sol <- run_sgl_linear(data, group.length, thresh = thresh, inner.iter = maxit, outer.iter = maxit, outer.thresh = thresh,
                          min.frac = min.frac, lambda1 = lambda1, lambda2=lambda2,
                          gamma = gamma, verbose = verbose, step = step, reset = reset,
                          groupW = groupW, beta.naught = betas.old, weights = weights)
    print(Sol)
    Sol = list(beta = Sol$beta, intercept = intercept)
    print(Sol)
  } else if (type == "logit")
  {
    Sol <- run_sgl_logit(data, group.length, thresh = thresh, inner.iter = maxit, outer.iter = maxit,
                         outer.thresh = thresh, min.frac = min.frac, lambda1 = lambda1, lambda2= lambda2,
                         gamma = gamma, verbose = verbose, step = step, reset = reset,
                         groupW = groupW, beta.naught = betas.old, weights = weights)
    Sol = list(beta = Sol$beta, intercept = Sol$intercept)
  }
  return(Sol)
}




SGL2 <- function (data, index, type = "linear", maxit = 1000, thresh = 0.001, 
                  min.frac = 0.1, nlam = 20, gamma = 0.8, standardize = TRUE, 
                  verbose = FALSE, step = 1, reset = 10, alpha = 0.95, lambdas = NULL,
                  weights = rep(1, NROW(data$x))) 
{
  
  temp = index2group.length(index)
  group.length = temp$group.length
  
  num.groups = length(group.length)
  range.group.ind = c(0, cumsum(group.length))
  index = group.length2index(group.length)
  
  groupW <- NULL
  
  weights <- weights / mean(weights)
  
  if( is.null(groupW) )
  {
    groupW = sqrt(group.length)
  }
  
  X.transform <- NULL
  if (standardize == TRUE) {
    X <- data$x
    means <- apply(X, 2, mean)
    X <- t(t(X) - means)
    var <- apply(X, 2, function(x) (sqrt(sum(x^2))))
    X <- t(t(X)/var)
    data$x <- X
    X.transform <- list(X.scale = var, X.means = means)
  }
  if (type == "linear") {
    if (standardize == TRUE) {
      intercept <- mean(data$y)
      data$y <- data$y - intercept
    }
    Sol <- oneDim2(data, index, thresh, inner.iter = maxit, 
                  outer.iter = maxit, outer.thresh = thresh, min.frac = min.frac, 
                  nlam = nlam, lambdas = lambdas, gamma = gamma, verbose = verbose, 
                  step = step, reset = reset, alpha = alpha, groupW = groupW, weights = weights)
    if (standardize == TRUE) {
      Sol <- list(beta = Sol$beta, lambdas = Sol$lambdas, 
                  type = type, intercept = intercept, X.transform = X.transform)
    }
    if (standardize == FALSE) {
      Sol <- list(beta = Sol$beta, lambdas = Sol$lambdas, 
                  type = type, X.transform = X.transform)
    }
  }
  if (type == "logit") {
    Sol <- oneDimLogit2(data, index, thresh = thresh, inner.iter = maxit, 
                       outer.iter = maxit, outer.thresh = thresh, min.frac = min.frac, 
                       nlam = nlam, lambdas = lambdas, gamma = gamma, verbose = verbose, 
                       step = step, alpha = alpha, reset = reset, groupW = groupW, weights = weights)
    Sol <- list(beta = Sol$beta, lambdas = Sol$lambdas, type = type, 
                intercept = Sol$intercept, X.transform = X.transform)
  }
  if (type == "cox") {
    Sol <- oneDimCox2(data, index, thresh = thresh, inner.iter = maxit, 
                     outer.iter = maxit, outer.thresh = thresh, min.frac = min.frac, 
                     nlam = nlam, lambdas = lambdas, gamma = gamma, verbose = verbose, 
                     step = step, alpha = alpha, reset = reset)
    Sol = list(beta = Sol$beta, lambdas = Sol$lambdas, type = type, 
               X.transform = X.transform)
  }
  class(Sol) = "SGL"
  return(Sol)
}



oneDimLogit2 <- function (data, index, thresh = 1e-04, lambdas = NULL, beta.naught = rep(0, ncol(data$x)), inner.iter = 100, outer.iter = 100, outer.thresh = 1e-04, 
                          gamma = 0.8, step = 1, reset = 10, alpha = 0.95, min.frac = 0.05, 
                          nlam = 20, verbose = FALSE, groupW, weights = rep(1, NROW(data$x)) ) 
{
  if (is.null(lambdas)) {
    lambdas <- SGL:::betterPathCalc(data = data, index = index, 
                              alpha = alpha, min.frac = min.frac, nlam = nlam, 
                              type = "logit")
  }
  X <- data$x
  y <- data$y
  n <- nrow(X)
  p <- ncol(X)
  ord <- order(index)
  index <- index[ord]
  X <- X[, ord]
  unOrd <- match(1:length(ord), ord)
  groups <- unique(index)
  num.groups <- length(groups)
  range.group.ind <- rep(0, (num.groups + 1))
  for (i in 1:num.groups) {
    range.group.ind[i] <- min(which(index == groups[i])) - 
      1
  }
  range.group.ind[num.groups + 1] <- ncol(X)
  group.length <- diff(range.group.ind)
  beta.naught <- rep(0, ncol(X))
  beta <- beta.naught
  beta.is.zero <- rep(1, num.groups)
  beta.old <- rep(0, ncol(X))
  beta <- matrix(0, nrow = ncol(X), ncol = nlam)
  eta <- rep(0, n)
  intercepts <- rep(log(sum(y)) - log(n - sum(y)), nlam)
  eta = eta + intercepts[1]
  beta.is.zero <- rep(1, num.groups)
  beta.old <- rep(0, ncol(X))
  for (i in 1:nlam) {
    junk <- .C("logitNest", X = as.double(as.vector(X)), 
               y = as.integer(y), index = as.integer(index), nrow = as.integer(nrow(X)), 
               ncol = as.integer(ncol(X)), numGroup = as.integer(num.groups), 
               rangeGroupInd = as.integer(range.group.ind), groupLen = as.integer(group.length), 
               lambda1 = as.double(alpha * lambdas[i]), lambda2 = as.double((1 - alpha) * lambdas[i]), beta = as.double(beta.old), 
               innerIter = as.integer(inner.iter), outerIter = as.integer(outer.iter), 
               thresh = as.double(thresh), outerThresh = as.double(outer.thresh), 
               eta = as.double(eta), gamma = as.double(gamma), betaIsZero = as.integer(beta.is.zero), 
               betaZero = as.double(intercepts[i]), step = as.double(step),
               weights = as.double(weights), PACKAGE = "sglfast")
    intercepts[i] = junk$betaZero
    if (i < nlam) {
      intercepts[i + 1] = intercepts[i]
    }
    beta.new <- junk$beta
    beta[, i] <- beta.new
    beta.is.zero <- junk$betaIsZero
    eta <- junk$eta
    beta.old <- beta.new
    if (verbose == TRUE) {
      write(paste("***Lambda", i, "***"), "")
    }
  }
  return(list(beta = beta[unOrd, ], lambdas = lambdas, intercepts = intercepts))
}


oneDim2 <- function (data, index, thresh = 1e-04, nlam = 20, lambdas = NULL, 
                     beta.naught = rep(0, ncol(data$x)), inner.iter = 100, outer.iter = 100, 
                     outer.thresh = 1e-04, gamma = 0.8, step = 1, reset = 10, 
                     alpha = 0.95, min.frac = 0.05, verbose = FALSE,
                     groupW,
                     weights = rep(1, NROW(data$x))) 
{
  if (is.null(lambdas)) 
  {
    lambdas <- betterPathCalc2(data = data, index = index, 
                               alpha = alpha, min.frac = min.frac, nlam = nlam, 
                               type = "linear")
  }
  X <- data$x
  y <- data$y
  n <- nrow(X)
  p <- ncol(X)
  ord <- order(index)
  index <- index[ord]
  X <- X[, ord]
  unOrd <- match(1:length(ord), ord)
  groups <- unique(index)
  num.groups <- length(groups)
  range.group.ind <- rep(0, (num.groups + 1))
  for (i in 1:num.groups) {
    range.group.ind[i] <- min(which(index == groups[i])) - 
      1
  }
  range.group.ind[num.groups + 1] <- ncol(X)
  group.length <- diff(range.group.ind)
  nlam = length(lambdas)
  beta.old <- rep(0, ncol(X))
  beta.is.zero <- rep(1, num.groups)
  beta <- array(0, c(ncol(X), nlam))
  eta <- rep(0, n)
  for (k in 1:nlam) {
    beta.is.zero <- rep(1, num.groups)
    beta.old <- rep(0, ncol(X))
    eta <- rep(0, n)
    junk <- .C("linNest", X = as.double(as.vector(X)), y = as.double(y), 
               index = as.integer(index), nrow = as.integer(nrow(X)), 
               ncol = as.integer(ncol(X)), numGroup = as.integer(num.groups), 
               rangeGroupInd = as.integer(range.group.ind), groupLen = as.integer(group.length), 
               lambda1 = as.double(lambdas[k] * alpha), lambda2 = as.double(lambdas[k] * (1 - alpha)), 
               beta = as.double(beta.old), innerIter = as.integer(inner.iter), 
               outerIter = as.integer(outer.iter), thresh = as.double(thresh), 
               outerThresh = as.double(outer.thresh), eta = as.double(eta), 
               gamma = as.double(gamma), betaIsZero = as.integer(beta.is.zero), 
               step = as.double(step), reset = as.integer(reset), 
               weights = as.double(weights), PACKAGE = "sglfast")
    beta.new <- junk$beta
    beta[, k] <- beta.new
    beta.is.zero <- junk$betaIsZero
    eta <- junk$eta
    beta.old <- beta.new
    if (verbose == TRUE) {
      write(paste("***Lambda", k, "***"), "")
    }
  }
  return(list(beta = beta[unOrd, ], lambdas = lambdas))
}



betterPathCalc2 <- function (data, index, alpha = 0.95, min.frac = 0.05, nlam = 20, 
                             type = "linear") 
{
  reset <- 10
  step <- 1
  gamma <- 0.8
  inner.iter <- 1000
  outer.iter <- 1000
  thresh = 10^(-3)
  outer.thresh = thresh
  n <- nrow(data$x)
  if (type == "linear") {
    X <- data$x
    resp <- data$y
    n <- nrow(X)
    p <- ncol(X)
    ord <- order(index)
    index <- index[ord]
    X <- X[, ord]
    unOrd <- match(1:length(ord), ord)
    groups <- unique(index)
    num.groups <- length(groups)
    range.group.ind <- rep(0, (num.groups + 1))
    for (i in 1:num.groups) {
      range.group.ind[i] <- min(which(index == groups[i])) - 
        1
    }
    range.group.ind[num.groups + 1] <- ncol(X)
    group.length <- diff(range.group.ind)
  }
  if (type == "logit") {
    X <- data$x
    y <- data$y
    n <- nrow(X)
    p <- ncol(X)
    ord <- order(index)
    index <- index[ord]
    X <- X[, ord]
    unOrd <- match(1:length(ord), ord)
    groups <- unique(index)
    num.groups <- length(groups)
    range.group.ind <- rep(0, (num.groups + 1))
    for (i in 1:num.groups) {
      range.group.ind[i] <- min(which(index == groups[i])) - 
        1
    }
    range.group.ind[num.groups + 1] <- ncol(X)
    group.length <- diff(range.group.ind)
    beta.naught <- rep(0, ncol(X))
    beta <- beta.naught
    beta.is.zero <- rep(1, num.groups)
    beta.old <- rep(0, ncol(X))
    betas <- matrix(0, nrow = ncol(X), ncol = nlam)
    eta <- rep(0, n)
    intercepts <- mean(y)
    eta = eta + intercepts
    m.y <- mean(y)
    resp <- m.y * m.y * (1 - m.y) - (y - m.y)
  }
  if (type == "cox") {
    covariates <- data$x
    n <- nrow(covariates)
    p <- ncol(covariates)
    time <- data$time
    status <- data$status
    death.order <- order(time)
    ordered.time <- sort(time)
    X <- covariates[death.order, ]
    ordered.status <- status[death.order]
    first.blood <- min(which(ordered.status == 1))
    X <- X[first.blood:n, ]
    ordered.status <- ordered.status[first.blood:n]
    ordered.time <- ordered.time[first.blood:n]
    death.order <- death.order[first.blood:n]
    n <- n - first.blood + 1
    death.times <- unique(ordered.time[which(ordered.status == 
                                               1)])
    risk.set <- rep(0, n)
    for (i in 1:n) {
      risk.set[i] <- max(which(death.times <= ordered.time[i]))
    }
    risk.set.ind <- rep(0, (length(death.times) + 1))
    for (i in 1:length(death.times)) {
      risk.set.ind[i] <- min(which(ordered.time >= death.times[i]))
    }
    risk.set.ind[length(risk.set.ind)] <- length(ordered.time) + 
      1
    num.deaths <- rep(0, length(death.times))
    for (i in 1:length(ordered.time)) {
      if (ordered.status[i] == 1) {
        num.deaths[which(death.times == ordered.time[i])] <- num.deaths[which(death.times == 
                                                                                ordered.time[i])] + 1
      }
    }
    death.index <- which(ordered.status == 1)
    total.deaths <- length(death.index)
    ord <- order(index)
    index <- index[ord]
    X <- X[, ord]
    unOrd <- match(1:length(ord), ord)
    groups <- unique(index)
    num.groups <- length(groups)
    range.group.ind <- rep(0, (num.groups + 1))
    for (i in 1:num.groups) {
      range.group.ind[i] <- min(which(index == groups[i])) - 
        1
    }
    range.group.ind[num.groups + 1] <- ncol(X)
    group.length <- diff(range.group.ind)
    beta.naught <- rep(0, ncol(X))
    beta <- beta.naught
    beta.is.zero <- rep(1, num.groups)
    beta.old <- rep(0, ncol(X))
    beta <- array(0, c(ncol(X), nlam, nlam))
    beta.is.zero <- rep(1, num.groups)
    eta <- rep(0, n)
    junk1 <- .C("Cox", riskSetInd = as.integer(risk.set.ind), 
                riskSet = as.integer(risk.set), numDeath = as.integer(num.deaths), 
                status = as.integer(ordered.status), ndeath = as.integer(length(death.times)), 
                nrow = as.integer(n), ncol = as.integer(p), beta = as.double(rep(0, 
                                                                                 p)), eta = as.double(rep(0, n)), y = as.double(rep(0, 
                                                                                                                                    n)), weights = as.double(rep(0, n)))
    resp <- junk1$y * junk1$weights
  }
  lambda.max <- rep(0, num.groups)
  if ((alpha != 0) * (alpha != 1)) {
    for (i in 1:num.groups) {
      ind <- groups[i]
      X.fit <- X[, which(index == ind)]
      cors <- t(X.fit) %*% resp
      ord.cors <- sort(abs(cors), decreasing = TRUE)
      if (length(ord.cors) > 1) {
        norms <- rep(0, length(cors) - 1)
        lam <- ord.cors/alpha
        for (j in 1:(length(ord.cors) - 1)) {
          norms[j] <- sqrt(sum((ord.cors[1:j] - ord.cors[j + 
                                                           1])^2))
        }
        if (norms[1] >= lam[2] * (1 - alpha) * sqrt(group.length[i])) {
          our.cors <- ord.cors[1]
          our.range <- c(ord.cors[2], ord.cors[1])/alpha
        }
        else {
          if (norms[length(ord.cors) - 1] <= lam[length(ord.cors)] * 
              (1 - alpha) * sqrt(group.length[i])) {
            our.cors <- ord.cors
            our.range <- c(0, ord.cors[length(ord.cors)])/alpha
          }
          else {
            my.ind <- max(which(norms[-length(norms)] <= 
                                  lam[2:(length(norms))] * (1 - alpha) * 
                                  sqrt(group.length[i]))) + 1
            our.cors <- ord.cors[1:my.ind]
            our.range <- c(ord.cors[my.ind + 1], ord.cors[my.ind])/alpha
          }
        }
        nn <- length(our.cors)
        if (alpha == 0.5) {
          alpha = 0.500001
        }
        A.term <- nn * alpha^2 - (1 - alpha)^2 * group.length[i]
        B.term <- -2 * alpha * sum(our.cors)
        C.term <- sum(our.cors^2)
        lams <- c((-B.term + sqrt(B.term^2 - 4 * A.term * 
                                    C.term))/(2 * A.term), (-B.term - sqrt(B.term^2 - 
                                                                             4 * A.term * C.term))/(2 * A.term))
        lambda.max[i] <- min(subset(lams, lams >= our.range[1] & 
                                      lams <= our.range[2]))
      }
      if (length(ord.cors) == 1) {
        lambda.max[i] <- ord.cors
      }
    }
  }
  if (alpha == 1) {
    lambda.max <- abs(t(X) %*% resp)
  }
  if (alpha == 0) {
    for (i in 1:num.groups) {
      ind <- groups[i]
      X.fit <- X[, which(index == ind)]
      cors <- t(X.fit) %*% resp
      lambda.max[i] <- sqrt(sum(cors^2))/sqrt(group.length[i])
    }
  }
  max.lam <- max(lambda.max)
  min.lam <- min.frac * max.lam
  lambdas <- exp(seq(log(max.lam), log(min.lam), (log(min.lam) - 
                                                    log(max.lam))/(nlam - 1)))
  return(lambdas/nrow(X))
}
