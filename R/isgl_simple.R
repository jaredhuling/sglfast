#' Two-parameter version of the iterative sparse-group lasso
#'
#' Computes the solution of the sparse-group lasso problem, defined as,
#' \deqn{\hat{\beta}(\lambda, \gamma)=\arg\min_{\beta\in B} \left\{ \hat{R}(\beta) + \lambda_2 \sum_{j=1}^{J} \gamma_j \|\beta^{(j)}\|_2 +  \lambda_1 \|\beta\|_1 \right\},}
#' where \eqn{\gamma_j} is the square root of j-th group's size.
#' Regularization parameter \eqn{\lambda} is automatically selected, so that the validation error is minimized.
#'
#' @param data.train A list with the matrix of observations \code{x}, and the response \code{y}, to train the model.
#' @param data.validate A list with the matrix of observations \code{x}, and the response \code{y}, to validate the model.
#' @param index A vector with group indices.
#' @param group.length A vector with group length, alternatively to \code{index}.
#' @param type String: \code{'linear'} or \code{'logit'}.
#' @param standardize Whether to stardardize \code{data.train} and \code{data.validate}. Recommended if \code{type} is \code{'logit'}.
#' @param momentum Acceleration rate. Should be > 1.
#' @return An object of class \code{isgl}.

isgl_simple = function( data.train, data.validate, index = NULL, group.length = NULL, type = "linear",
                        standardize = FALSE, momentum = 2, weights = rep(1, NROW(data.train))) 
{

  # We tranform the initial data
  if (standardize){
    temp = standarize_inital_data(data.train, data.validate, type)
    data.train = temp$data.train
    data.validate =  temp$data.validate
    X.transform = temp$X.transform
    intercept = temp$intercept
  }
  else{
    X.transform = list(X.means = rep(0, ncol(data.train$x)), X.scale = rep(1, ncol(data.train$x)))
    intercept = 0
  }

  # We sort and compute group lengths
  if( is.null(group.length) )
  {
    if( is.null(index) )
    {
      stop("Error 1: You must provide valid group indices or lengths")
    }
    temp = index2group.length(index)
    group.length = temp$group.length
    ord = temp$ord
    data.train$x = data.train$x[, ord]
    data.validate$x = data.validate$x[, ord]
    unord = match(1:length(ord),ord)
  } else
  {
    unord = 1:ncol(data.train$x)
  }
  

  # Compute initial lambdas
  lambda.max = c(0, 0)
  gamma = rep(0, length(group.length))
  lambda.max[1] = get_lambda1.max(data.train, group.length, type = type, weights = weights)

  lambda.init = c(lambda.max[1]*0.1, 1)
  for (i in 1:length(gamma)) {
    gamma[i] = get_gammak.max(data.train, group.length, i, type, lambda.init, weights = weights)
  }

  lambda.max[2] = max(gamma/sqrt(group.length))
  
  # Start the iterative search
  nparams = 2
  best_lambdas <- lambda.max*0.1
  num_solves <- 0
  max_solves <- 5000

  # Compute initial model
  model_params = solve_inner_problem(data.train, group.length, best_lambdas, type, simple = TRUE, weights = weights)
  best_cost = get_validation_cost(data.validate$x, data.validate$y, model_params, type)
  best_beta = model_params
  num_solves = num_solves+1
  
  # Set initial coordinate
  coord <- 1
  fixed = 0

  # Main loop
  while ((num_solves < max_solves) &&(fixed < nparams) ) {
    old_lambdas <- best_lambdas

    # No need to compute lambda1,2 max?

    curr_lambdas = best_lambdas
    t0 = best_lambdas[coord]*0.1

    # Direction: ->
    dir = 1
    t = t0
    if(t == 0){t = 0.01}
    while (dir >= -1) {
      curr_lambdas[coord] = best_lambdas[coord] + dir*runif(1, 0.1*t, t)
      model_params <- solve_inner_problem(data.train, group.length, curr_lambdas, type, simple = TRUE, weights = weights)
      num_solves <- num_solves + 1
      cost <- get_validation_cost( data.validate$x, data.validate$y, model_params, type)

      if(best_cost - cost > 0.00001 ){
        best_cost = cost
        best_lambdas <- curr_lambdas
        best_beta = model_params
        if(dir==1){t = momentum*t}else{t = min(momentum*t, curr_lambdas[coord])}
      }else{
        dir = dir - 2
        t = t0
      }
    }

    if(old_lambdas[coord] == best_lambdas[coord]){
      fixed = fixed + 1
    }else{
      fixed = 0
    }
    coord <- coord%%(nparams) + 1
  }

  solution = list(best_lambdas=best_lambdas, num_solves=num_solves, best_cost=best_cost,
                  beta = best_beta$beta[unord],
                  intercept = best_beta$intercept + intercept,
                  X.transform = X.transform,
                  type = type)
  class(solution) = "isgl"
  return( solution )
}
