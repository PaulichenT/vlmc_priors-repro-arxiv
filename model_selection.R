# ------------------------------------------------------------
# Bayesian model comparison
# ------------------------------------------------------------
# This script implements the model selection methods
# ------------------------------------------------------------

#' @description
#' Evaluates a list of priors in the Bayesian model and returns a vector of
#' log10 marginal likelihood values.
#'
#' @details
#' For each weight function provided in `priors`, the log10 marginal likelihood
#' is computed and stored with the corresponding prior. 
#'
#' @param data A character vector representing the observed sequence.
#' @param maximalDepth Maximum depth (integer). 
#' @param alpha Dirichlet hyperparameter.
#' @param priors A list of weight functions.
#'
#' @returns
#' A numeric vector containing the base-10 log marginal likelihood values corresponding 
#' to each weight function in `priors`.
log_evidence <- function(data, maximalDepth, alpha, priors) {
  sapply(priors, function(w) {
    bt <- baConTree$new(data, maximalDepth = maximalDepth, alpha = alpha, priorWeights = w)
    bt$getMarginalLikelihood() / log(10)
  })
}

#' @description
#' Evaluates a range of candidate maximal depths for the Bayesian context tree
#' model and selects the depth that maximizes the log posterior score.
#'
#' @details
#' The function iterates through candidate maximal tree depths from 1 to `L_max`.
#' For each depth L, it computes the log marginal likelihood and adds the log 
#' prior probability derived from `prior_L(L)`, obtaining the log posterior score. 
#' The optimal depth `best_L` corresponds to the candidate L achieving the maximum 
#' total log score.
#'
#' @param data A vector representing the observed sequence data.
#' @param L_max The upper limit of candidate maximal depths (integer). Default is 10.
#' @param priorWeights Weight function `w`. Default is the unit weight `u`.
#' @param prior_L Prior probability function over {1, 2, ..., L_max}. 
#' Defaults to the uniform prior `function(L) 1 / L_max`.
#'
#' @returns
#' A named list containing:
#' - `best_L`: The depth L that maximizes the score. 
#' - `log_values`: A named numeric vector containing the computed log posterior
#'  scores for each candidate depth.
#'  
select_depth <- function(data, L_max = 10, priorWeights = u, prior_L = function(L) 1 / (L_max)) {
  log_values <- numeric(L_max)
  
  for (L in 1:L_max) {
    # Compute the evidence for a fixed L
    bt <- baConTree$new(data, maximalDepth = L, alpha = 0.5, priorWeights = priorWeights)
    
    # Compute the log posterior score
    log_values[L] <- bt$getMarginalLikelihood() + log(prior_L(L))
  }
  
  # Find the L that maximizes the log posterior score
  best_L <- (1:L_max)[which.max(log_values)]
  
  return(list(
    best_L = best_L,
    log_values = setNames(log_values, 1:L_max)
  ))
}

