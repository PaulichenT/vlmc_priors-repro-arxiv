# ------------------------------------------------------------
# Bayes-factor based methods
# ------------------------------------------------------------
# This script implements the Bayes factor for context trees and
# presents both methods for selecting the maximal depth and
# performing model selection.
# ------------------------------------------------------------


#' @description
#' Computes the base-10 logarithm of the Bayes factor comparing two models
#' specified by context-tree functions
#'
#' @details
#' Given two context-tree functions `f` and `g`, the Bayes factor is
#' defined as the ratio of their marginal likelihoods under the same data and
#' hyperparameters. This function returns the logarithm (base 10) of this ratio,
#' computed as the difference between the log10 marginal likelihoods.
#' The marginal likelihoods are obtained via the method `getMarginalLikelihood`
#' of the BaContrees package.
#'
#' @param data A character vector representing the observed sequence.
#' @param maximalDepth Integer. Maximum depth allowed in the tree space.
#' @param alpha Numeric. Hyperparameter of the Bayesian context tree model.
#' @param f A context-tree function.
#' @param g A context-tree function.
#'
#' @returns
#' A numeric value corresponding to the base-10 logarithm of the Bayes factor
#' in favor of model induced by `f` over model induced by `g`.
#'
log_bayes_factor <- function(data, maximalDepth, alpha, f, g) {

  btf <- baConTree$new(data, maximalDepth, alpha, f)
  btg <- baConTree$new(data, maximalDepth, alpha, g)

  log_evidence_f <- btf$getMarginalLikelihood(TRUE) / log(10)
  log_evidence_g <- btg$getMarginalLikelihood(TRUE) / log(10)

  log_evidence_f - log_evidence_g
}


#' @description
#' Selects the maximal depth of a Bayesian context tree model using a sequential
#' Bayes factor-based procedure.
#'
#' @details
#' Starting from the maximal depth `maximalDepth`, the function compares
#' nested models of decreasing depth using the log10 Bayes factor. At each step,
#' the model with current depth `current_L` is compared against a model
#' with smaller depth `k`. If the log Bayes factor is below the specified
#' `threshold`, the smaller model is preferred and `current_L` is updated.
#' The procedure continues until all depths are evaluated, returning the
#' selected depth.
#'
#' @param data A character vector representing the observed sequence.
#' @param maximalDepth Integer. Maximum depth considered in the model.
#' @param threshold Numeric. Threshold for the log10 Bayes factor (default is 0).
#'
#' @returns
#' A list with two elements:
#' - `sample_size` The length of the input data.
#' - `selected_depth` The selected maximal depth.
#'
select_depth_bf <- function(data, maximalDepth = 10, threshold = 0) {

  current_L <- maximalDepth

  for (k in seq(maximalDepth - 1, 0, by = -1)) {
    bf <- log_bayes_factor(data, maximalDepth, alpha = 0.5, d_l(current_L), d_l(k))
    if (bf < threshold) {
      current_L <- k
    }
  }

  list(
    sample_size = length(data),
    selected_depth = current_L
  )
}


#' @description
#' Constructs the product of two context-tree functions.
#'
#' @details
#' Given two context-tree functions `f` and `g`, this function
#' returns a new context-tree function whose value at each node is the
#' product of the values returned by `f` and `g` at that node.
#'
#' @param `f` A context-tree function.
#' @param `g` A context-tree function.
#'
#' @returns
#' A context-tree function representing the product `fg`
#'
multiply_ctfunctions <- function(f, g) {
  function(node) {
    f(node) * g(node)
  }
}

#' @description
#' Selects the maximal depth of a Bayesian context tree model associated
#' with a given context-tree function `f`, using a sequential Bayes factor-based
#' procedure.
#'
#' @details
#' Starting from the maximal depth `maximalDepth`, the function compares
#' nested models of decreasing depth using the log10 Bayes factor. At each step,
#' the model with current depth `current_L` is compared against a model
#' with smaller depth `k`. If the log Bayes factor is below the specified
#' `threshold`, the smaller model is preferred and `current_L` is updated.
#'
#' The comparison is performed by combining the l-depth context-tree function
#' `d_l` with the context-tree function `f`.
#'
#' @param data A character vector representing the observed sequence.
#' @param maximalDepth Integer. Maximum depth considered in the model.
#' @param alpha Numeric. Hyperparameter of the Bayesian context tree model.
#' @param f A context-tree function.
#' @param threshold Numeric. Threshold for the log10 Bayes factor (default is 0).
#'
#' @returns
#' An integer corresponding to the selected maximal depth.
#'
select_depth_f_bf <- function(data, maximalDepth = 10, alpha = 0.5, f) {

  current_L <- maximalDepth

  for (k in seq(maximalDepth - 1, 0, by = -1)) {

    f_current   <- multiply_ctfunctions(d_l(current_L), f)
    f_candidate <- multiply_ctfunctions(d_l(k), f)

    bf <- log_bayes_factor(data, maximalDepth, alpha, f_current, f_candidate)
    if (bf < 0) {current_L <- k}
  }

  current_L
}


#' @description
#' Performs model selection among a set of context-tree functions using
#' a two-stage Bayes factor procedure.
#'
#' @details
#' The procedure consists of two steps:
#'
#' - Step 1 (depth selection):
#' For each prior (context-tree function) in `priors`, the optimal maximal
#' depth is selected using `select_depth_f_bf`, resulting in a pair
#' `(f, l*)` for each context-tree function.
#'
#' - Step 2 (model comparison):
#' The resulting models are compared via pairwise log10 Bayes factors. Each model
#' is defined by combining the l-depth function at the selected depth `d_l*`
#' with the corresponding context-tree function `f`. The procedure sequentially
#' updates the current best model whenever another model yields a larger marginal
#' likelihood.
#'
#' @param data A character vector representing the observed sequence.
#' @param maximalDepth Integer. Maximum depth considered in the model.
#' @param alpha Numeric. Hyperparameter of the Bayesian context tree model.
#' @param priors A named list of context-tree functions.
#'
#' @returns
#' A list with two elements:
#' - `best_f` The selected context-tree function (prior).
#' - `best_depth` The selected maximal depth associated with the prior.
#'
model_selection <- function(data, maximalDepth = 10, alpha = 0.5, priors) {

  m <- length(priors)

  # Part 1:
  depth_results <- vector("list", m)

  for (i in seq_len(m)) {
    res <- select_depth_f_bf(data, maximalDepth, alpha, f = priors[[i]])

    depth_results[[i]] <- list(
      prior  = priors[[i]],
      l_star = res
    )
  }

  # Part 2:
  best_idx <- 1

  for (i in 2:m) {

    f_best <- multiply_ctfunctions(
      d_l(depth_results[[best_idx]]$l_star),
      depth_results[[best_idx]]$prior
    )

    f_i <- multiply_ctfunctions(d_l(depth_results[[i]]$l_star),
                           depth_results[[i]]$prior)

    bf <- log_bayes_factor(data, maximalDepth, alpha, f_best, f_i)
    if (bf < 0) {best_idx <- i}
  }

  list(
    best_f = priors[best_idx],
    best_depth = depth_results[[best_idx]]$l_star
  )
}







