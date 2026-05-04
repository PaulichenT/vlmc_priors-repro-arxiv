# ------------------------------------------------------------
# Auxiliary functions
# ------------------------------------------------------------
# This script implements auxiliary functions for the simulation
# studies of the paper.
# ------------------------------------------------------------

#' @description
#' Generates a list of samples from a variable-length Markov chain model
#' for a set of specified sample sizes.
#'
#' @details
#' For each value in `sample_sizes`, a sequence of length `n` is
#' generated using the `rvlmc()` function, according to the provided
#' context tree.
#'
#' @param sample_sizes A numeric vector of sample sizes.
#' @param true_contexts A vector specifying the contexts of the true model.
#' @param true_probs A list of transition probability vectors associated
#' with each context.
#'
#' @returns
#' A named list of samples, where each element corresponds to a generated
#' sequence of length given by `sample_sizes`. The names of the list
#' elements are of the form "n_<size>".
#'
generate_samples <- function(sample_sizes, true_contexts, true_probs) {

  samples <- lapply(sample_sizes, function(n) {

    set.seed(1)
    rvlmc(
      n = n,
      alphabet = c("0", "1"),
      context_list = true_contexts,
      context_probs = true_probs
    )
  })

  names(samples) <- paste0("n_", sample_sizes)

  return(samples)
}


#' @description
#' Fits the Bayesian context tree model to a sample and computes evaluation metrics.
#'
#' @details
#' Given a sample, a context-tree function `f`, and the true contexts, the model is
#' fitted using the bacontrees package. The function computes:
#' - the distance between the true tree and the MAP tree;
#' - the prior and posterior probabilities of the true tree; and
#' - the log10 marginal likelihood of the model.
#'
#' @param sample A character vector representing the observed sequence.
#' @param f A context-tree function (prior weights).
#' @param true_contexts A vector specifying the contexts of the true tree.
#'
#' @returns
#' A list containing:
#' - `delta_map` distance between true and MAP tree;
#' - `prior_true_tree` prior probability of the true tree;
#' - `post_true_tree` posterior probability of the true tree;
#' - `log_10_evidence` log10 marginal likelihood.
#'
fit_single_model <- function(sample, f, true_contexts) {

  # Fit Bayesian context tree model
  bt <- baConTree$new(sample, 10, 0.5, f)

  # Obtain the value of evidence
  evidence <- bt$getMarginalLikelihood(TRUE)

  # Activate the true tree and obtain its probabilities
  bt$activateFromContexts(true_contexts)
  probs <- bt$activeTreeProbabilities()
  prior_true_tree <- probs[1]
  post_true_tree  <- probs[2]

  # Obtain the contexts of the true tree
  true_tree_contexts <- bt$getActiveNodes()

  # Activate the MAP tree and obtain its contexts
  map_tree_contexts  <- bt$activateMap()$getActiveNodes()

  # Function for computing the structural distance
  tree_distance <- function(contexts1, contexts2) {
    inner_nodes1 <- bt$activateFromContexts(contexts1)$getInnerNodes()
    inner_nodes2 <- bt$activateFromContexts(contexts2)$getInnerNodes()

    length(setdiff(inner_nodes1, inner_nodes2)) +
      length(setdiff(inner_nodes2, inner_nodes1))
  }

  # Distance between the true tree and the MAP tree
  delta_map <- tree_distance(true_tree_contexts, map_tree_contexts)

  list(
    delta_map = delta_map,
    prior_true_tree = prior_true_tree,
    post_true_tree = post_true_tree,
    log_10_evidence = evidence / log(10)
  )
}


#' @description
#' Evaluates a set of prior specifications on a given sample using the
#' Bayesian context tree model.
#'
#' @details
#' For each context-tree function in `priors`, the model is fitted to the provided
#' `sample` via `fit_single_model()`. The function returns, for each prior,
#' the corresponding evaluation metrics.
#'
#' @param sample A character vector representing the observed sequence.
#' @param priors A named list of prior (context-tree) functions.
#' @param true_contexts A vector specifying the contexts of the true tree.
#'
#' @returns
#' A list of results, one for each context-tree function in `priors`.
#' Each element is itself a list containing the output of `fit_single_model()`.
#'
fit_priors_on_sample <- function(sample, priors, true_contexts) {

  lapply(priors, function(f) {
    fit_single_model(sample = sample, f = f, true_contexts = true_contexts)
  })
}

#' @description
#' Executes the simulation study using a pre-generated list of samples and
#' multiple context-tree functions, returning the results as a data frame.
#'
#' @details
#' For each sample in `sample_list`, corresponding to different sample sizes,
#' the Bayesian context tree model is fitted under each context-tree function
#' provided in `priors`. For every combination of sample size and prior,
#' the function computes the evaluation metrics.
#'
#' @param sample_list A list of samples generated from `generate_samples` function.
#' @param sample_sizes A numeric vector of sample sizes (must match `sample_list`).
#' @param priors A list of context-tree functions.
#' @param true_contexts A vector specifying the contexts of the true tree.
#'
#' @returns A data frame with one row per (sample size, prior), containing the evaluation
#' metrics for each fitted model.
#'
run_scenario_from_samples <- function(sample_list, sample_sizes, priors, true_contexts) {

  simulation_results <- lapply(sample_list, function(sample) {
    fit_priors_on_sample(
      sample = sample,
      priors = priors,
      true_contexts = true_contexts
    )
  })

  purrr::map2_dfr(
    simulation_results,
    sample_sizes,
    function(res_n, n_val) {
      purrr::imap_dfr(res_n, function(res_prior, prior_name) {
        tibble::tibble(
          n = n_val,
          prior = prior_name,
          delta_map = res_prior$delta_map,
          prior_true_tree = res_prior$prior_true_tree,
          post_true_tree = res_prior$post_true_tree,
          log_evidence = res_prior$log_10_evidence
        )
      })
    }
  )
}


