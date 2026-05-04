# ------------------------------------------------------------
# Context-tree functions
# ------------------------------------------------------------
# This script implements the context-tree functions presented
# in Table 1 of the paper. Each function is defined as a function
# on nodes, f(s), corresponding to the third column of Table 1.
# ------------------------------------------------------------

#' @description
#' Unity function
#'
u <- function(node) {1}

#' @description
#' l-depth indicator function
#' @param l depth of interest.
#'
d_l <- function(l) {
  function(node) {
    if (node$getDepth() <= l) 1 else 0
  }
}

#' @description
#' a-renewal indicator function:
#' @param a a symbol from alphabet.
#'
i_a <- function(a) {
  function(node) {
    node <- node$getPath()
    symbols <- unlist(strsplit(node, "\\."))
    if (length(symbols) > 1 && any(symbols[-length(symbols)] == a)) {
      return(0)}
    else {return(1)}
    }
}

#' @description
#' Exponential function
#' @param alpha a positive real number.
#'
e_alpha <- function(alpha) {function(node) {exp(-alpha)}}

#' @description
#' l(s)-exponential function
#'
e_ls <- function(node) {exp(-node$getDepth())}

#' @description
#' CTW function
#' @param maximalDepth The maximal depth `L` considered in the model.
#'
ctw <- function(maximalDepth) {
  function(node) {
    if (node$getDepth() < maximalDepth) {return(1/4)
    } else {return(1/2)}
  }
}

#' @description
#' BCT function
#' @param beta Stopping probability of the nodes
#' @param maximalDepth The maximal depth `L` considered in the model
#' @param m Length of alphabet `A` considered.
#'
b_beta <- function(beta, maximalDepth, m) {
  function(node) {
    if (node$getDepth() < maximalDepth) {
      return((1 - beta)^(1 / (m - 1)) * beta)
    } else {return((1 - beta)^(1 / (m - 1)))}
  }
}

#' @description
#' Target l-depth function
#' @param l Depth of interest
#' @param k Concentration parameter
#'
tcl <- function(l, k) {
  function(node) {
    return(k^(-abs(node$getDepth() - l)))
  }
}

