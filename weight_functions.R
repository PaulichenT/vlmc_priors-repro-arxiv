# ------------------------------------------------------------
# Weight functions
# ------------------------------------------------------------
# This script implements the weight functions presented
# in Table 1 of the paper. 
# ------------------------------------------------------------

#' @description
#' Unity function
#'
u <- function(node) {1}


#' @description
#' beta-constant function
#' @param beta a constant 
#'
k_beta <- function(beta) {function(node) {beta}}


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
#' l-m-depth indicator function
#' @param l a depth (integer)
#' @param m a depth (integer)
#'
d_l_m <- function(l, m) {
  function(node) {
    if (node$getDepth() >= l && node$getDepth() <= m) 1 else 0
  }
}


#' @description
#' beta-exponential function
#' @param beta a constant
#'
e_beta <- function(beta) {
  function(node) {exp(beta * node$getDepth())}}


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
bct_beta <- function(beta, maximalDepth, m) {
  function(node) {
    if (node$getDepth() < maximalDepth) {
      return((1 - beta)^(1 / (m - 1)) * beta)
    } else {return((1 - beta)^(1 / (m - 1)))}
  }
}


#' @description
#' beta-target l-depth function
#' @param l Depth of interest
#' @param beta Concentration parameter
#'
t_beta_l <- function(l, beta) {
  function(node) {
    return(beta^(-abs(node$getDepth() - l)))
  }
}
