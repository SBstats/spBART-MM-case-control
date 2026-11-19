## BART: Bayesian Additive Regression Trees - Stateful Sampler
## Modified for the semi-parametric probit BART model (spBART)
## Version 2.9.9.1: probit error s.d. fixed at 1, RNG synchronised with R,
## bartModelMatrix() cutpoints passed to C++, leaf prior default matches pbart(),
## DART switched on explicitly via wbart_start_dart() (pbart() uses burn/2).
##
## This file provides R wrapper functions for the stateful BART sampler
## that can be used within Gibbs samplers for semi-parametric models.

#' Create Stateful BART Sampler
#'
#' Creates a stateful BART sampler object that maintains tree state across iterations.
#' This is designed for use within Gibbs samplers for semi-parametric models.
#'
#' @param x.train Training predictor matrix (n x p)
#' @param y.train Training response vector (length n)
#' @param x.test Optional test predictor matrix (n_test x p)
#' @param sparse Logical; if TRUE, use Dirichlet sparsity prior (DART)
#' @param theta DART parameter (not used when sparse=TRUE with Dirichlet)
#' @param omega DART parameter (not used when sparse=TRUE with Dirichlet)
#' @param a Dirichlet prior parameter (sparse splitting rule a_tau)
#' @param b Dirichlet prior parameter (sparse splitting rule b_tau)
#' @param rho Dirichlet concentration parameter (default: p)
#' @param augment Logical; use data augmentation for DART
#' @param xinfo Matrix of cutpoints (if pre-specified)
#' @param usequants Logical; if TRUE, use quantiles for cutpoints
#' @param numcut Number of cutpoints for each predictor
#' @param cont Logical vector indicating continuous variables
#' @param rm.const Logical; if TRUE, remove constant predictors
#' @param power BART tree prior parameter (beta in Chipman et al. 2010)
#' @param base BART tree prior parameter (alpha in Chipman et al. 2010)
#' @param k Leaf prior scale; with tau = NULL, tau = 3/(k*sqrt(ntree)), the pbart() default
#' @param tau Prior standard deviation for leaf parameters (default: 3/(k*sqrt(ntree)))
#' @param ntree Number of trees in ensemble
#' @param transposed Logical; if TRUE, x.train is already transposed (p x n)
#'
#' @return An external pointer to the stateful BART sampler object
#'
#' @details
#' This function creates a BART sampler that can be updated iteratively within
#' a Gibbs sampler. Unlike wbart() which runs a complete MCMC chain, this:
#' 1. Initializes trees ONCE
#' 2. Allows updating the response without reinitializing (wbart_update_response)
#' 3. Runs one iteration at a time (wbart_run_iteration)
#' 4. Must be explicitly destroyed when done (wbart_destroy)
#'
#' @examples
#' \dontrun{
#' # Create sampler
#' sampler <- wbart_create(X, y, sparse=TRUE, rho=ncol(X))
#'
#' # Update response and run iteration
#' wbart_update_response(sampler, y_new)
#' result <- wbart_run_iteration(sampler)
#'
#' # Clean up
#' wbart_destroy(sampler)
#' }
#'
#' @export
wbart_create <- function(
  x.train, y.train, x.test = matrix(0.0, 0, 0),
  sparse = FALSE, theta = 0, omega = 1,
  a = 0.5, b = 1, rho = NULL, augment = FALSE,
  xinfo = matrix(0.0, 0, 0), usequants = FALSE,
  numcut = 100L, cont = FALSE, rm.const = TRUE,
  power = 2.0, base = 0.95, k = 2.0, tau = NULL,
  ntree = 200L, transposed = FALSE
) {

  #--------------------------------------------------------------------------
  # Process training data
  #--------------------------------------------------------------------------
  n <- length(y.train)

  if (!transposed) {
    temp <- bartModelMatrix(x.train, numcut, usequants = usequants,
                           cont = cont, xinfo = xinfo, rm.const = rm.const)
    x.train <- t(temp$X)
    numcut <- temp$numcut
    xinfo <- temp$xinfo

    if (length(x.test) > 0) {
      x.test <- bartModelMatrix(x.test)
      x.test <- t(x.test[, temp$rm.const])
    }

    rm.const <- temp$rm.const
    rm(temp)
  } else {
    rm.const <- NULL
  }

  if (n != ncol(x.train)) {
    stop('The length of y.train and the number of rows in x.train must be identical')
  }

  p <- nrow(x.train)
  np <- ncol(x.test)

  # Set default rho
  if (is.null(rho)) rho <- p
  if (length(rm.const) == 0) rm.const <- 1:p

  #--------------------------------------------------------------------------
  # Set default tau if not provided
  #--------------------------------------------------------------------------
  if (is.null(tau)) {
    # Same leaf prior as pbart(): sigma_mu = 3 / (k * sqrt(ntree)), k = 2
    # (Chipman, George and McCulloch 2010, binary-response default)
    tau <- 3 / (k * sqrt(ntree))
  }

  #--------------------------------------------------------------------------
  # Call C++ function to create sampler
  #--------------------------------------------------------------------------
  sampler <- .Call("cwbart_create",
                  as.integer(n),        # n
                  as.integer(p),        # p
                  as.integer(np),       # np
                  x.train,              # ix
                  y.train,              # iy
                  x.test,               # ixp
                  as.integer(ntree),    # m
                  as.integer(numcut),   # numcut
                  power,                # mybeta
                  base,                 # alpha
                  tau,                  # tau
                  as.integer(sparse),   # dart
                  theta,                # theta
                  omega,                # omega
                  a,                    # a
                  b,                    # b
                  rho,                  # rho
                  as.integer(augment),  # aug
                  xinfo,                # Xinfo
                  PACKAGE = "BART")

  # Store metadata as attributes
  attr(sampler, "n") <- n
  attr(sampler, "p") <- p
  attr(sampler, "np") <- np
  attr(sampler, "ntree") <- ntree
  attr(sampler, "rm.const") <- rm.const

  return(sampler)
}

#' Update Response in Stateful BART Sampler
#'
#' Updates the response vector in the stateful BART sampler without
#' reinitializing the trees.
#'
#' @param sampler External pointer from wbart_create()
#' @param y.train.new New response vector (same length as original)
#'
#' @details
#' This function updates the response vector that the BART sampler uses,
#' but keeps all tree structures intact. This is essential for Gibbs sampling
#' where the response changes each iteration (e.g., U - Z'beta in semi-parametric models).
#'
#' @export
wbart_update_response <- function(sampler, y.train.new) {
  .Call("cwbart_update_response",
        sampler,
        as.numeric(y.train.new),
        PACKAGE = "BART")

  invisible(NULL)
}

#' Run One Iteration of Stateful BART Sampler
#'
#' Runs one iteration of the BART MCMC sampler from the current tree state.
#'
#' @param sampler External pointer from wbart_create()
#'
#' @return A list with components:
#' \item{yhat.train}{Fitted values for training data (sum of all trees)}
#' \item{yhat.test}{Fitted values for test data (if test data was provided)}
#' \item{varcount}{Integer vector of length p indicating how many times each variable was used as a splitting rule}
#' \item{sigma}{Error standard deviation (fixed at 1, probit latent scale)}
#' \item{dart_on}{Whether the Dirichlet update is switched on}
#'
#' @details
#' This function runs ONE iteration of the BART MCMC sampler:
#' 1. For each tree, propose birth/death/change/swap moves
#' 2. Update tree structure and leaf parameters
#' 3. If DART has been switched on (wbart_start_dart), update splitting
#'    probabilities from their Dirichlet full conditional
#' 4. Return current state (fitted values and variable usage)
#'
#' The trees maintain their structure across calls, allowing proper MCMC mixing.
#'
#' @export
wbart_run_iteration <- function(sampler) {
  result <- .Call("cwbart_run_iteration",
                  sampler,
                  PACKAGE = "BART")

  return(result)
}

#' Switch On the Dirichlet (DART) Update
#'
#' Turns on the Dirichlet splitting-probability update for a sampler created
#' with sparse = TRUE. Calling it more than once has no further effect.
#' pbart() switches DART on at iteration burn/2; call this at the same point
#' to use the same schedule.
#'
#' @param sampler External pointer from wbart_create()
#' @return Logical: TRUE if DART is now on
#'
#' @export
wbart_start_dart <- function(sampler) {
  .Call("cwbart_start_dart", sampler, PACKAGE = "BART")
}

#' Destroy Stateful BART Sampler
#'
#' Frees memory associated with the stateful BART sampler.
#'
#' @param sampler External pointer from wbart_create()
#'
#' @details
#' This function MUST be called when you are done with the sampler to avoid
#' memory leaks. It deletes all C++ objects and frees allocated memory.
#'
#' @export
wbart_destroy <- function(sampler) {
  .Call("cwbart_destroy",
        sampler,
        PACKAGE = "BART")

  invisible(NULL)
}
