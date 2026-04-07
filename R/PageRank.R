#' (Personalized) PageRank
#'
#' Computes the PageRank vector for a given adjacency or transition matrix
#' using power iteration.
#'
#' @param A Square non-negative matrix (adjacency or weight matrix).
#'   Columns are normalized internally to form a column-stochastic matrix.
#' @param damping Damping factor (default: 0.85). Probability of following a link.
#' @param personalization Optional personalization vector of length \code{nrow(A)}.
#'   If NULL, uniform distribution is used (standard PageRank).
#'   Personalized PageRank biases the random jumps toward this distribution.
#' @param max_iter Maximum number of power iterations (default: 1000)
#' @param tol Convergence tolerance on L1 norm change (default: 1e-8)
#' @return A list with components:
#'   \describe{
#'     \item{vector}{PageRank vector (sums to 1)}
#'     \item{iter}{Number of iterations until convergence}
#'     \item{converged}{Logical, whether the algorithm converged}
#'   }
#' @export
#' @examples
#' # Simple web graph
#' A <- matrix(c(0,1,1,0, 1,0,0,1, 0,1,0,1, 1,0,1,0), 4, 4)
#' pr <- PageRank(A)
#' pr$vector
PageRank <- function(A, damping = 0.85, personalization = NULL,
                     max_iter = 1000L, tol = 1e-8) {
  if (!is.matrix(A)) stop("A must be a matrix")
  n <- nrow(A)
  if (n != ncol(A)) stop("A must be square")

  # Personalization vector
  if (is.null(personalization)) {
    v <- rep(1 / n, n)
  } else {
    if (length(personalization) != n) {
      stop("personalization must have length nrow(A)")
    }
    v <- personalization / sum(personalization)
  }

  # Column-normalize A to get transition matrix
  col_sums <- colSums(A)
  # Handle dangling nodes (columns with zero sum)
  dangling <- which(col_sums == 0)
  col_sums[dangling] <- 1  # avoid division by zero
  M <- sweep(A, 2, col_sums, "/")

  # Power iteration: pi = d * M %*% pi + (1-d) * v + d * dangling_contribution
  pi <- rep(1 / n, n)
  converged <- FALSE

  for (iter in seq_len(max_iter)) {
    pi_new <- damping * (M %*% pi) + (1 - damping) * v
    # Dangling node contribution: redistribute their mass
    if (length(dangling) > 0) {
      dangling_mass <- damping * sum(pi[dangling])
      pi_new <- pi_new + dangling_mass * v
    }
    pi_new <- as.vector(pi_new)
    pi_new <- pi_new / sum(pi_new)

    if (sum(abs(pi_new - pi)) < tol) {
      converged <- TRUE
      pi <- pi_new
      break
    }
    pi <- pi_new
  }

  list(vector = pi, iter = iter, converged = converged)
}
