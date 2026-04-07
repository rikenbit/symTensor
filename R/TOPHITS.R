#' TOPHITS: Tensor HITS Algorithm
#'
#' Computes hub and authority scores for higher-order link analysis using
#' Tucker decomposition of a tensor. Generalization of HITS (Hyperlink-Induced
#' Topic Search) to tensors.
#'
#' For a 3-mode tensor \eqn{T \in R^{I \times J \times K}},
#' TOPHITS computes Tucker decomposition \eqn{T \approx C \times_1 U_1 \times_2 U_2 \times_3 U_3}
#' and derives hub/authority scores from the factor matrices.
#'
#' For a symmetric tensor (all modes equal), a single set of scores is returned.
#'
#' @param x An rTensor Tensor object (2-mode or 3-mode)
#' @param R Integer or integer vector specifying the number of components per mode.
#'   If scalar, same rank is used for all modes.
#' @param max_iter Maximum iterations for Tucker decomposition (default: 100)
#' @param tol Convergence tolerance (default: 1e-8)
#' @return A list with components:
#'   \describe{
#'     \item{scores}{List of score matrices, one per mode (each n_i x R_i).
#'       For symmetric tensors, a single matrix.}
#'     \item{core}{Core tensor from Tucker decomposition}
#'     \item{rankings}{List of integer vectors giving node rankings per mode
#'       (ordered by dominant score)}
#'     \item{decomp}{Full Tucker decomposition result from rTensor}
#'   }
#' @references
#' Kolda, T. G., & Bader, B. W. (2006).
#' The TOPHITS Model for Higher-Order Web Link Analysis.
#' Workshop on Link Analysis, Counterterrorism and Security.
#' @export
#' @examples
#' library(rTensor)
#' # Create a 3-mode tensor
#' arr <- array(abs(rnorm(125)), dim = c(5, 5, 5))
#' tnsr <- as.tensor(arr)
#' result <- TOPHITS(tnsr, R = 2)
#' result$rankings
TOPHITS <- function(x, R = 2L, max_iter = 100L, tol = 1e-8) {
  if (!inherits(x, "Tensor")) {
    stop("x must be an rTensor Tensor object")
  }

  d <- x@num_modes
  modes <- x@modes

  # Set ranks per mode
  if (length(R) == 1) {
    ranks <- rep(R, d)
  } else {
    if (length(R) != d) stop("R must be scalar or length equal to num_modes")
    ranks <- R
  }
  # Ensure ranks don't exceed mode dimensions
  ranks <- pmin(ranks, modes)

  # Tucker decomposition (HOOI via rTensor)
  tuck <- rTensor::tucker(x, ranks = ranks, max_iter = max_iter, tol = tol)

  # Extract factor matrices (scores)
  U_list <- tuck$U
  core <- tuck$Z

  # Compute scores: |U_i| weighted by core tensor norms along each mode
  # For each mode, the score of node j is the norm of the j-th row of U_i
  # weighted by corresponding core slices
  scores <- list()
  rankings <- list()

  is_symmetric <- length(unique(modes)) == 1

  for (m in seq_len(d)) {
    U <- U_list[[m]]
    # Primary score: magnitude of each row in factor matrix
    # Weight by squared singular values (core tensor energy along this mode)
    core_unfold <- rTensor::k_unfold(core, m)@data
    weights <- rowSums(core_unfold^2)
    weights <- weights / sum(weights)

    # Weighted score for each node
    score_mat <- sweep(U^2, 2, weights, "*")
    node_scores <- rowSums(score_mat)

    scores[[m]] <- U
    rankings[[m]] <- order(node_scores, decreasing = TRUE)
  }

  if (is_symmetric) {
    names(scores) <- paste0("mode", seq_len(d))
    names(rankings) <- paste0("mode", seq_len(d))
  } else {
    names(scores) <- paste0("mode", seq_len(d))
    names(rankings) <- paste0("mode", seq_len(d))
  }

  list(
    scores = scores,
    core = core,
    rankings = rankings,
    decomp = tuck
  )
}
