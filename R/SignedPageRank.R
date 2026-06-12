#' Signed PageRank
#'
#' Computes PageRank on a signed graph using the doubled state-space approach.
#' The key property is that a negative edge followed by a negative edge
#' produces a positive contribution ("enemy of my enemy is my friend").
#'
#' The algorithm constructs an extended transition matrix:
#' \deqn{M = \begin{pmatrix} P^+ & P^- \\ P^- & P^+ \end{pmatrix}}
#' where \eqn{P^+} and \eqn{P^-} are column-normalized positive and negative
#' transition matrices, respectively. The state vector is doubled into
#' positive and negative components, and standard power iteration is applied.
#'
#' @param A_pos Square non-negative matrix of positive edge weights (n x n)
#' @param A_neg Square non-negative matrix of negative edge weights (n x n)
#' @param damping Damping factor (default: 0.85)
#' @param personalization Optional personalization vector of length n.
#'   If NULL, uniform distribution is used.
#' @param max_iter Maximum number of power iterations (default: 1000)
#' @param tol Convergence tolerance on L1 norm change (default: 1e-8)
#' @return A list with components:
#'   \describe{
#'     \item{positive}{Positive PageRank scores (length n)}
#'     \item{negative}{Negative PageRank scores (length n)}
#'     \item{net}{Net scores: positive - negative (length n)}
#'     \item{total}{Total scores: positive + negative (length n)}
#'     \item{iter}{Number of iterations until convergence}
#'     \item{converged}{Logical, whether the algorithm converged}
#'   }
#' @export
#' @examples
#' # A --(-)--> B --(-)--> C: net effect A->C should be positive
#' A_pos <- matrix(0, 3, 3)
#' A_neg <- matrix(0, 3, 3)
#' A_neg[2, 1] <- 1  # A inhibits B
#' A_neg[3, 2] <- 1  # B inhibits C
#' result <- SignedPageRank(A_pos, A_neg)
#' result$net
SignedPageRank <- function(A_pos, A_neg, damping = 0.85,
                           personalization = NULL,
                           max_iter = 1000L, tol = 1e-8) {
  if (!is.matrix(A_pos) || !is.matrix(A_neg))
    stop("A_pos and A_neg must be matrices")
  n <- nrow(A_pos)
  if (n != ncol(A_pos) || nrow(A_neg) != n || ncol(A_neg) != n)
    stop("A_pos and A_neg must be square matrices of the same size")
  if (any(A_pos < 0) || any(A_neg < 0))
    stop("A_pos and A_neg must be non-negative")

  # Out-degree normalization: D_i = sum_j (A_pos[j,i] + A_neg[j,i])
  out_deg <- colSums(A_pos) + colSums(A_neg)
  dangling <- which(out_deg == 0)
  out_deg[dangling] <- 1  # avoid division by zero

  P_pos <- sweep(A_pos, 2, out_deg, "/")
  P_neg <- sweep(A_neg, 2, out_deg, "/")

  # Extended transition matrix M (2n x 2n), stored as block operations
  # M = [P+ P-]
  #     [P- P+]

  # Personalization vector for extended space
  if (is.null(personalization)) {
    v <- rep(1 / n, n)
  } else {
    if (length(personalization) != n)
      stop("personalization must have length nrow(A_pos)")
    v <- personalization / sum(personalization)
  }
  # Start in positive state
  v_ext <- c(v, rep(0, n))

  # Initial state: uniform in positive half
  r <- c(rep(1 / n, n), rep(0, n))
  converged <- FALSE

  for (iter in seq_len(max_iter)) {
    r_pos <- r[1:n]
    r_neg <- r[(n + 1):(2 * n)]

    # M^T %*% r (column-stochastic, so we use M %*% r for row perspective)
    # r_pos_new = P+^T %*% r_pos + P-^T %*% r_neg
    # r_neg_new = P-^T %*% r_pos + P+^T %*% r_neg
    new_pos <- as.vector(P_pos %*% r_pos + P_neg %*% r_neg)
    new_neg <- as.vector(P_neg %*% r_pos + P_pos %*% r_neg)

    # Dangling node contribution
    if (length(dangling) > 0) {
      dangling_mass <- sum(r_pos[dangling]) + sum(r_neg[dangling])
      new_pos <- new_pos + dangling_mass * v
    }

    r_new <- c(new_pos, new_neg)
    r_new <- damping * r_new + (1 - damping) * v_ext

    # Normalize
    r_total <- sum(r_new)
    if (r_total > 0) r_new <- r_new / r_total

    if (sum(abs(r_new - r)) < tol) {
      converged <- TRUE
      r <- r_new
      break
    }
    r <- r_new
  }

  pos_scores <- r[1:n]
  neg_scores <- r[(n + 1):(2 * n)]

  list(
    positive = pos_scores,
    negative = neg_scores,
    net = pos_scores - neg_scores,
    total = pos_scores + neg_scores,
    iter = iter,
    converged = converged
  )
}
