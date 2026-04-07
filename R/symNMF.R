#' Symmetric Non-negative Matrix Factorization (symNMF)
#'
#' Decomposes a symmetric non-negative matrix \eqn{\Omega \approx M S M^\top}
#' using multiplicative update (MU) rules with Beta-divergence.
#'
#' Two models are supported:
#' \itemize{
#'   \item \strong{MSM}: \eqn{\Omega \approx M S M^\top} where M is N x K and S is K x K
#'   \item \strong{MM}: \eqn{\Omega \approx M M^\top} where M is N x K (special case with S = I)
#' }
#'
#' The loss function is the Beta-divergence \eqn{D_\beta(\Omega \| \hat{\Omega})}:
#' \itemize{
#'   \item \eqn{\beta = 2}: Frobenius (Euclidean) distance
#'   \item \eqn{\beta = 1}: Kullback-Leibler divergence
#'   \item \eqn{\beta = 0}: Itakura-Saito divergence
#'   \item Other values of \eqn{\beta} interpolate/extrapolate these cases
#' }
#'
#' @param Omega Symmetric non-negative matrix (N x N)
#' @param K Number of components (rank)
#' @param model Character, either \code{"MSM"} (default) or \code{"MM"}
#' @param Beta Beta parameter for Beta-divergence (default: 2).
#'   2 = Frobenius, 1 = KL, 0 = IS.
#' @param initM Initial M matrix (N x K). If NULL, random non-negative initialization.
#' @param initS Initial S matrix (K x K). If NULL, random non-negative initialization.
#'   Ignored when \code{model = "MM"}.
#' @param lambda Regularization parameter for \eqn{|\det(S)|} penalty (default: 0).
#'   Only used with \code{model = "MSM"}.
#' @param num.iter Maximum number of iterations (default: 100)
#' @param thr Convergence threshold on relative change (default: 1e-10)
#' @param verbose Logical, print iteration info (default: FALSE)
#' @return A list with components:
#'   \describe{
#'     \item{M}{Factor matrix (N x K)}
#'     \item{S}{Coefficient matrix (K x K). Identity matrix when \code{model = "MM"}.}
#'     \item{RecError}{Vector of reconstruction errors at each iteration}
#'     \item{RelChange}{Vector of relative changes at each iteration}
#'     \item{converged}{Logical}
#'     \item{iter}{Number of iterations performed}
#'   }
#' @references
#' Kuang, D., Park, H., & Ding, C. H. Q. (2012).
#' Symmetric Nonnegative Matrix Factorization for Graph Clustering.
#' SDM 2012.
#'
#' Fevotte, C. & Idier, J. (2011).
#' Algorithms for Nonnegative Matrix Factorization with the
#' Beta-Divergence. Neural Computation, 23(9).
#' @export
#' @examples
#' set.seed(42)
#' M_true <- matrix(c(rep(c(1,0), 5), rep(c(0,1), 5)), 10, 2)
#' Omega <- M_true %*% t(M_true) + matrix(runif(100, 0, 0.1), 10, 10)
#' Omega <- (Omega + t(Omega)) / 2
#'
#' # Frobenius (Beta=2, default)
#' res_fro <- symNMF(Omega, K = 2, Beta = 2)
#'
#' # KL divergence (Beta=1)
#' res_kl <- symNMF(Omega, K = 2, Beta = 1)
#'
#' # Itakura-Saito (Beta=0)
#' res_is <- symNMF(Omega, K = 2, Beta = 0)
symNMF <- function(Omega, K, model = c("MSM", "MM"),
                   Beta = 2,
                   initM = NULL, initS = NULL,
                   lambda = 0,
                   num.iter = 100L, thr = 1e-10,
                   verbose = FALSE) {
  model <- match.arg(model)
  N <- nrow(Omega)
  eps <- .Machine$double.eps

  # Initialize M
  if (is.null(initM)) {
    M <- matrix(runif(N * K, 0, 1), N, K)
  } else {
    M <- initM
  }

  # Initialize S
  if (model == "MSM") {
    if (is.null(initS)) {
      S <- matrix(runif(K * K, 0, 1), K, K)
      S <- (S + t(S)) / 2
    } else {
      S <- initS
    }
  } else {
    S <- diag(K)
  }

  RecError <- numeric(num.iter)
  RelChange <- numeric(num.iter)
  converged <- FALSE

  for (iter in seq_len(num.iter)) {
    # Reconstruction
    Approx <- M %*% S %*% t(M)
    Approx_safe <- Approx + eps

    # Compute Beta-divergence
    RecError[iter] <- .betaDiv(Omega, Approx_safe, Beta)

    # Relative change
    if (iter > 1) {
      RelChange[iter] <- abs(RecError[iter - 1] - RecError[iter]) /
        (abs(RecError[iter]) + eps)
      if (RelChange[iter] < thr) {
        converged <- TRUE
        if (verbose) message("Converged at iteration ", iter)
        break
      }
    }

    if (verbose && (iter %% 10 == 0 || iter == 1)) {
      message("Iter ", iter, ": RecError = ", round(RecError[iter], 6))
    }

    # General Beta-divergence MU updates
    # Weight matrix: W = Omega^(Beta-2) * Approx^(-(Beta-2)) for numerator part
    # Numerator weight: Omega^(Beta-1) * Approx^(-(Beta-1)) = (Omega/Approx)^(Beta-1)
    # Denominator weight: Approx^(Beta-1) * Approx^(-(Beta-1)) = 1 (but with structure)
    #
    # General MU rule for V in X ≈ W H:
    #   H <- H * ((W^T (X * Approx^(Beta-2))) / (W^T Approx^(Beta-1)))^gamma
    # where gamma = 1 for Beta in [1,2], specific exponents otherwise
    #
    # For symmetric case Omega ≈ M S M^T:
    numer_weight <- Omega * Approx_safe^(Beta - 2)
    denom_weight <- Approx_safe^(Beta - 1)

    if (model == "MSM") {
      # Update M
      numer_M <- numer_weight %*% M %*% S
      denom_M <- denom_weight %*% M %*% S + eps
      gamma <- .muExponent(Beta)
      M <- M * (numer_M / denom_M)^gamma

      # Recompute Approx after M update
      Approx <- M %*% S %*% t(M)
      Approx_safe <- Approx + eps
      numer_weight <- Omega * Approx_safe^(Beta - 2)
      denom_weight <- Approx_safe^(Beta - 1)

      # Update S
      numer_S <- t(M) %*% numer_weight %*% M
      denom_S <- t(M) %*% denom_weight %*% M + eps
      if (lambda > 0) {
        S_inv <- tryCatch(solve(S + eps * diag(K)), error = function(e) NULL)
        if (!is.null(S_inv)) {
          denom_S <- denom_S + lambda * abs(det(S)) * t(S_inv)
        }
      }
      S <- S * (numer_S / denom_S)^gamma
    } else {
      # model == "MM": Omega ≈ M M^T
      numer_M <- numer_weight %*% M
      denom_M <- denom_weight %*% M + eps
      gamma <- .muExponent(Beta)
      M <- M * (numer_M / denom_M)^gamma
    }

    # Ensure non-negativity
    M[M < eps] <- eps
    if (model == "MSM") S[S < eps] <- eps
  }

  RecError <- RecError[seq_len(iter)]
  RelChange <- RelChange[seq_len(iter)]

  list(
    M = M,
    S = S,
    RecError = RecError,
    RelChange = RelChange,
    converged = converged,
    iter = iter
  )
}

#' Compute Beta-divergence D_beta(X || Y)
#'
#' @param X Target matrix
#' @param Y Approximation matrix
#' @param beta Beta parameter
#' @return Scalar divergence value
#' @keywords internal
.betaDiv <- function(X, Y, beta) {
  eps <- .Machine$double.eps
  if (abs(beta - 2) < eps) {
    # Frobenius (Euclidean)
    sum((X - Y)^2) / 2
  } else if (abs(beta - 1) < eps) {
    # KL divergence
    sum(X * log(X / Y + eps) - X + Y)
  } else if (abs(beta) < eps) {
    # Itakura-Saito
    sum(X / Y - log(X / Y + eps) - 1)
  } else {
    # General Beta
    sum(X^beta / (beta * (beta - 1)) +
        Y^beta / beta -
        X * Y^(beta - 1) / (beta - 1))
  }
}

#' MU exponent for Beta-divergence
#'
#' Returns the exponent gamma for the multiplicative update rule.
#' gamma = 1 for Beta in [1,2], otherwise adjusted for convergence guarantee.
#'
#' @param beta Beta parameter
#' @return Exponent gamma
#' @references
#' Fevotte, C. & Idier, J. (2011). Neural Computation, 23(9).
#' @keywords internal
.muExponent <- function(beta) {
  if (beta >= 1 && beta <= 2) {
    1
  } else {
    1 / (3 - beta)
  }
}
