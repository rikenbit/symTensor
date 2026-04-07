#' Higher-order Power Method for Symmetric Tensor Decomposition
#'
#' Decomposes a symmetric tensor into a sum of rank-1 symmetric terms
#' using the higher-order power method (S-HOPM) with deflation.
#'
#' For a d-th order symmetric tensor \eqn{T \in R^{n \times n \times \cdots \times n}},
#' finds \eqn{T \approx \sum_{r=1}^{R} \lambda_r v_r \otimes v_r \otimes \cdots \otimes v_r}.
#'
#' @param x A symmetric rTensor Tensor object (all modes must be equal)
#' @param R Number of rank-1 components to extract (default: 1)
#' @param max_iter Maximum iterations per component (default: 1000)
#' @param tol Convergence tolerance (default: 1e-8)
#' @param init Initialization method: \code{"random"} (default) or \code{"svd"}
#' @return A list with components:
#'   \describe{
#'     \item{lambdas}{Numeric vector of length R with eigenvalues}
#'     \item{vectors}{Matrix (n x R) where each column is a unit eigenvector}
#'     \item{residual}{Residual tensor after deflation}
#'     \item{converged}{Logical vector of length R}
#'     \item{iters}{Integer vector of iterations per component}
#'   }
#' @references
#' Kolda, T. G., & Mayo, J. R. (2011).
#' Shifted Power Method for Computing Tensor Eigenpairs.
#' SIAM Journal on Matrix Analysis and Applications.
#'
#' De Lathauwer, L., De Moor, B., & Vandewalle, J. (2000).
#' On the Best Rank-1 and Rank-(R1, R2, ..., RN) Approximation of
#' Higher-Order Tensors. SIAM Journal on Matrix Analysis and Applications.
#' @export
#' @examples
#' library(rTensor)
#' # Create a symmetric 3-mode tensor
#' n <- 5
#' v1 <- rnorm(n); v1 <- v1 / sqrt(sum(v1^2))
#' # T = 3.0 * v1 o v1 o v1
#' arr <- 3.0 * outer(outer(v1, v1), v1)
#' tnsr <- as.tensor(arr)
#' result <- HigherOrderPower(tnsr, R = 1)
#' result$lambdas  # should be close to 3.0
HigherOrderPower <- function(x, R = 1L, max_iter = 1000L, tol = 1e-8,
                             init = c("random", "svd")) {
  init <- match.arg(init)

  if (!inherits(x, "Tensor")) {
    stop("x must be an rTensor Tensor object")
  }
  modes <- x@modes
  if (length(unique(modes)) != 1) {
    stop("All modes must be equal (symmetric tensor required)")
  }
  n <- modes[1]
  d <- x@num_modes

  lambdas <- numeric(R)
  vectors <- matrix(0, n, R)
  converged_vec <- logical(R)
  iters_vec <- integer(R)

  residual <- x@data

  for (r in seq_len(R)) {
    # Initialize
    if (init == "svd" && d >= 2) {
      # Use leading left singular vector of mode-1 unfolding
      unf <- matrix(residual, n, n^(d - 1))
      v <- svd(unf, nu = 1, nv = 0)$u[, 1]
    } else {
      v <- rnorm(n)
    }
    v <- v / sqrt(sum(v^2))

    conv <- FALSE
    for (iter in seq_len(max_iter)) {
      # Contract residual tensor with v along all modes except the first
      # T(I, v, v, ..., v)
      w <- .tensorVectorContract(residual, v, d)
      lambda <- sum(w * v)
      w_norm <- sqrt(sum(w^2))

      if (w_norm < .Machine$double.eps) {
        v_new <- v
      } else {
        v_new <- w / w_norm
      }

      # Check convergence (allow sign flip)
      change <- min(sqrt(sum((v_new - v)^2)), sqrt(sum((v_new + v)^2)))
      if (change < tol) {
        conv <- TRUE
        v <- v_new
        break
      }
      v <- v_new
    }

    # Final eigenvalue
    lambda <- .tensorVectorContractAll(residual, v, d)

    lambdas[r] <- lambda
    vectors[, r] <- v
    converged_vec[r] <- conv
    iters_vec[r] <- iter

    # Deflation: subtract rank-1 component
    rank1 <- lambda * .outerRepeat(v, d)
    residual <- residual - rank1
  }

  list(
    lambdas = lambdas,
    vectors = vectors,
    residual = as.tensor(residual),
    converged = converged_vec,
    iters = iters_vec
  )
}

#' Contract a tensor with a vector along all modes except the first
#' @param arr Array of dimension (n, n, ..., n)
#' @param v Numeric vector of length n
#' @param d Number of modes
#' @return Numeric vector of length n
#' @keywords internal
.tensorVectorContract <- function(arr, v, d) {
  # Contract modes 2, 3, ..., d with v
  result <- arr
  for (m in seq(d, 2, by = -1)) {
    # Contract mode m with v
    result <- .contractMode(result, v, m)
  }
  as.vector(result)
}

#' Full contraction of tensor with vector along all modes
#' @keywords internal
.tensorVectorContractAll <- function(arr, v, d) {
  result <- arr
  for (m in seq(d, 1, by = -1)) {
    result <- .contractMode(result, v, m)
  }
  as.numeric(result)
}

#' Contract an array along a specific mode with a vector
#' @param arr Multi-dimensional array
#' @param v Vector
#' @param mode Mode to contract
#' @return Array with one fewer dimension
#' @keywords internal
.contractMode <- function(arr, v, mode) {
  dims <- dim(arr)
  if (is.null(dims)) {
    # Already a vector
    return(sum(arr * v))
  }
  d <- length(dims)
  n <- dims[mode]

  # Permute so that the target mode is the last
  perm <- c(seq_len(d)[-mode], mode)
  arr_perm <- aperm(arr, perm)

  # Reshape to matrix: (product of other dims) x n
  other_dims <- dims[-mode]
  arr_mat <- matrix(arr_perm, prod(other_dims), n)

  # Multiply
  result <- arr_mat %*% v

  # Reshape back
  if (length(other_dims) == 1) {
    return(as.vector(result))
  }
  array(result, dim = other_dims)
}

#' Compute repeated outer product of a vector
#' @param v Numeric vector
#' @param d Number of times to take outer product
#' @return d-dimensional array
#' @keywords internal
.outerRepeat <- function(v, d) {
  n <- length(v)
  if (d == 1) return(v)
  result <- outer(v, v)
  if (d == 2) return(result)
  for (i in 3:d) {
    result <- outer(result, v)
  }
  result
}
