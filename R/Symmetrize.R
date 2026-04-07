#' Symmetrize a matrix or tensor
#'
#' For a matrix, computes (A + A^T) / 2.
#' For a tensor, averages over all permutations of indices.
#'
#' @param x A matrix or rTensor Tensor object (must have equal modes)
#' @return Symmetrized matrix or Tensor
#' @export
#' @examples
#' A <- matrix(c(1, 2, 3, 4), 2, 2)
#' Symmetrize(A)
#'
#' library(rTensor)
#' arr <- array(rnorm(27), dim = c(3, 3, 3))
#' S <- Symmetrize(as.tensor(arr))
#' isSymmetricTensor(S)  # TRUE
Symmetrize <- function(x) {
  if (is.matrix(x)) {
    if (nrow(x) != ncol(x)) stop("Matrix must be square")
    return((x + t(x)) / 2)
  }
  if (inherits(x, "Tensor")) {
    modes <- x@modes
    if (length(unique(modes)) != 1) {
      stop("All modes must be equal for symmetrization")
    }
    d <- x@num_modes
    if (d == 1) return(x)
    arr <- x@data
    perms <- .allPermutations(d)
    result <- array(0, dim = modes)
    for (p in perms) {
      result <- result + aperm(arr, p)
    }
    result <- result / length(perms)
    return(as.tensor(result))
  }
  stop("x must be a matrix or rTensor Tensor object")
}
