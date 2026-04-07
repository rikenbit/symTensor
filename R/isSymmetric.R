#' Check if a matrix or tensor is symmetric
#'
#' For a matrix, checks if A[i,j] == A[j,i] for all i,j.
#' For a tensor, checks if the tensor is invariant under all permutations
#' of its indices (i.e., fully symmetric).
#'
#' @param x A matrix or rTensor Tensor object
#' @param tol Numeric tolerance for comparison (default: 1e-10)
#' @return Logical indicating whether x is symmetric
#' @export
#' @examples
#' # Matrix
#' A <- matrix(c(1,2,2,3), 2, 2)
#' isSymmetricTensor(A)  # TRUE
#'
#' # 3-mode symmetric tensor
#' library(rTensor)
#' arr <- array(0, dim = c(3, 3, 3))
#' for (i in 1:3) for (j in 1:3) for (k in 1:3) {
#'   arr[i, j, k] <- i + j + k
#' }
#' isSymmetricTensor(as.tensor(arr))  # TRUE
isSymmetricTensor <- function(x, tol = 1e-10) {
  if (is.matrix(x)) {
    if (nrow(x) != ncol(x)) return(FALSE)
    return(all(abs(x - t(x)) <= tol))
  }
  if (inherits(x, "Tensor")) {
    modes <- x@modes
    # All modes must be equal for full symmetry
    if (length(unique(modes)) != 1) return(FALSE)
    d <- x@num_modes
    if (d == 1) return(TRUE)
    arr <- x@data
    # Check all pairwise transpositions
    perms <- .allPermutations(d)
    for (p in perms) {
      if (identical(p, seq_len(d))) next
      if (any(abs(arr - aperm(arr, p)) > tol)) return(FALSE)
    }
    return(TRUE)
  }
  stop("x must be a matrix or rTensor Tensor object")
}

#' Generate all permutations of 1:n
#' @param n Integer
#' @return List of integer vectors
#' @keywords internal
.allPermutations <- function(n) {
  if (n == 1) return(list(1L))
  if (n == 2) return(list(1:2, 2:1))
  result <- list()
  sub <- .allPermutations(n - 1L)
  for (perm in sub) {
    for (pos in seq_len(n)) {
      new_perm <- append(perm, n, after = pos - 1)
      result <- c(result, list(new_perm))
    }
  }
  result
}
