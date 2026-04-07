#' Label Propagation via PageRank
#'
#' Semi-supervised label propagation on a graph. Uses Personalized PageRank
#' internally to propagate labels from seed nodes to the rest of the graph.
#'
#' @param A Square non-negative adjacency/weight matrix (n x n)
#' @param seeds Named integer vector or factor. Names correspond to node indices
#'   (1-based), values are the class labels. Unlabeled nodes will be classified.
#' @param damping Damping factor passed to \code{\link{PageRank}} (default: 0.85)
#' @param max_iter Maximum iterations for PageRank (default: 1000)
#' @param tol Convergence tolerance for PageRank (default: 1e-8)
#' @return A list with components:
#'   \describe{
#'     \item{labels}{Integer vector of predicted labels for all nodes}
#'     \item{scores}{Matrix (n x K) of PageRank-based scores per class}
#'     \item{classes}{Unique class labels}
#'   }
#' @export
#' @examples
#' # 6-node graph with 2 communities
#' A <- matrix(0, 6, 6)
#' A[1,2] <- A[2,1] <- 1; A[1,3] <- A[3,1] <- 1; A[2,3] <- A[3,2] <- 1
#' A[4,5] <- A[5,4] <- 1; A[4,6] <- A[6,4] <- 1; A[5,6] <- A[6,5] <- 1
#' A[3,4] <- A[4,3] <- 0.5  # weak bridge
#' seeds <- c("1" = 1, "6" = 2)  # node 1 -> class 1, node 6 -> class 2
#' result <- LabelPropagation(A, seeds)
#' result$labels
LabelPropagation <- function(A, seeds, damping = 0.85,
                             max_iter = 1000L, tol = 1e-8) {
  if (!is.matrix(A)) stop("A must be a matrix")
  n <- nrow(A)

  seed_indices <- as.integer(names(seeds))
  seed_labels <- as.integer(seeds)
  classes <- sort(unique(seed_labels))
  K <- length(classes)

  # Run Personalized PageRank for each class
  scores <- matrix(0, nrow = n, ncol = K)
  colnames(scores) <- as.character(classes)

  for (k in seq_len(K)) {
    # Personalization: uniform over seed nodes of class k
    v <- rep(0, n)
    class_seeds <- seed_indices[seed_labels == classes[k]]
    v[class_seeds] <- 1 / length(class_seeds)

    pr <- PageRank(A, damping = damping, personalization = v,
                   max_iter = max_iter, tol = tol)
    scores[, k] <- pr$vector
  }

  # Assign each node to the class with highest score
  labels <- classes[apply(scores, 1, which.max)]

  list(labels = labels, scores = scores, classes = classes)
}
