#' Signed Path Contribution
#'
#' Extracts the contribution of paths (1 to \code{max_hop} hops) in a signed
#' graph. For each path, computes the sign (product of edge signs along the
#' path) and the contribution (product of transition probabilities).
#'
#' Key property: a path with an even number of negative edges has a net
#' positive effect (e.g., "enemy of my enemy is my friend").
#'
#' @param A_pos Square non-negative matrix of positive edge weights (n x n)
#' @param A_neg Square non-negative matrix of negative edge weights (n x n)
#' @param max_hop Maximum path length to enumerate (default: 3)
#' @param min_contribution Minimum contribution threshold to report (default: 0).
#'   Paths with contribution below this are omitted.
#' @return A data.frame with columns:
#'   \describe{
#'     \item{source}{Source node index}
#'     \item{target}{Target node index}
#'     \item{hop}{Number of hops}
#'     \item{path}{Character string of node indices (e.g., "1->3->2")}
#'     \item{signs}{Character string of edge signs along path (e.g., "-,-")}
#'     \item{net_sign}{Final sign of path: +1 or -1}
#'     \item{contribution}{Product of transition probabilities along path}
#'   }
#' @export
#' @examples
#' # A --(-)--> B --(-)--> C: path A->B->C has net positive sign
#' A_pos <- matrix(0, 3, 3)
#' A_neg <- matrix(0, 3, 3)
#' A_neg[2, 1] <- 1  # A -> B negative
#' A_neg[3, 2] <- 1  # B -> C negative
#' paths <- SignedPathContribution(A_pos, A_neg, max_hop = 2)
#' paths
SignedPathContribution <- function(A_pos, A_neg, max_hop = 3L,
                                   min_contribution = 0) {
  if (!is.matrix(A_pos) || !is.matrix(A_neg))
    stop("A_pos and A_neg must be matrices")
  n <- nrow(A_pos)
  if (n != ncol(A_pos) || nrow(A_neg) != n || ncol(A_neg) != n)
    stop("A_pos and A_neg must be square matrices of the same size")

  # Out-degree normalization
  out_deg <- colSums(A_pos) + colSums(A_neg)
  out_deg[out_deg == 0] <- 1
  P_pos <- sweep(A_pos, 2, out_deg, "/")
  P_neg <- sweep(A_neg, 2, out_deg, "/")

  # Use environment to accumulate results across recursive calls
  acc <- new.env(parent = emptyenv())
  acc$results <- list()

  for (hop in seq_len(max_hop)) {
    for (i in seq_len(n)) {
      .enumeratePaths(acc, i, i, c(), c(), 1, hop, n,
                      P_pos, P_neg, min_contribution)
    }
  }

  results <- acc$results

  if (length(results) == 0) {
    return(data.frame(
      source = integer(0), target = integer(0), hop = integer(0),
      path = character(0), signs = character(0),
      net_sign = integer(0), contribution = numeric(0),
      stringsAsFactors = FALSE
    ))
  }

  do.call(rbind, lapply(results, as.data.frame, stringsAsFactors = FALSE))
}

#' Recursively enumerate paths in signed graph
#' @param acc Environment with a \code{results} list for accumulation
#' @param source Starting node of the path
#' @param current Current node in traversal
#' @param path_nodes Nodes visited so far (excluding source)
#' @param path_signs Signs of edges traversed so far
#' @param current_prob Accumulated transition probability
#' @param remaining_hops Number of hops remaining
#' @param n Number of nodes
#' @param P_pos Positive transition matrix
#' @param P_neg Negative transition matrix
#' @param min_contribution Minimum contribution threshold
#' @keywords internal
.enumeratePaths <- function(acc, source, current, path_nodes, path_signs,
                            current_prob, remaining_hops, n,
                            P_pos, P_neg, min_contribution) {
  if (remaining_hops == 0) return(invisible(NULL))

  for (j in seq_len(n)) {
    for (sign_type in c("+", "-")) {
      prob <- if (sign_type == "+") P_pos[j, current] else P_neg[j, current]
      if (prob == 0) next

      new_prob <- current_prob * prob
      if (min_contribution > 0 && new_prob < min_contribution) next

      new_path <- c(path_nodes, j)
      new_signs <- c(path_signs, sign_type)

      if (remaining_hops == 1L) {
        net <- prod(ifelse(new_signs == "+", 1L, -1L))
        acc$results[[length(acc$results) + 1L]] <- list(
          source = source, target = j,
          hop = length(new_path),
          path = paste(c(source, new_path), collapse = "->"),
          signs = paste(new_signs, collapse = ","),
          net_sign = as.integer(net),
          contribution = new_prob
        )
      } else {
        .enumeratePaths(acc, source, j, new_path, new_signs, new_prob,
                        remaining_hops - 1L, n, P_pos, P_neg,
                        min_contribution)
      }
    }
  }
  invisible(NULL)
}
