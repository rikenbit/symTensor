library(symTensor)
library(rTensor)

context("symTensor")

test_that("isSymmetricTensor works for matrices", {
  A <- matrix(c(1, 2, 2, 3), 2, 2)
  expect_identical(isSymmetricTensor(A), TRUE)

  B <- matrix(c(1, 2, 3, 4), 2, 2)
  expect_identical(isSymmetricTensor(B), FALSE)
})

test_that("isSymmetricTensor works for 3-mode tensors", {
  arr <- array(0, dim = c(3, 3, 3))
  for (i in 1:3) for (j in 1:3) for (k in 1:3) {
    arr[i, j, k] <- i + j + k
  }
  expect_identical(isSymmetricTensor(as.tensor(arr)), TRUE)
})

test_that("Symmetrize works for matrices", {
  A <- matrix(c(1, 2, 3, 4), 2, 2)
  S <- Symmetrize(A)
  expect_identical(isSymmetricTensor(S), TRUE)
})

test_that("Symmetrize works for tensors", {
  arr <- array(rnorm(27), dim = c(3, 3, 3))
  S <- Symmetrize(as.tensor(arr))
  expect_identical(isSymmetricTensor(S), TRUE)
})

test_that("PageRank basic functionality", {
  A <- matrix(c(0, 1, 1, 0, 1, 0, 0, 1, 0, 1, 0, 1, 1, 0, 1, 0), 4, 4)
  pr <- PageRank(A)
  expect_lt(abs(sum(pr$vector) - 1), 1e-10)
  expect_identical(pr$converged, TRUE)
})

test_that("PageRank with personalization", {
  A <- matrix(c(0, 1, 1, 0, 1, 0, 0, 1, 0, 1, 0, 1, 1, 0, 1, 0), 4, 4)
  v <- c(1, 0, 0, 0)
  pr <- PageRank(A, personalization = v)
  expect_identical(pr$converged, TRUE)
  expect_identical(which.max(pr$vector), 1L)
})

test_that("LabelPropagation basic", {
  A <- matrix(0, 6, 6)
  A[1, 2] <- A[2, 1] <- 1
  A[1, 3] <- A[3, 1] <- 1
  A[2, 3] <- A[3, 2] <- 1
  A[4, 5] <- A[5, 4] <- 1
  A[4, 6] <- A[6, 4] <- 1
  A[5, 6] <- A[6, 5] <- 1
  A[3, 4] <- A[4, 3] <- 0.1
  seeds <- c("1" = 1, "6" = 2)
  result <- LabelPropagation(A, seeds)
  expect_identical(result$labels[1:3], c(1L, 1L, 1L))
  expect_identical(result$labels[4:6], c(2L, 2L, 2L))
})

test_that("symNMF MM Frobenius (Beta=2) error decreases", {
  set.seed(123)
  M_true <- matrix(0, 10, 2)
  M_true[1:5, 1] <- runif(5, 0.5, 1)
  M_true[6:10, 2] <- runif(5, 0.5, 1)
  Omega <- M_true %*% t(M_true)
  Omega <- Omega + matrix(runif(100, 0, 0.01), 10, 10)
  Omega <- (Omega + t(Omega)) / 2

  result <- symNMF(Omega, K = 2, model = "MM", Beta = 2, num.iter = 200)
  expect_gt(length(result$RecError), 1)
  expect_lte(tail(result$RecError, 1), result$RecError[1])
})

test_that("symNMF MSM Frobenius (Beta=2) error decreases", {
  set.seed(42)
  M_true <- matrix(runif(20, 0, 1), 10, 2)
  S_true <- matrix(c(2, 0.5, 0.5, 3), 2, 2)
  Omega <- M_true %*% S_true %*% t(M_true)

  result <- symNMF(Omega, K = 2, model = "MSM", Beta = 2, num.iter = 200)
  expect_gt(length(result$RecError), 1)
  expect_lte(tail(result$RecError, 1), result$RecError[1])
})

test_that("symNMF KL divergence (Beta=1) error decreases", {
  set.seed(42)
  M_true <- matrix(runif(20, 0.1, 1), 10, 2)
  Omega <- M_true %*% t(M_true)
  Omega <- Omega + matrix(runif(100, 0, 0.01), 10, 10)
  Omega <- (Omega + t(Omega)) / 2

  result <- symNMF(Omega, K = 2, model = "MM", Beta = 1, num.iter = 200)
  expect_gt(length(result$RecError), 1)
  expect_lte(tail(result$RecError, 1), result$RecError[1])
})

test_that("symNMF IS divergence (Beta=0) error decreases", {
  set.seed(42)
  M_true <- matrix(runif(20, 0.1, 1), 10, 2)
  Omega <- M_true %*% t(M_true)
  Omega <- Omega + matrix(runif(100, 0.01, 0.1), 10, 10)
  Omega <- (Omega + t(Omega)) / 2

  result <- symNMF(Omega, K = 2, model = "MM", Beta = 0, num.iter = 200)
  expect_gt(length(result$RecError), 1)
  expect_lte(tail(result$RecError, 1), result$RecError[1])
})

test_that("HigherOrderPower recovers rank-1 tensor", {
  set.seed(42)
  n <- 5
  v1 <- rnorm(n)
  v1 <- v1 / sqrt(sum(v1^2))
  lam <- 3.0
  arr <- lam * outer(outer(v1, v1), v1)
  tnsr <- as.tensor(arr)

  result <- HigherOrderPower(tnsr, R = 1)
  expect_lt(abs(abs(result$lambdas[1]) - 3.0), 0.1)
  cosine <- abs(sum(result$vectors[, 1] * v1))
  expect_gt(cosine, 0.95)
})

test_that("TOPHITS runs without error", {
  set.seed(42)
  arr <- array(abs(rnorm(125)), dim = c(5, 5, 5))
  tnsr <- as.tensor(arr)
  result <- TOPHITS(tnsr, R = 2)
  expect_identical(length(result$rankings), 3L)
  expect_identical(length(result$rankings[[1]]), 5L)
})
