[![GitHub Actions](https://github.com/rikenbit/symTensor/actions/workflows/build_test_push.yml/badge.svg)](https://github.com/rikenbit/symTensor/actions/workflows/build_test_push.yml)

# symTensor
R package for Symmetric Matrix and Tensor Decomposition

## Installation

~~~~
git clone https://github.com/rikenbit/symTensor/
R CMD INSTALL symTensor
~~~~
or type the code below in the R console window
~~~~
library(devtools)
devtools::install_github("rikenbit/symTensor")
~~~~

## Functions

- **isSymmetricTensor**: Check whether a tensor is symmetric
- **Symmetrize**: Symmetrize a tensor
- **PageRank**: PageRank algorithm
- **LabelPropagation**: Label Propagation algorithm
- **symNMF**: Symmetric Non-negative Matrix Factorization
- **HigherOrderPower**: Higher-Order Power Method for tensor decomposition
- **TOPHITS**: Tensor-based algorithm for hub and authority analysis

## References

- **Symmetric NMF**: Da Kuang, Chris Ding, and Haesun Park, "Symmetric Nonnegative Matrix Factorization for Graph Clustering", SDM, 2012
- **Beta-divergence**: Cedric Fevotte and Jerome Idier, "Algorithms for Nonnegative Matrix Factorization with the Beta-Divergence", Neural Computation, 2011
- **Beta-divergence**: Mikkel N. Schmidt, Masayuki Nakano, et al., "Nonnegative Matrix Factorization with Beta-divergence", ISMIR, 2010
- **PageRank**: Sergey Brin and Lawrence Page, "The Anatomy of a Large-Scale Hypertextual Web Search Engine", Computer Networks and ISDN Systems, 1998
- **Higher-Order Power Method**: Tamara G. Kolda and Jackson R. Mayo, "Shifted Power Method for Computing Tensor Eigenpairs", Journal of Computational and Applied Mathematics, 2011
- **Higher-Order Power Method**: Lieven De Lathauwer, Bart De Moor, and Joos Vandewalle, "On the Best Rank-1 and Rank-(R1,R2,...,RN) Approximation of Higher-Order Tensors", SIAM Journal on Matrix Analysis and Applications, 2000
- **TOPHITS**: Tamara G. Kolda and Brett W. Bader, "The TOPHITS Model for Higher-Order Web Link Analysis", Workshop on Link Analysis, Counterterrorism and Security, 2006
- **Label Propagation**: Xiaojin Zhu and Zoubin Ghahramani, "Learning from Labeled and Unlabeled Data with Label Propagation", CMU-CALD-02-107, 2002

## Contributing
If you have suggestions for how `symTensor` could be improved, or want to report a bug, open an issue! We'd love all and any contributions.

For more details, check the [CONTRIBUTING.md](CONTRIBUTING.md).

## Authors
- Koki Tsuyuzaki
