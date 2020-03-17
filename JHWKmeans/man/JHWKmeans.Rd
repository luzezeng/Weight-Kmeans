\name{JHWKmeans}
\alias{JHWKmeans}
\title{Hello, World!}
\usage{
library(JHWKmeans)
result = kmeans.weighting(x, k, beta = 8, iter.max = 10)
U = result$U
Z = result$Z
W = result$W
}
\description{
X is the input data as a n * m matrix, n is the sample total count and m is the variable count.
U is an n * k partition matrix, U[i, l] indicates that object i is allocated to cluster l;
Z = {Z[1], Z[2], ..., Z[k]} is a set of k vectors representing the centroids of the k clusters;
W = {W[1], W[2], ..., W[m]} is the weight for m variables.
}
\examples{
library(JHWKmeans)
set.seed(11)
x <- matrix(rnorm(50 * 70), ncol = 70)
x[1 : 25, 1 : 20] <- x[1 : 25, 1 : 20] + 1
x <- scale(x, TRUE, TRUE)
result = kmeans.weighting(x, k, beta = 8, iter.max = 10)
print(result)
}
