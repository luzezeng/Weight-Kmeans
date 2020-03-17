library(JHWKmeans)
set.seed(11)
x <- matrix(rnorm(50 * 70), ncol = 70)
x[1 : 25, 1 : 20] <- x[1 : 25, 1 : 20] + 1
x <- scale(x, TRUE, TRUE)
result = kmeans.weighting(x, 2)
print(result)
