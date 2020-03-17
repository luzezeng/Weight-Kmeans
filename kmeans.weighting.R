self.dist <- function(x, y) {
  
  x = as.matrix(as.numeric(x))
  y = as.matrix(as.numeric(y))
  
  if(nrow(x) != nrow(y))
    stop("can't calculate distance between 2 different vector")
  
  if(ncol(x) > 1 | ncol(y) > 1)
    stop("column number should be one !")
  
  return((x - y)*(x - y))
  
}

check.UZW <- function(U, Z, W) {
  if(missing(U) | missing(Z) | missing(W))
    stop("U or Z or W is missing !")
  
  k = ncol(U)
  n = nrow(U)
  m = ncol(Z)
  
  if(nrow(Z) != k)
    stop("row of Z need to be equal with k !")
  
  if(nrow(W) != 1 | ncol(W) != m)
    stop("row of W need to be one, and col of W need to be equal with m !")
  
  if(sum(W) != 1)
    stop("sum of W should be one !")
  
  for (i in 1 : n) {
    if(sum(U[i, ]) != 1)
      stop("sum of every row for U should be one !")
    if(is.null(which(U[i, ] != 0 & U[i, ] != 1)))
      stop("item in U should be 0 or 1 !")
  }
  
  return(TRUE)
}

calculate.P <- function(X, U, Z, W, beta = 8) {
  
  if(!check.UZW(U, Z, W))
    stop("NOT pass the check on U, Z, W !")
  
  k = ncol(U)
  n = nrow(U)
  m = ncol(Z)
  
  P = 0;
  for (l in 1 : k) {
    for (i in 1 : n) {
      for (j in 1 : m) {
        P = P + U[i, l] * (W[j] ^ beta) * (X[i, j] * Z[l, j])
      }
    }
  }
  
  return(P)
  
}

decide.U <- function(X, U, Z, W, beta = 8) {
  
  if(!check.UZW(U, Z, W))
    stop("NOT pass the check on U, Z, W !")
  
  k = ncol(U)
  n = nrow(U)
  m = ncol(Z)
  
  distance = matrix(data = 0, nrow = n, ncol = k)
  
  for (i in 1 : n) {
    for (l in 1 : k) {
      temp.dist = 0
      for (j in 1 : m) {
        temp.dist = temp.dist + W[j] * ((X[i, j] - Z[l, j])^2)
      }
      distance[i, l] = temp.dist
    }
  }
  
  for (i in 1 : n) {
    min.index = which.min(distance[i, ])
    U[i, min.index] = 1
    for (l in 1 : k) {
      if(l != min.index)
        U[i, l] = 0
    }
  }
  
  return(U)
  
}

decide.Z <- function(X, U, Z, W, beta = 8) {
  
  if(!check.UZW(U, Z, W))
    stop("NOT pass the check on U, Z, W !")
  
  k = ncol(U)
  n = nrow(U)
  m = ncol(Z)
  
  for (l in 1 : k) {
    for (j in 1 : m) {
      temp.nume = 0
      temp.deno = 0
      for (i in 1 : n) {
        temp.nume = temp.nume + U[i, l] * X[i, j]
        temp.deno = temp.deno + U[i, l]
      }
      Z[l, j] = temp.nume / temp.deno
    }
  }
  
  return(Z)
  
}

decide.W <- function(X, U, Z, W, beta = 8) {
  
  if(!check.UZW(U, Z, W))
    stop("NOT pass the check on U, Z, W !")
  
  k = ncol(U)
  n = nrow(U)
  m = ncol(Z)
  
  D = matrix(0, nrow = 1, ncol = m)
  
  for (j in 1 : m) {
    temp.d = 0;
    for (l in 1 : k) {
      for (i in 1 : n) {
        temp.d = temp.d + U[i, l]*((X[i, j] - Z[l, j]) ^ 2)
      }
    }
    D[j] = temp.d;
  }
  
  h = 0
  for (j in 1 : m) {
    if(D[j] != 0)
      h = h + 1
  }
  
  for (j in 1 : m) {
    if(D[j] == 0) {
      W[j] = 0
    } else {
      temp.w = 0
      for (t in 1 : h) {
        temp.w = temp.w + ((D[j] / D[t]) ^ (1 / (beta - 1)))
      }
      W[j] = 1 / temp.w
    }
  }
  
  return(W)
  
}

kmeans.weighting = function(X, k, beta = 8, iter.max = 10) {
  
  if(missing(X))
    stop("the input data frame X must be provided !")
  
  if(missing(k))
    stop("the excepted cluster count k must be provided !")
  
  # reform the type as integer
  n = as.integer(nrow(X)) # n samples
  m = as.integer(ncol(X)) # m features for each sample
  k = as.integer(k) # expect to converge to k clusters
  beta = as.integer(beta) # parameter beta for W
  iter.max = as.integer(iter.max) # maximum iterator times
  
  if(k > n)
    stop("the row of data frame X must be larger than k !")
  
  if(iter.max <= 0)
    stop("the maximum iterator count must be larger than 0 !")
  
  # initialize Z as center point for K clusters
  Z = X[sample(x = n, size = k, replace = F), ]
  
  # initialize W as weight for M features
  random.serial = runif(m, min = 0, max = 1)
  W = t(as.matrix(random.serial / sum(random.serial)))
  
  # initialize U as one partition
  U = matrix(data = 0, nrow = n, ncol = k)
  U[, 1] = 1
  
  P = calculate.P(X, U, Z, W, beta)
  
  # begin to iterate
  for (step in 1 : iter.max) {
    
    # 1. update for U
    U = decide.U(X, U, Z, W, beta)
    Pt = calculate.P(X, U, Z, W, beta)
    if(Pt == P) {
      result <- list(U, Z, W)
      names(result) <- c("U", "Z", "W")
      return(result)
    }
    
    
    # 2. update for Z
    P = Pt
    Z = decide.Z(X, U, Z, W, beta)
    Pt = calculate.P(X, U, Z, W, beta)
    if(Pt == P) {
      result <- list(U, Z, W)
      names(result) <- c("U", "Z", "W")
      return(result)
    }
    
    # 3. update for W
    P = Pt
    W = decide.W(X, U, Z, W, beta)
    Pt = calculate.P(X, U, Z, W, beta)
    if(Pt == P) {
      result <- list(U, Z, W)
      names(result) <- c("U", "Z", "W")
      return(result)
    }
    P = Pt
    
  }
  
}

set.seed(11)
x <- matrix(rnorm(50 * 70), ncol = 70)
x[1 : 25, 1 : 20] <- x[1 : 25, 1 : 20] + 1
x <- scale(x, TRUE, TRUE)
result = kmeans.weighting(x, 2)