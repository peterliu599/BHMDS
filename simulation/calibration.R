library(igraph)
library(mvtnorm)
library(LaplacesDemon)
library(truncnorm)
library(hydra)
library(smacof)

# Lorentzian product
Lorentz <- function (x, y) {
  res <- -x[1] * y[1] + sum(x[-1] * y[-1])
  return(res)
}

# Convert V to X with respect to mu0
VtoX <- function(V) {
  X <- matrix(rep(0, n * (d + 1)), ncol = (d + 1))
  mu0 <- c(1, rep(0, d))
  for (i in 1:n) {
    vtemp <- c(0, V[i,])
    lor <- sqrt(Lorentz(vtemp, vtemp))
    X[i,] <- cosh(lor) * mu0 + sinh(lor) / lor * vtemp
  }
  return(X)
}

# Compute delta
compute_delta <- function(X, kappa) {
  temp = X
  temp[, 1] <- -temp[, 1]
  temp = -temp %*% t(X); diag(temp) = 1
  temp = acosh(temp) * sqrt(1 / kappa)
  return(temp)
}

# Compute SSR
compute_SSR <- function(D, delta) {
  tempSSR = sum((D[upper.tri(D)] - delta[upper.tri(delta)]) ^ 2)
  return (tempSSR)
}

# Compute SSR via matrix
compute_SSR_mat <- function(D, V, kappa) {
  kappa = kappa
  Xtemp <- matrix(rep(0, n * (d + 1)), ncol = (d + 1))
  Xtemp = VtoX(V)
  deltatemp = compute_delta(Xtemp, kappa)
  res = compute_SSR(D, deltatemp)
  return(res)
}

# Compute delta for Euclidean
compute_delta_E <- function (X) {
  return (as.matrix(dist(X), ncol = n))
}

# Compute SSR via matrix for Euclidean
compute_SSR_mat_E <- function(D, X) {
  deltatemp = compute_delta_E(X)
  res = compute_SSR(D, deltatemp)
  return(res)
}

# Centering and rotating for bmds
crotX <- function(X) {
  xbar = colMeans(X)
  for (i in 1:n) {
    X[i,] = X[i,] - xbar
  }
  Xcov = t(X) %*% X / n
  Xout = X %*% eigen(Xcov)$vector
  return(Xout)
}

# Single update step for bmds
update_xvec <- function(D, X, id, sigma2, constant, lambdas) {
  Xold = Xtgt = X
  stepsize = sqrt(sigma2 * constant / (n - 1))
  for (i in 1:d) {
    Xtgt[id, i] = Xtgt[id, i] + rnorm(1, 0, stepsize)
  }
  sigma = sqrt(sigma2)
  xold = Xold[id,]
  xtgt = Xtgt[id,]
  deltgt = delold = rep(0, n)
  for (i in 1:n) {
    deltgt[i] = sqrt(sum((xtgt - Xtgt[i,])^2))
    delold[i] = sqrt(sum((xold - Xold[i,])^2))
  }
  Q1tgt = sum((D[id,] - deltgt)^2) / sigma2
  Q1old = sum((D[id,] - delold)^2) / sigma2
  Q2tgt = sum(xtgt^2 * (1/lambdas))
  Q2old = sum(xold^2 * (1/lambdas))
  t3tgt = sum(log(pnorm(deltgt / sigma)))
  t3old = sum(log(pnorm(delold / sigma)))
  ftgt = -(Q1tgt+Q2tgt)/2.0 - t3tgt
  fold = -(Q1old+Q2old)/2.0 - t3old 
  fratio = min(exp(ftgt - fold), 1)
  if (runif(1) <= fratio) {
    return(xtgt)
  }
  return(xold)
}

# bmds
main_bmds <- function(n, d, D, V0, sigg0, a, alpha, maxiter, constant, betas, burnin) {
  tracedelta = matrix(rep(0, n^2), ncol = n)
  
  m = choose(n, 2)
  Xold = crotX(X0)
  Xnew = matrix(rep(0, n * d), ncol = d)
  
  SSRnew = 0
  SSRold = compute_SSR_mat_E(D, Xold)
  SSRmin = SSRold
  
  sigma2 = sigg0
  sigtmp = 0
  vecs = rep(0, d)
  lambdas = rep(0, d)
  tmprow = rep(0, d)
  
  varalpha = 0.0
  varbeta  = 0.0
  varvar = 0.0
  varratio = 0.0
  
  for (i in 1:maxiter) {
    print(i)
    for (j in 1:d) {
      vecs[j] = var(Xold[,j]) * n
    }
    for (j in 1:d) {
      lambdas[j] = rinvgamma(1, alpha + n/2, betas[j] + vecs[j] / 2)
    }
    Xnew = Xold
    for (j in 1:n) {
      tmprow = update_xvec(D, Xnew, j, sigma2, constant, lambdas)
      Xnew[j,] = tmprow
    }
    SSRnew = compute_SSR_mat_E(D, Xnew)
    Xnew = crotX(Xnew)
    
    varalpha = m/2 + a
    varbeta  = SSRnew/2 + b
    varvar   = (varbeta*varbeta)/((varalpha-1)*(varalpha-1)*(varalpha-2))
    sigtmp = sigma2 + rnorm(1, 0, sqrt(constant*varvar))
    if (sigtmp > 0){ 
      varratio = dinvgamma(sigtmp,varalpha,varbeta)/dinvgamma(sigma2,varalpha,varbeta)
      if (varratio > 1){
        varratio = 1.0
      }
      if (runif(1) <= varratio){
        sigma2 = sigtmp
      }
    }
    if (i > burnin) {
      if (SSRmin > SSRnew) {
        tracedelta <- compute_delta_E(Xnew)
        SSRmin = SSRnew
      }
    }
    SSRold = SSRnew
    Xold   = Xnew
    if (i == burnin) {
      SSRmin = compute_SSR_mat_E(D, Xold)
    }
  }
  return (tracedelta)
}

# Single update step for bhmds
update_vvec <- function (D, V, X, id, sigma2, constant, lambdas) {
  Vnew = V
  stepsize = sqrt(sigma2 * constant / (n - 1))
  # proposal 
  for (i in 1:d) {
    Vnew[id, i] = Vnew[id, i] + rnorm(1, 0, stepsize) 
  }
  sigma = sqrt(sigma2)
  
  Xnew = X
  lornew = sqrt(Lorentz(c(0, Vnew[id,]), c(0, Vnew[id,])))
  Xnew[id,] = cosh(lornew) * mu0 + sinh(lornew) / lornew * c(0, Vnew[id,])
  xold = X[id,]
  xnew = Xnew[id,]
  xold[1] = -xold[1]; xnew[1] = -xnew[1]
  delnew = -xnew %*% t(Xnew); delold = -xold %*% t(X)
  delnew[id] = delold[id] = 1 # Just in case
  delnew = acosh(delnew) * (1 / sqrt(kappa)); delold = acosh(delold) * (1 / sqrt(kappa))
  
  Q1new = sum((D[id,] - delnew)^2) / sigma2
  Q1old = sum((D[id,] - delold)^2) / sigma2
  Q2new = sum(Vnew[id,]^2 * (1/lambdas))
  Q2old = sum(V[id,]^2 * (1/lambdas))
  t3new = sum(log(pnorm(delnew / sigma)))
  t3old = sum(log(pnorm(delold / sigma)))
  fnew = -(Q1new+Q2new)/2.0 - t3new
  fold = -(Q1old+Q2old)/2.0 - t3old
  fratio = min(exp(fnew - fold), 1)
  if (runif(1) <= fratio) {
    return(c(Vnew[id,], Xnew[id,]))
  } else {
    return(c(V[id,], X[id,]))
  }
}

# bhmds
main_bhmds <- function(n, d, D, V0, X0, sigg0, a, alpha, maxiter, constant, betas, burnin, kappa) {
  kappa = kappa
  tracedelta = matrix(rep(0, n^2), ncol = n)
  m = choose(n, 2)
  V = V0
  X = X0
  Vnew = matrix(rep(0, n * d), ncol = d)
  
  SSRnew = 0
  SSRold = compute_SSR_mat(D, V, kappa)
  
  sigma2 = sigg0
  sigtmp = 0
  vecs = rep(0, d)
  lambdas = rep(0, d)
  tmprow = rep(0, d)
  
  varalpha = 0
  varbeta  = 0
  varvar = 0
  varratio = 0
  flag = 0
  
  for (i in 1:maxiter) {
    for (j in 1:d) {
      vecs[j] = var(V[,j]) * n
    }
    for (j in 1:d) {
      lambdas[j] = rinvgamma(1, alpha + n/2, betas[j] + vecs[j] / 2)
    }
    for (j in 1:n) {
      tmprow = update_vvec(D, V, X, j, sigma2, constant, lambdas)
      V[j,] = tmprow[1:d]
      X[j,] = tmprow[-1:-d]
    }
    SSRnew = compute_SSR_mat(D, V, kappa)
    
    varalpha = m/2 + a
    varbeta  = SSRnew/2 + b
    varvar   = (varbeta*varbeta)/((varalpha-1)*(varalpha-1)*(varalpha-2))
    sigtmp   = sigma2 + rnorm(1, 0, sqrt(constant*varvar))
    if (sigtmp > 0){ 
      varratio = dinvgamma(sigtmp,varalpha,varbeta)/dinvgamma(sigma2,varalpha,varbeta)
      if (varratio > 1){
        varratio = 1.0
      }
      if (runif(1) <= varratio){
        sigma2 = sigtmp
      }
    }
    if (i > burnin) {
      if (SSRmin > SSRnew) {
        tracedelta = compute_delta(X, kappa)
        SSRmin = SSRnew
      } 
    }
    SSRold = SSRnew
    if (i == burnin) {
      SSRmin = compute_SSR_mat(D, V, kappa)
      tracedelta = compute_delta(X, kappa)
    }
    print(i)
  }
  return (tracedelta)
}

# Number of data 
n = 50
## Paramter variation
d_vec = c(2, 3, 4, 5)
# Dimension
d = sample(d_vec, 1)
# mean of wrapped normal distribution
mu_vec <- rnorm(d, 0, 2)
# variance of wrapped normal distribution
sigma_vec <- runif(d, min = 5, max = 10)
# hyperbolic curvature
kappa <- runif(1, min = 0.2, max = 2)

mu0 <- c(1, rep(0, d))
m = choose(n, 2)
V <- matrix(rep(0, n * d), ncol = d)

record_true <- rep(0, 500)
record_bhmds <- record_bmds <- rep(0, 500)

for (simindex in 1:20) {
  for (i in 1:n) {
    V[i,] <- rmvnorm(1, mu_vec, diag(sigma_vec))
  }
  
  delta <- compute_delta(VtoX(V), kappa)
  
  tempvec <- rep(0, 25)
  for (i in 1:25) {
    tempvec[i] = delta[(2 * i - 1), (2 * i)]
  }
  record_true[((simindex - 1) * 25 + 1):(simindex * 25)] = tempvec
  # True sigma
  sigma <- 1
  
  # Observed D matrix
  D <- matrix(rep(0, n^2), ncol = n)
  for (i in 1:n) {
    for (j in 1:n) { 
      if (i < j) {
        D[i, j] = rtruncnorm(1, a = 0, b = Inf, delta[i, j], sigma)
      }
    }
  }
  D = D + t(D)
  
  mu0 <- c(1, rep(0, d))
  m = choose(n, 2)
  
  ## bhmds
  # Initialization
  hydra <- hydraPlus(D, dim = d, curvature = kappa, control = list(return.lorentz = T), maxit = 5000)
  nancount = 0
  tempx0 <- (1 + hydra$r^2) / (1 - hydra$r^2)
  tempx <- hydra$directional * sqrt(tempx0^2 - 1)
  tempX <- cbind(tempx0, tempx)
  Vini <- matrix(rep(0, n*(d+1)), ncol = (d + 1))
  nancount = 0
  for (i in 1:n) {
    alpha = -Lorentz(mu0, tempX[i,])
    if (is.nan(acosh(alpha)) == TRUE) {
      alpha = 1.01
      nancount = nancount + 1
    }
    Vini[i,] = acosh(alpha) / sqrt(alpha^2 - 1) * (tempX[i,] - alpha * mu0)
  }
  V0 <- Vini[, -1]
  X0 <- VtoX(V0)
  maxiter = 20000
  delta0 <- compute_delta(VtoX(V0), kappa)
  sigg0 = compute_SSR(D, delta0) / m
  a = 5
  b = (a - 1) * sigg0
  alpha = 0.5
  beta0 = apply(V0, 2, var) / 2
  constant = 1
  burnin = 3000
  
  result_bhmds <- main_bhmds(n, d, D, V0, X0, sigg0, a, alpha, maxiter, constant, beta0, burnin, kappa)
  tempvec <- rep(0, 25)
  for (i in 1:25) {
    tempvec[i] = result_bhmds[(2 * i - 1), (2 * i)]
  }
  record_bhmds[((simindex - 1) * 25 + 1):(simindex * 25)] = tempvec
  
  ## bmds
  # Initialization
  coef <- sqrt(sum(D^2)) / sqrt(dim(D)[1] ^ 2)
  smacof <- mds(D, ndim = d)
  X0 <- smacof$conf * coef
  m = choose(n, 2)
  delta0 <- compute_delta_E(X0)
  sigg0 = compute_SSR(D, delta0) / m
  a = 5
  b = (a - 1) * sigg0
  alpha = 0.5
  beta0 = apply(X0, 2, var) / 2
  constant = 1
  
  maxiter = 20000
  burnin = 3000
  result_bmds <- main_bmds(n, d, D, V0, sigg0, a, alpha, maxiter, constant, beta0, burnin)
  tempvec <- rep(0, 25)
  for (i in 1:25) {
    tempvec[i] = result_bmds[(2 * i - 1), (2 * i)]
  }
  record_bmds[((simindex - 1) * 25 + 1):(simindex * 25)] = tempvec
}

file_name_true = paste("true_iter_num_", slurm_arrayid, ".csv", sep = "")
file_name_bhmds = paste("bhmds_iter_num_", slurm_arrayid, ".csv", sep = "")
file_name_bmds = paste("bmds_iter_num_", slurm_arrayid, ".csv", sep = "")
write.csv(record_true, file_name_true)
write.csv(record_bhmds, file_name_bhmds)
write.csv(record_bmds, file_name_bmds)

# true_res <- matrix(rep(0, 500 * 1000), ncol = 500)
# for (i in 1:1000) {
#   file_name <- paste("true_iter_num_", i, ".csv", sep = "")
#   file <- read.csv(file_name)[,-1]
#   true_res[i, ] = file
# }
# write.csv(true_res, "comb_true.csv")
# 
# bhmds_res <- matrix(rep(0, 500 * 1000), ncol = 500)
# for (i in 1:1000) {
#   file_name <- paste("bhmds_iter_num_", i, ".csv", sep = "")
#   file <- read.csv(file_name)[,-1]
#   bhmds_res[i, ] = file
# }
# write.csv(bhmds_res, "comb_bhmds.csv")
# 
# bmds_res <- matrix(rep(0, 500 * 1000), ncol = 500)
# for (i in 1:1000) {
#   file_name <- paste("bmds_iter_num_", i, ".csv", sep = "")
#   file <- read.csv(file_name)[,-1]
#   bmds_res[i, ] = file
# }
# write.csv(bmds_res, "comb_bmds.csv")

# bhmds <- as.matrix(read.csv("comb_bhmds.csv")[,-1])
# bmds <- as.matrix(read.csv("comb_bmds.csv")[,-1])
# true <- as.matrix(read.csv("comb_true.csv")[,-1])
# 
# random_draw <- as.matrix(read.csv("random_draw.csv")[,-1])
# max(true); min(true)
# seq <- seq(from = 0, to = 20, length.out = 2000)
# record_h <- record_e <- rep(0, 2000)
# for (i in 1:2000) {
#   print(i)
#   x = seq[i]
#   tempsum_h <- tempsum_e <- 0
#   for (j in 1:1000) {
#     ecdf_h <- ecdf(bhmds[j,])
#     ecdf_e <- ecdf(bmds[j,])
#     tempsum_h = tempsum_h + ecdf_h(x)
#     tempsum_e = tempsum_e + ecdf_e(x)
#   }
#   record_h[i] <- tempsum_h / 1000 - mean(random_draw <= x)
#   record_e[i] <- tempsum_e / 1000 - mean(random_draw <= x)
# }

# write.csv(record_h, "record_h.csv")
# write.csv(record_e, "record_e.csv")

# record_e <- as.matrix(read.csv("record_e.csv")[,-1])
# record_h <- as.matrix(read.csv("record_h.csv")[,-1])
# 
# dat <- data.frame(y1 = record_h, y2 = record_e, x = seq)
# library(tidyverse)
# ggplot(data = dat) + geom_line(aes(x = x, y = y1), color = "blue") + geom_line(aes(x = x, y = y2)) + geom_hline(yintercept = 0, color = "red") + xlab("Threshold Value")+ 
#   ylab("F.fcast - F.obs") + xlim(0, 17.5) + 
#   theme(text = element_text(size=20), axis.text=element_text(size=20))+ 
#   ggtitle("Marginal calibration") + theme(plot.title = element_text(hjust = 0.5))
