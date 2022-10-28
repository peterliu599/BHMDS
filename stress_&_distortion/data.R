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

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# karate
data(karate)
g <- graph_from_adjacency_matrix(karate$adjacency)
D <- distances(g)
kappa = 1

# phylogenetic tree
# g <- read.graph("phylo_tree.edges")
# D <- distances(g)
# kappa = 0.14

# cs phd
# g <- read.graph("ca-CSphd.edges")
# D <- distances(g)
# kappa = 0.55

# wordnet mammal subtree
# D <- as.matrix(read.csv("dist_wordnet_mammal.csv")[,-1])
# kappa = 2.06

## grid search for curvature
# kappa_vec <- seq(from = 0.1, to = 15, by = 0.01)
# l = length(kappa_vec)
# d = 2
# stress <- rep(0, l)
# for (i in 1:l) {
#   kappa = kappa_vec[i]
#   stress[i] = hydraPlus(D, dim = d, curvature = kappa, maxit = 5000)$stress
# }
# kappa_vec[which.min(stress)]


# Number of data 
n = dim(D)[1]
# Dimension 
d = 2

mu0 <- c(1, rep(0, d))
m = choose(n, 2)

## bhmds
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
# maxiter = 20000
# delta0 <- compute_delta(VtoX(V0), kappa)
# sigg0 = compute_SSR(D, delta0) / m
# a = 5
# b = (a - 1) * sigg0
# alpha = 0.5
# beta0 = apply(V0, 2, var) / 2
# constant = 1
# burnin = 3000
# 
# result_bhmds <- main_bhmds(n, d, D, V0, X0, sigg0, a, alpha, maxiter, constant, beta0, burnin, kappa)

## bmds
# coef <- sqrt(sum(D^2)) / sqrt(dim(D)[1] ^ 2)
# smacof <- mds(D, ndim = d)
# X0 <- smacof$conf * coef
# m = choose(n, 2)
# delta0 <- compute_delta_E(X0)
# sigg0 = compute_SSR(D, delta0) / m
# a = 5
# b = (a - 1) * sigg0
# alpha = 0.5
# beta0 = apply(X0, 2, var) / 2
# constant = 1
# 
# maxiter = 20000
# burnin = 3000
# result_bmds <- main_bmds(n, d, D, V0, sigg0, a, alpha, maxiter, constant, beta0, burnin)

file_name_bhmds = "bhmds_karate.csv"
file_name_bmds = "bmds_karate.csv"

# file_name_bhmds = "bhmds_tree.csv"
# file_name_bmds = "bmds_tree.csv"
# 
# file_name_bhmds = "bhmds_csphd.csv"
# file_name_bmds = "bmds_csphd.csv"
# 
# file_name_bhmds = "bhmds_wordnet.csv"
# file_name_bmds = "bmds_wordnet.csv"

result_bhmds <- as.matrix(read.csv(file_name_bhmds)[,-1])
result_bmds <- as.matrix(read.csv(file_name_bmds)[,-1])

# write.csv(result_bhmds, file_name_bhmds)
# write.csv(result_bmds, file_name_bmds)

## Stress evaluation
# bhmds stress
sqrt(sum((result_bhmds - D)^2)) / sqrt(sum(D^2))

# bmds stress
sqrt(sum((result_bmds - D)^2)) / sqrt(sum(D^2))

# hydraPlus stress
hydra$stress / sqrt(sum(D^2))

# hydra stress
hydra(D, dim = d, curvature = kappa)$stress / sqrt(sum(D^2))

# smacof stress
smacof$stress

# cmds stress
Xtemp = cmdscale(D, k =2)
sqrt(sum((compute_delta_E(Xtemp) - D)^2)) / sqrt(sum(D^2))

## Distortion evaluation
# bhmds distortion
sum(abs(result_bhmds - D) / (D + diag(n))) / m

# bmds distortion
sum(abs(result_bmds - D) / (D + diag(n))) / m

# hydraPlus distortion
sum(abs(compute_delta(VtoX(V0), 1) - D) / (D + diag(n))) / m

# hydra distortion
hydra <- hydra(D, dim = d, curvature = kappa, control = list(return.lorentz = T))
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
sum(abs(compute_delta(VtoX(V0), 1) - D) / (D + diag(n))) / m

# smacof distortion
coef <- sqrt(sum(D^2)) / sqrt(dim(D)[1] ^ 2)
smacof <- mds(D, ndim = d)
X0 <- smacof$conf * coef
m = choose(n, 2)
delta0 <- compute_delta_E(X0)
sum(abs(delta0 - D) / (D + diag(n))) / m

# cmds distortion
Xtemp = cmdscale(D, k =2)
sum(abs(compute_delta_E(Xtemp)- D) / (D + diag(n))) / m



