library(igraph)
library(mvtnorm)
library(LaplacesDemon)
library(truncnorm)
library(hydra)

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

# Complete update, after optimization
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
  delnew[id] = delold[id] = 1
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

main_bhmds <- function(n, d, D, V0, X0, sigg0, a, alpha, maxiter, constant, betas, burnin, kappa) {
  kappa = kappa
  tracedelta = matrix(rep(0, n^2), ncol = n)
  trace_12 <- rep(0, maxiter - burnin)
  trace_sigma <- rep(0, maxiter - burnin)
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
      temp = compute_delta(X, kappa)
      if (SSRmin > SSRnew) {
        tracedelta = temp
        SSRmin = SSRnew
      }
      trace_12[i - burnin] = temp[1, 2]
      trace_sigma[i - burnin] = sigma2
    }
    SSRold = SSRnew
    if (i == burnin) {
      SSRmin = compute_SSR_mat(D, V, kappa)
      tracedelta = compute_delta(X, kappa)
    }
    print(i)
  }
  # return trace for delta_12, sigma, and Bayesian estimate of delta_ij
  result <- c(I(list(tracedelta)), I(list(trace_12)), I(list(trace_sigma)))
  return (result)
}

# Sample size n
n = 50

# Hyperbolic dimension
d = 2

# Observational error 
sigma = 1

# Curvature
kappa = 1

mu0 <- c(1, rep(0, d))
m = choose(n, 2)

# Generate true matrix
V <- matrix(rep(0, n * d), ncol = d)
for (i in 1:n) {
  V[i,] <- rmvnorm(1, rep(0, d), diag(rep(3, d)))
}
delta <- compute_delta(VtoX(V), kappa)

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

result <- main_bhmds(n, d, D, V0, X0, sigg0, a, alpha, maxiter, constant, beta0, burnin, kappa)
