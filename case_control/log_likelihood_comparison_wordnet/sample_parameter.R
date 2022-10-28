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

# Compute true likelihood

# Testing set
# V_true = matrix(c(1:4), ncol = 2)
# D = V_true + 0.3
# d = n = 2
# sigma2 = 4
# constant = 1
# lambdas = c(1, 1)
# V = V_true
# true_proposal = c(-1.034, 1.59)
# id = 1
# X_true = VtoX(V_true)
# X = X_true
# mu0 = c(1, 0, 0)
# kappa = 1

# Complete update, after optimization
update_vvec <- function (D, V, X, id, sigma2, constant, lambdas) {
  Vnew = V
  stepsize = sqrt(sigma2 * constant / (n - 1))
  # proposal 
  for (i in 1:d) {
    Vnew[id, i] = Vnew[id, i] + rnorm(1, 0, stepsize) #true_proposal[i]
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

main_bhmds <- function(n, d, D, V0, X0, sigg0, a, alpha, maxiter, constant, betas, burnin, kappa) {
  kappa = kappa
  tracesigma <- rep(0, maxiter - burnin)
  trace_12 <- trace_34 <- trace_56 <- trace_78 <- trace_910 <- rep(0, maxiter - burnin)
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
  
  varalpha = 0.0
  varbeta  = 0.0
  varvar = 0.0
  varratio = 0.0
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
    sigtmp = sigma2 + rnorm(1, 0, sqrt(constant*varvar))
    if (sigtmp > 0){ 
      varratio = dinvgamma(sigtmp,varalpha,varbeta)/dinvgamma(sigma2,varalpha,varbeta)
      # varratio = (sigma2/sigtmp)^(varalpha + 1) * exp(-varbeta/sigtmp  + varbeta/sigma2)
      if (varratio > 1){
        varratio = 1.0
      }
      if (runif(1) <= varratio){
        sigma2 = sigtmp
      }
    }
    if (i > burnin) {
      temp = compute_delta(X, kappa)
      tracesigma[i - burnin] = sigma2
      trace_12[i - burnin] = temp[1, 2]
      temp = compute_delta(X, kappa)
      file_name_delta = paste("sample_delta_", (i - burnin), ".csv", sep = "")
      file_name_v = paste("sample_v_", (i - burnin), ".csv", sep = "")
      file_name_param = paste("sample_param_", (i - burnin), ".csv", sep = "")
      write.csv(temp, file_name_delta)
      write.csv(V, file_name_v)
      write.csv(c(sigma2, lambdas), file_name_param)
      if (SSRmin > SSRnew) {
        tracedelta = temp
        SSRmin = SSRnew
        # print(SSRmin)
      } 
    }
    SSRold = SSRnew
    if (i == burnin) {
      SSRmin = compute_SSR_mat(D, V, kappa)
      tracedelta = compute_delta(X, kappa)
    }
    print(i)
  }
  return (rbind(trace_12, tracesigma))
}

n = 1141
d = 2
mat <- as.matrix(read.csv("graph_wordnet_mammal.csv")[,-1])
mat = mat + t(mat)
colnames(mat) <- rownames(mat) <- 1:1141
g <- graph_from_adjacency_matrix(mat)
dist <- distances(g)
D = dist
d = 2
mu0 <- c(1, rep(0, d))
m = choose(n, 2)
kappa = 2.06

# Initialization
V0 <- as.matrix(read.csv("hydrav0.csv")[,-1])
X0 <- VtoX(V0)
maxiter = 100
delta0 <- compute_delta(VtoX(V0), kappa)
sigg0 = compute_SSR(D, delta0) / m
a = 5
b = (a - 1) * sigg0
alpha = 0.5
beta0 = apply(V0, 2, var) / 2
constant = 1
burnin = 50

result <- main_bhmds(n, d, D, V0, X0, sigg0, a, alpha, maxiter, constant, beta0, burnin, kappa)
