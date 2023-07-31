library(igraph)
library(mvtnorm)
library(LaplacesDemon)
library(truncnorm)
library(hydra)
dat <- read.csv("dist.csv", header = F)
data <- as.matrix(dat)
n = dim(data)[1]
temp <- data[upper.tri(data)]
temp_mat <- matrix(rep(0, n^2), ncol = n)
temp_mat[upper.tri(temp_mat)] = temp
temp_mat = temp_mat + t(temp_mat)
D = temp_mat / 10
(isSymmetric.matrix(D))
# hydra <- hydraPlus(data, dim = 2, curvature = 1, control = list(return.lorentz = T), maxit = 3000)
# saveRDS(hydra, "hydra_10.RData")
hydra <- readRDS("hydra_10_5.RData")

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

pilot_update_vvec <- function (D, V, X, id, sigma2, constant, lambdas, nearest, random_sample, strata) {
  Vnew = V
  stepsize = sqrt(sigma2 * constant / (n - 1))
  # proposal 
  for (i in 1:d) {
    Vnew[id, i] = Vnew[id, i] + rnorm(1, 0, stepsize) #true_proposal[i]
  }
  sigma = sqrt(sigma2)
  sample_size = param_r * mean_length
  
  # Compute the updated vector
  Xnew = X
  lornew = sqrt(Lorentz(c(0, Vnew[id,]), c(0, Vnew[id,])))
  Xnew[id,] = cosh(lornew) * mu0 + sinh(lornew) / lornew * c(0, Vnew[id,])
  xold = X[id,]; xnew = Xnew[id,]
  xold[1] = -xold[1]; xnew[1] = -xnew[1]
  Q1new = Q1old = Q2new = Q2old = Q3new = Q3old = 0
  tempnew = 0
  tempold = 0
  
  ## Only proportion that is not sampled
  Q2new = sum(Vnew[id,]^2 * (1/lambdas))
  Q2old = sum(V[id,]^2 * (1/lambdas))
  
  if (length(nearest[[id]]) != 0) {
    # Compute the loglik of nearest points
    if (length(nearest[[id]]) == 1) {
      delnew = -xnew %*% t(t(Xnew[nearest[[id]],]))
      delold = -xold %*% t(t(X[nearest[[id]],]))
    } else {
      delnew = -xnew %*% t(Xnew[nearest[[id]],])
      delold = -xold %*% t(X[nearest[[id]],]) 
    }
    delnew = acosh(delnew) * (1 / sqrt(kappa))
    delold = acosh(delold) * (1 / sqrt(kappa))
    Q1new = sum((D[id, nearest[[id]]] - delnew)^2) / sigma2
    Q1old = sum((D[id, nearest[[id]]] - delold)^2) / sigma2
    Q3new = sum(log(pnorm(delnew / sigma)))
    Q3old = sum(log(pnorm(delold / sigma)))
    
    # Sample in all strata
    delnew_sample = -xnew %*% t(Xnew[random_sample[id,],])
    delold_sample = -xold %*% t(X[random_sample[id,],])
    delnew_sample = acosh(delnew_sample) * (1 / sqrt(kappa))
    delold_sample = acosh(delold_sample) * (1 / sqrt(kappa))
    Q1new_sample = (D[id, random_sample[id,]] - delnew_sample)^2 / sigma2
    Q1old_sample = (D[id, random_sample[id,]] - delold_sample)^2 / sigma2
    Q3new_sample = log(pnorm(delnew_sample / sigma))
    Q3old_sample = log(pnorm(delold_sample / sigma)) 
  } else {
    # Sample in all strata
    Q1new = Q1old = Q3new = Q3old = 0
    delnew_sample = -xnew %*% t(Xnew[random_sample[id,],])
    delold_sample = -xold %*% t(X[random_sample[id,],])
    delnew_sample = acosh(delnew_sample) * (1 / sqrt(kappa))
    delold_sample = acosh(delold_sample) * (1 / sqrt(kappa))
    Q1new_sample = (D[id, random_sample[id,]] - delnew_sample)^2 / sigma2
    Q1old_sample = (D[id, random_sample[id,]] - delold_sample)^2 / sigma2
    Q3new_sample = log(pnorm(delnew_sample / sigma))
    Q3old_sample = log(pnorm(delold_sample / sigma)) 
  }
  
  # Compute the likelihod change in each strata
  strata_vec = strata[[id]]
  Qnew_vec = -Q1new_sample/2 - Q3new_sample
  Qold_vec = -Q1old_sample/2 - Q3old_sample
  Q_diff_vec = Qnew_vec - Qold_vec
  w_vec_temp = abs(Q_diff_vec / sum(Q_diff_vec))
  w_vec = rep(0, length(strata_vec))
  for (i in 1:length(strata_vec)) {
    d_temp = D[id, random_sample[id,]]
    w_vec[i] = sum(w_vec_temp[as.vector(which(d_temp == strata_vec[i]))])
  }
  
  Q1new_sample = sum(Q1new_sample)
  Q1old_sample = sum(Q1old_sample)
  Q3new_sample = sum(Q3new_sample)
  Q3old_sample = sum(Q3old_sample)
  
  # fnew = (-Q1new / 2 - Q3new) * n / sample_size + (-Q2new / 2)
  # fold = (-Q1old / 2 - Q3old) * n / sample_size + (-Q2old / 2)
  
  f_influential_new = -(Q1new+Q2new)/2 - Q3new
  f_influential_old = -(Q1old+Q2old)/2 - Q3old
  f_estimate_new = (-Q1new_sample/2 - Q3new_sample) * (n - length(nearest[[id]]) - 1) / sample_size
  f_estimate_old = (-Q1old_sample/2 - Q3old_sample) * (n - length(nearest[[id]]) - 1) / sample_size
  fnew = f_influential_new + f_estimate_new
  fold = f_influential_old + f_estimate_old
  fratio = min(exp(fnew - fold), 1)
  if (runif(1) <= fratio) {
    return(c(Vnew[id,], Xnew[id,], w_vec))
  } else {
    return(c(V[id,], X[id,], w_vec))
  }
}

main_bhmds <- function(n, d, D, V0, X0, sigg0, a, alpha, maxiter, constant, betas, burnin, kappa, nearest, random_sample, strata) {
  kappa = kappa
  tracedelta = matrix(rep(0, n^2), ncol = n)
  m = choose(n, 2)
  V = V0
  X = X0
  w_list <- list()
  for (i in 1:n) {
    w_list[[i]] = rep(0, length(strata[[i]]))
  }
  # Vnew = matrix(rep(0, n * d), ncol = d)
  
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
      # print(j)
      tmprow = pilot_update_vvec(D, V, X, j, sigma2, constant, lambdas, nearest, random_sample, strata)
      #tmprow = update_vvec(D, V, X, j, sigma2, constant, lambdas)
      V[j,] = tmprow[1:d]
      X[j,] = tmprow[(d + 1): (2*d + 1)]
      w_list[[j]] = w_list[[j]] + tmprow[-1:-(2*d + 1)]
    }
    SSRnew = compute_SSR_mat(D, V, kappa)
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
        temp = compute_delta(X, kappa)
        tracedelta = temp
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
  return (w_list)
}

length_res <- rep(0, n)
for (i in 1:n) {
  length_res[i] = length(which(D[i,] < 7)) - 1
}

(mean_length = round(mean(length_res)))
param_r = 10
random_sample <- matrix(rep(0, param_r * mean_length * n), nrow = n)
nearest <- list()
strata <- list()
for (i in 1:n) {
  print(i)
  pool = as.vector(which(D[i, ] < 7))
  pool = pool[pool != i]
  nearest[[i]] <- pool
  if (length(pool) != 0) {
    random_sample[i,] = sample(c(1:n)[-c(pool, i)], size = mean_length * param_r, replace = F)
    strata[[i]] = sort(unique(D[i, random_sample[i,]]))
  } else {
    random_sample[i,] = sample(c(1:n)[-i], size = mean_length * param_r, replace = F)
    strata[[i]] = sort(unique(D[i, random_sample[i,]]))
  }
}

d = 5
mu0 <- c(1, rep(0, d))
m = choose(n, 2)

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
V0 = Vini[,-1]
X0 <- VtoX(V0)
maxiter = 3000
kappa = 1
delta0 <- compute_delta(VtoX(V0), kappa)
sigg0 = compute_SSR(D, delta0) / m
a = 5
b = (a - 1) * sigg0
alpha = 0.5
beta0 = apply(V0, 2, var) / 2
constant = 1
burnin = 1000

w_list <- main_bhmds(n, d, D, V0, X0, sigg0, a, alpha, maxiter, constant, beta0, burnin, kappa, nearest, random_sample, strata) 
rescale_w_list <- w_list
for (i in 1:n) {
  rescale_w_list[[i]] = rescale_w_list[[i]] / (maxiter - burnin)  
}
write.csv(random_sample, "rand_sample.csv")
saveRDS(w_list, "w_list.RData")
w_list <- readRDS("w_list.Rdata")
temp <- as.numeric(as.matrix(read.csv("rand_sample.csv")[,-1]))
random_sample <- matrix(temp, ncol = 310)
n = dim(random_sample)[1]
rescale_w_list <- w_list
for (i in 1:n) {
  rescale_w_list[[i]] = rescale_w_list[[i]] / (maxiter - burnin)  
}
hist(D)
strata_count <- list()
for (i in 1:n) {
  tempsamp <- random_sample[i,]
  tempstrata <- c()
  for (j in 1:17) {
    tempstrata <- c(tempstrata, sum((D[i, tempsamp] > (j + 6)) & (D[i, tempsamp] <= (j + 7))))
  }
  tempstrata <- c(tempstrata, sum(D[i, tempsamp] > 24))
  strata_count[[i]] <- tempstrata
}
saveRDS(strata_count, "strata_count.Rdata")
strata_list <- list()
for (i in 1:n) {
  tempstrata <- strata_count[[i]]
  strata_list[[i]] <- which(tempstrata != 0) + 6
}
saveRDS(strata_list, "strata_list.Rdata")
strata <- matrix(rep(0, param_r * mean_length * n), nrow = n)
for (i in 1:n) {
  tempsamp <- random_sample[i,]
  for (j in 1: (param_r * mean_length)) {
    res <- floor(D[i, tempsamp[j]])
    if (res != 25) {
      strata[i, j] = res
    } else {
      strata = 24
    }
  }
}
write.csv(strata, "strata.csv")
