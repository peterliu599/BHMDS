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

pilot_update_vvec <- function (D, V, X, tempvec, id, sigma2, constant, lambdas, nearest, random_sample, coef) {
  Vnew = V
  Xnew = X
  Vnew[id, ] = tempvec
  stepsize = sqrt(sigma2 * constant / (n - 1))
  sigma = sqrt(sigma2)
  # sample_size = param_r * mean_length
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
  
  # Compute the loglik of nearest points
  delnew = -xnew %*% t(Xnew[nearest,])
  delold = -xold %*% t(X[nearest,])
  delnew = acosh(delnew) * (1 / sqrt(kappa))
  delold = acosh(delold) * (1 / sqrt(kappa))
  Q1new = sum((D[id, nearest] - delnew)^2) / sigma2
  Q1old = sum((D[id, nearest] - delold)^2) / sigma2
  Q3new = sum(log(pnorm(delnew / sigma)))
  Q3old = sum(log(pnorm(delold / sigma)))
  
  # Sample in all strata
  delnew_sample = -xnew %*% t(Xnew[random_sample,])
  delold_sample = -xold %*% t(X[random_sample,])
  delnew_sample = acosh(delnew_sample) * (1 / sqrt(kappa))
  delold_sample = acosh(delold_sample) * (1 / sqrt(kappa))
  Q1new_sample = sum((D[id, random_sample] - delnew_sample)^2 / sigma2 * coef)
  Q1old_sample = sum((D[id, random_sample] - delold_sample)^2 / sigma2 * coef)
  Q3new_sample = sum(log(pnorm(delnew_sample / sigma)) * coef)
  Q3old_sample = sum(log(pnorm(delold_sample / sigma)) * coef)
  
  f_influential_new = -(Q1new+Q2new)/2 - Q3new
  f_influential_old = -(Q1old+Q2old)/2 - Q3old
  f_estimate_new = (-Q1new_sample/2 - Q3new_sample)
  f_estimate_old = (-Q1old_sample/2 - Q3old_sample)
  fnew = f_influential_new + f_estimate_new
  fold = f_influential_old + f_estimate_old
  
  return((fnew - fold))
}

# Complete update, after optimization
update_vvec <- function (D, V, X, tempvec, id, sigma2, constant, lambdas) {
  Vnew = V
  Xnew = X
  Vnew[id, ] = tempvec
  sigma = sqrt(sigma2)
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
  return((fnew - fold))
}

# Using the intialization values are defaults

# Read in subtree
n = 1141
d = 2
mat <- as.matrix(read.csv("graph_wordnet_mammal.csv")[,-1])
mat = mat + t(mat)
colnames(mat) <- rownames(mat) <- 1:1141
g <- graph_from_adjacency_matrix(mat)
dist <- distances(g)
D = dist

mu0 <- c(1, rep(0, d))
m = choose(n, 2)
kappa = 2.06

# Saved hydraPlus result
V0 <- as.matrix(read.csv("hydrav0.csv")[,-1])
X0 <- VtoX(V0)
delta0 <- compute_delta(VtoX(V0), kappa)
sigg0 = compute_SSR(D, delta0) / m
a = 10
b = (a - 1) * sigg0
alpha = 0.5
beta0 = apply(V0, 2, var) / 2
constant = 1

nearest <- list()
for (i in 1:1141) {
  pool = c(as.vector(which(D[i,] == 1)), as.vector(which(D[i,] == 2)))
  nearest[[i]] <- pool
}

length_nearest <- rep(0, 1141)
for (i in 1:1141) {
  length_nearest[i] = length(which(D[i,] == 1)) + length(which(D[i,] == 2))
}
mean_length_nearest = round(mean(length_nearest))
param_r = 5
# Weight list 
w_list <- readRDS("weight_5.RData")
# Strata list
strata_list <- readRDS("strata_5.RData")
# How long is each strata, population value, BIG N
strata_length <- list()
# How many sample for each strata, SMALL N
sample_list <- list()
# The length for each i sample, no longer need.
sample_length <- rep(0, n)
for (i in 1:n) {
  sample_list[[i]] = strata_length[[i]] = rep(0, length(w_list[[i]]))
}
for (i in 1:n) {
  tempsum = sum(w_list[[i]])
  sample_list[[i]] = ceiling(param_r * mean_length_nearest * w_list[[i]] / tempsum)
  for (j in 1:length(sample_list[[i]])) {
    strata_length[[i]][j] = d_count = length(which(D[i, ] == strata_list[[i]][j]))
    if (sample_list[[i]][j] > d_count) {
      sample_list[[i]][j] = d_count
    }
  }
  sample_length[i] = sum(sample_list[[i]])
}
# Random samples
rand_sample <- list()
# Coef for the log-likelihood
coef <- list()
for (i in 1:n) {
  rand_sample[[i]] = coef[[i]] = rep(0, sample_length[i])
}
for (i in 1:n) {
  tempsam = tempcoef = c()
  for (j in 1:length(sample_list[[i]])) {
    dtemp = strata_list[[i]][j] # sample d
    pool = as.vector(which(D[i,] == dtemp)) # sample d pool
    tempsam = c(tempsam, sample(pool, size = sample_list[[i]][j], replace = F))
    coef_id = strata_length[[i]][j] / sample_list[[i]][j]
    tempcoef = c(tempcoef, rep(coef_id, sample_list[[i]][j]))
  }
  rand_sample[[i]] = tempsam
  coef[[i]] = tempcoef
}


V = V0
X = X0
sigma2 = sigg0
vecs <- rep(0, d)
lambdas = rep(0, d)
betas = beta0
for (i in 1:d) {
  vecs[i] = var(V[,i]) * n
}
for (i in 1:d) {
  lambdas[i] = rinvgamma(1, alpha + n/2, betas[i] + vecs[i] / 2)
}

likelihood_result <- matrix(rep(0, 10000 * 2), ncol = 2)

for (i in 1:100) {
  file_name_delta = paste("sample_delta_", i, ".csv", sep = "")
  file_name_v = paste("sample_v_", i, ".csv", sep = "")
  file_name_param = paste("sample_param_", i, ".csv", sep = "")
  temp = as.matrix(read.csv(file_name_delta)[,-1])
  data = matrix(rep(0, n^2), ncol = n)
  data[upper.tri(data)] = temp[upper.tri(temp)]
  data = data + t(data)
  V = as.matrix(read.csv(file_name_v)[,-1])
  X = VtoX(V)
  object_id = sample(1:n, 1)
  param = as.matrix(read.csv(file_name_param)[,-1])
  sigma2 = param[1]
  lambdas = param[-1]
  stepsize = sqrt(sigma2 * constant / (n - 1))
  id_vec <- sample(1:n, 100)
  for (k in 1:100) {
    object_id = id_vec[k]
    tempvec <- V[object_id,]
    for (j in 1:d) {
      tempvec[j] = tempvec[j] + rnorm(1, 0, stepsize)
    }
    likelihood_result[((i - 1) * 100 + k), 1] = update_vvec(D, V, X, tempvec, object_id, sigma2, constant, lambdas)
    likelihood_result[((i - 1) * 100 + k), 2] = pilot_update_vvec(D, V, X, tempvec, object_id, sigma2, constant, lambdas, nearest[[object_id]], rand_sample[[object_id]], coef[[object_id]])
  }
  print(i)
}

# write.csv(likelihood_result, "likelihood_result.csv")

likelihood_result <- as.matrix(read.csv("likelihood_result.csv")[,-1])
y = likelihood_result[,1]
x = likelihood_result[,2]

cor(y, x)
library(tidyverse)
dat <- data.frame(
  x = likelihood_result[,1],
  y = likelihood_result[,2]
)
ggplot(data = dat, aes(y = y, x = x)) + geom_point() + 
  geom_abline(slope = 1, intercept = 0, col = "red") + 
  xlab("Full Log-lik") +
  ylab("Approx Log-lik") + 
  theme(text = element_text(size=20), axis.text=element_text(size=15)) + 
  ggtitle("Full Log-lik vs. Case-controled  Log-lik") + theme(plot.title = element_text(hjust = 0.5)) 

