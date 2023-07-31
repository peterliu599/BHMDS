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

pilot_update_vvec <- function (D, V, X, id, sigma2, constant, lambdas, nearest, random_sample, coef) {
  Vnew = V
  
  stepsize = sqrt(sigma2 * constant / (n - 1))
  # proposal 
  for (i in 1:d) {
    Vnew[id, i] = Vnew[id, i] + rnorm(1, 0, stepsize) 
  }
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
  
  if (length(nearest) != 0) {
    # Compute the loglik of nearest points
    if (length(nearest) == 1) {
      delnew = -xnew %*% t(t(Xnew[nearest,]))
      delold = -xold %*% t(t(X[nearest,]))
    } else {
      delnew = -xnew %*% t(Xnew[nearest,])
      delold = -xold %*% t(X[nearest,]) 
    }
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
  } else {
    # Sample in all strata
    Q1new = Q1old = Q3new = Q3old = 0
    delnew_sample = -xnew %*% t(Xnew[random_sample,])
    delold_sample = -xold %*% t(X[random_sample,])
    delnew_sample = acosh(delnew_sample) * (1 / sqrt(kappa))
    delold_sample = acosh(delold_sample) * (1 / sqrt(kappa))
    Q1new_sample = sum((D[id, random_sample] - delnew_sample)^2 / sigma2 * coef)
    Q1old_sample = sum((D[id, random_sample] - delold_sample)^2 / sigma2 * coef)
    Q3new_sample = sum(log(pnorm(delnew_sample / sigma)) * coef)
    Q3old_sample = sum(log(pnorm(delold_sample / sigma)) * coef)
  }
  f_influential_new = -(Q1new+Q2new)/2 - Q3new
  f_influential_old = -(Q1old+Q2old)/2 - Q3old
  f_estimate_new = (-Q1new_sample/2 - Q3new_sample)
  f_estimate_old = (-Q1old_sample/2 - Q3old_sample)
  fnew = f_influential_new + f_estimate_new
  fold = f_influential_old + f_estimate_old
  fratio = min(exp(fnew - fold), 1)
  if (runif(1) <= fratio) {
    return(c(Vnew[id,], Xnew[id,]))
  } else {
    return(c(V[id,], X[id,]))
  }
}


main_bhmds <- function(n, d, D, V0, X0, sigg0, a, alpha, maxiter, constant, 
                       betas, burnin, kappa, nearest, rand_sample, coef, grouplabel) {
  kappa = kappa
  
  # Trace record for delta_12 and sigma
  # trace12 <- tracesigma <- rep(0, maxiter - burnin)
  # tracehdist <- matrix(rep(0, (n * (maxiter - burnin))), ncol = n)
  # tracedelta = matrix(rep(0, n^2), ncol = n)
  # sum_count = matrix(rep(0, ((maxiter - burnin) * 15)), ncol = 15)
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
  iter = 0
  
  grouplength <- rep(0, 15)
  for (i in 1:15) {
    grouplength[i] = length(which(grouplabel == i))
  }
  
  for (i in 1:maxiter) {
    for (j in 1:d) {
      vecs[j] = var(V[,j]) * n
    }
    for (j in 1:d) {
      lambdas[j] = rinvgamma(1, alpha + n/2, betas[j] + vecs[j] / 2)
    }
    acceptrate <- rep(0, n)
    for (j in 1:n) {
      tmprow = pilot_update_vvec(D, V, X, j, sigma2, constant, lambdas, nearest[[j]], rand_sample[[j]], coef[[j]])
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
      flag = i - burnin
      temp = compute_delta(X, kappa)
      if (SSRmin > SSRnew) {
        traceV = V
        tracedelta = temp
        SSRmin = SSRnew
        iter = i
      }
      # record delta_12
      write.csv(temp[1, 2], paste("delta_12_", flag, ".csv", sep = ""))
      # print("done")
      # record sigma
      write.csv(sigma2, paste("sigma2_", flag, ".csv", sep = ""))
      # print("done")
      # record d_h(X)
      tracehdist = acosh(X[,1])
      write.csv(tracehdist, paste("tracehdist_", flag, ".csv", sep = ""))
      # print("done")
      # record cluster distance
      tracecdist <- matrix(rep(0, 15^2), ncol = 15)
      for (j in 1:15) {
        for (k in 1:15) {
          if (j < k) {
            tracecdist[j, k] = sum(temp[which(grouplabel == j), which(grouplabel == k)]) / 
              (grouplength[j] * grouplength[k])
            # print(c(j, k))
          }
        }
      }
      write.csv(tracecdist, paste("tracecdist_", flag, ".csv", sep = ""))
    }
    SSRold = SSRnew
    if (i == burnin) {
      SSRmin = compute_SSR_mat(D, V, kappa)
      tracedelta = compute_delta(X, kappa)
      traceV = V
      iter = i
    }
    print(i)
  }
  result <- c(I(list(tracedelta)), I(list(traceV)), iter)
  return (result)
}

dat <- read.csv("dist.csv", header = F)
data <- as.matrix(dat)
n = dim(data)[1]
temp <- data[upper.tri(data)]
temp_mat <- matrix(rep(0, n^2), ncol = n)
temp_mat[upper.tri(temp_mat)] = temp
temp_mat = temp_mat + t(temp_mat)
D = temp_mat / 10
(isSymmetric.matrix(D))
hydra <- readRDS("hydra_10_5.RData")

d = 5
mu0 <- c(1, rep(0, d))
m = choose(n, 2)
kappa = 1

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
delta0 <- compute_delta(VtoX(V0), kappa)
sigg0 = compute_SSR(D, delta0) / m
a = 5
b = (a - 1) * sigg0
alpha = 0.5
beta0 = apply(V0, 2, var) / 2
constant = 1

nearest <- list()
for (i in 1:n) {
  pool = as.vector(which(D[i, ] < 7))
  pool = pool[pool != i]
  nearest[[i]] <- pool
}

length_nearest <- rep(0, n)
for (i in 1:n) {
  length_nearest[i] = length(nearest[[i]])
}

mean_length_nearest = round(mean(length_nearest))
param_r = 10
# Weight list
w_list <- readRDS("w_list.RData")
for (i in 1:n) {
  w_list[[i]] = w_list[[i]] / 2000  
}
# Strata list: what are the strata?
strata_list <- readRDS("strata_list.RData")
# strata id for each sample element, n * p
strata <- as.numeric(as.matrix(read.csv("strata.csv")[,-1]))
strata <- matrix(strata, ncol = 310)
# How long is each strata, population value, BIG N
strata_length <- list()
for (i in 1:n) {
  tempD <- D[i, ]
  tempstrata <- c()
  for (j in 1:17) {
    tempstrata <- c(tempstrata, sum((tempD > (j + 6)) & (tempD <= (j + 7))))
  } 
  tempstrata <- c(tempstrata, sum(tempD > 24))
  strata_length[[i]] <- tempstrata
}
# How many sample for each strata, SMALL N
strata_count <- readRDS("strata_count.RData")
rescale_w_list <- list()
for (i in 1:n) {
  tempres <- rep(0, length(strata_count[[i]]))
  for (j in 1: (param_r * mean_length_nearest)) {
    flag = strata[i, j] - 6
    tempres[flag] = tempres[flag] + w_list[[i]][j]
  }
  rescale_w_list[[i]] = tempres
}
sample_list <- list()
sample_length <- rep(0, n)
for (i in 1:n) {
  tempsum = sum(w_list[[i]])
  sample_list[[i]] = ceiling(param_r * mean_length_nearest * rescale_w_list[[i]] / tempsum)
  for (j in 1:length(sample_list[[i]])) {
    d_count = strata_length[[i]][j]
    if (sample_list[[i]][j] > d_count) {
      sample_list[[i]][j] = d_count
    }
  }
  sample_length[i] = sum(sample_list[[i]])
}
tempsum = 0
for (i in 1:n) {
  tempsum = tempsum + (length(sample_list[[i]]) == length(sort(unique(strata[i,]))))
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
    print(c(i, j))
    dtemp = j + 6 #sort(unique(strata[i,]))[j] # sample d
    if (is.na(dtemp) == F) {
      if (dtemp != 24) {
        tag <- (D[i,] > dtemp) & (D[i,] <= dtemp + 1)
        pool = as.vector(which(tag))#D[i, tag]# sample d pool
      } else {
        tag <- (D[i, ] > 24)
        pool = as.vector(which(tag))#D[i, tag]
      }
      if (length(pool) != 0) {
        tempsam = c(tempsam, sample(pool, size = sample_list[[i]][j], replace = F))
        coef_id = strata_length[[i]][j] / sample_list[[i]][j]
        tempcoef = c(tempcoef, rep(coef_id, sample_list[[i]][j])) 
      }
    }
  }
  rand_sample[[i]] = tempsam
  coef[[i]] = tempcoef
}

grouplabel_temp <- as.matrix(read.csv("GroupLabelTrim3.csv", header = F))
label <- unique(t(grouplabel_temp))
grouplabel <- rep(0, n)
for (i in 1:n) {
  for (j in 1:15) {
    if (grouplabel_temp[i] == label[j]) {
      grouplabel[i] = j
    }
  }
}

maxiter = 25000
burnin = 10000

result <- main_bhmds(n, d, D, V0, X0, sigg0, a, alpha, maxiter,
                     constant, beta0, burnin, kappa, nearest, rand_sample, coef, grouplabel)

saveRDS(result, "final_result.RData")

## Plots

# tracecdist <- matrix(rep(0, 15000 * 105), ncol = 105)
# for (i in 1:15000) {
#   file = paste("tracecdist_", i , ".csv", sep = "")
#   temp = as.matrix(read.csv(file, header = T)[, -1])
#   tracecdist[i, ] = temp[upper.tri(temp)]
#   print(i)
# }
# 
# med_tracecdist <- matrix(rep(0, 225), ncol = 15)
# med_temp <- rep(0, 105)
# for (i in 1:105) {
#   # plot(tracecdist[,i], type ="l")
#   med_temp[i] = median(tracecdist[,i])
# }
# med_tracecdist[upper.tri(med_tracecdist)] = med_temp
# med_tracecdist = med_tracecdist + t(med_tracecdist)
# colnames(med_tracecdist) = label
# rownames(med_tracecdist) = label
# label_3 <- c("blood neoplasm cell line","non leukemic blood neoplasm", "leukemia","normal blood", 
#              "blood non neoplastic disease","nervous system neoplasm","solid tissue non neoplastic disease",
#              "normal solid tissue","solid tissue neoplasm cell line", "non neoplastic cell line",
#              "non breast carcinoma","breast cancer", "germ cell neoplasm","sarcoma", "other neoplasm")
# temp_med_tracecdist <- med_tracecdist[label_3, label_3]
# diag(temp_med_tracecdist) = 10
# df <- melt(temp_med_tracecdist)
# colnames(df) <- c("Cell Type 1", "Cell Type 2", "Pseudotime")
# ggplot(df, aes(y = `Cell Type 1`, x = `Cell Type 2`, fill = Pseudotime)) +
#   geom_tile() + coord_fixed() +  scale_fill_distiller(palette = "OrRd", direction = 1, name ="Cluster distance") +
#   theme(legend.title = element_text(size=20), legend.text = element_text(size=15),
#         axis.text=element_text(size=16), axis.title=element_text(size=1)) + 
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))  +
#   xlab("") + ylab("")
# 
# tracehdist <- matrix(rep(0, 15000 * n), ncol = n)
# for (i in 1:15000) {
#   file = paste("tracehdist_", i , ".csv", sep = "")
#   temp = read.csv(file, header = T)[, 2]
#   tracehdist[i, ] = temp
#   print(i)
# }
# 
# hdist_list <- vector(mode = "list", length = 15)
# for (i in 1:n) {
#   citemp <- ci(tracehdist[,i], method = "ETI")
#   templab <- which((tracehdist[,i] > citemp$CI_low) & (tracehdist[,i] < citemp$CI_high))
#   hdist_list[[grouplabel[i]]] <- c(hdist_list[[grouplabel[i]]], tracehdist[templab, i])
# }
# 
# s = 10000
# sam <- ord <- matrix(rep(0, 15 * s), ncol = 15)
# for (i in 1:15) {
#   sam[,i] = sample(hdist_list[[i]], s, replace = T)
# }
# for (i in 1:s) {
#   ord[i,] = order(sam[i,])
# }
# res <- matrix(rep(0, 225), ncol = 15)
# for (i in 1:15){
#   tempsam <- ord[,i]
#   for (j in 1:15) {
#     res[j, i] = length(which(tempsam == j))
#   }
# }
# res <- res / 10000
# rownames(res) <- label
# colnames(res) <- c("1st", "2nd", "3rd", "4th","5th","6th","7th","8th","9th","10th","11th","12th","13th","14th","15th")
# label_3 <- c("solid tissue non neoplastic disease", "normal solid tissue","normal blood", "leukemia",
#              "blood non neoplastic disease", "blood neoplasm cell line",
#              "nervous system neoplasm", "non neoplastic cell line", 
#              "solid tissue neoplasm cell line" ,
#              "sarcoma","non leukemic blood neoplasm", 
#              "non breast carcinoma", "other neoplasm",
#              "breast cancer" ,"germ cell neoplasm")
# res <- res[label_3, ]
# res <- round(res, 3)
# df <- melt(res)
# colnames(df) <- c("Cell Type", "Rank", "value")
# ggplot(df, aes(y = `Cell Type`, x = Rank, fill = value)) +
#   geom_tile(color = "white",
#             lwd = 1.5,
#             linetype = 1) + scale_fill_distiller(palette = "YlOrRd", direction = 1, name = "Frequency")+
#   coord_fixed() + geom_text(aes(label = value), color = "black", size = 3.5) +
#   theme(legend.title = element_text(size=25), legend.text = element_text(size=20),
#         axis.text=element_text(size = 18), axis.title=element_text(size=25))