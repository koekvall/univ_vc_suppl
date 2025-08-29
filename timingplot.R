## 2matrices timing
set.seed(145)
library(dplyr)
library(ggplot2)
library(patchwork)

TOL <- 1e-8
N <- 30
nseq <- c(200, 500, 1000, 2000, 3000)

# data simulation with upper triangular factor of the Cholesky decomposition of K
simData <- function(h21, h22, sigma2=1, R1, R2, V) { # R are sqrt(eigen of K)
  n <- length(R1)
  e = rnorm(n, mean = 0, sd = sqrt(sigma2 * (1-h21-h22)))
  u1 = rnorm(n, mean = 0, sd = sqrt(sigma2 * h21)) # mvrnorm(1, mu = rep(0, n), Sigma = sigma2 * h21 * K1)
  u2 = rnorm(n, mean = 0, sd = sqrt(sigma2 * h22)) # mvrnorm(1, mu = rep(0, n), Sigma = sigma2 * h22 * K2)
  y = R1 * u1 + R2 * u2 + e              #crossprod(R2, u2) + e
  return(V%*%y)
}

alpha <- 0.05
h22 <- 0.2
sigma2 <- 1
h21 <- 0
M = 2



time_naive = matrix(0, nrow = N, ncol = length(nseq))
time_test0 = matrix(0, nrow = N, ncol = length(nseq))
time_ortho = matrix(0, nrow = N, ncol = length(nseq))

test_naive = matrix(0, nrow = N, ncol = length(nseq))
test_test0 = matrix(0, nrow = N, ncol = length(nseq))
test_ortho = matrix(0, nrow = N, ncol = length(nseq))

for (i in 1:length(nseq)) {
  n <- nseq[i]
  K <- 0.5^abs(outer(1:n, 1:n, '-'))
  eoK <- eigen(K)
  V <- eoK$vectors
  lambda <- eoK$val

  n_each <- floor(n / (M + 1))
  idx_mat <- matrix(sample(n, size = M * n_each), ncol = M)
  Klist <- list()
  for(ii in 1:M){
    Klist[[ii]] <- V[, idx_mat[, ii]] %*% (lambda[idx_mat[, ii]] * t(V[, idx_mat[, ii]]))
  }
  #Rlist <- lapply(Klist, function(K) { sqrt(pmax(diag(t(V) %*% K %*% V), 0))})
  R1 <- sqrt(pmax(diag(t(V) %*% Klist[[1]] %*% V), 0))
  R2 <- sqrt(pmax(diag(t(V) %*% Klist[[2]] %*% V), 0))

  # eoK <- eigen(K)
  # V <- eoK$vec
  # lambda <- eoK$val
  # shuffled <- sample(length(lambda))
  #
  # lambda1 = lambda2 = lambda
  # lambda1[shuffled[(floor(length(lambda)/3)+1):length(lambda)]] = 0
  # lambda2[1:floor(length(lambda)/3*2)] = 0
  # R1 <- sqrt(lambda1)
  # R2 <- sqrt(lambda2)    # used to simulate Y
  # K1 <- V %*% diag(lambda1) %*% t(V)
  # K2 <- V %*% diag(lambda2) %*% t(V)          # form K1 and K2 for each n
  for (k in 1:N) {
    Y <- simData(h21, h22, sigma2, R1, R2, V)
    index0 <- sample(1:n, floor(n/2), replace = F)
    index1 <- (1:n)[-index0]

    # naive method
    start <- Sys.time()
    test_naive[k, i]<- varcomp::slrt_naive(y=Y, K1=Klist[1], K2=Klist[2],
                                           i1=index1, i0=index0, rho = h21, parameter.set = "h2")
    time_naive[k, i] <- difftime(Sys.time(), start, units = "secs")
    print(paste("n =", nseq[i], "; k =", k, "; Naive time used:", time_naive[k, i]))


    # test0
    start <- Sys.time()
    test_test0[k, i] <- varcomp::slrt_test0(y=Y, K1=Klist[1], K2=Klist[2],
                                            i1=index1, i0=index0, parameter.set = "h2")
    time_test0[k, i] <- difftime(Sys.time(), start, units = "secs")
    print(paste("n =", nseq[i], "; k =", k, "; Test0 time used:", time_test0[k, i]))

    # orthogonal K
    start <- Sys.time()
    eoK1 <- eigen(Klist[[1]])
    Vnew <- eoK1$vectors
    idx <- which(abs(eoK1$values) < TOL)  # find indices of 0
    basis <- Vnew[, idx]
    eoK2sub <- eigen(t(basis) %*% Klist[[2]] %*% basis)
    Vnew[, idx] <- basis %*% eoK2sub$vectors # transform the basis
    K1new <- eoK1$values
    K2new <- c(rep(0, length(K1new)-length(idx)), eoK2sub$values)
    Ynew <- crossprod(Vnew, Y)

    test_ortho[k, i] <- varcomp::slrt_ortho(y=Ynew, K1_list = list(K1new), K2_list = list(K2new),
                                            rho=h21, i1=index1, i0=index0,
                                            parameter.set = 'h2')
    time_ortho[k, i] <- difftime(Sys.time(), start, units = "secs")
    print(paste("n =", nseq[i], "; k =", k, "; Orthogonal time used:", time_ortho[k, i]))
  }
}
print("#")
#saveRDS(list(time_naive, time_test0, time_ortho), "time_M=2")



# Create the boxplot
p_2 <- data.frame(
  n = as.factor(c(rep(nseq, each = N),rep(nseq, each = N),rep(nseq, each = N))),
  Method = c(rep("Naive", times = N*length(nseq)),
             rep("Diagonalize under null", times = N*length(nseq)),
             rep("Diagonalize", times = N*length(nseq))),
  Value = (c(as.vector(time_naive), as.vector(time_test0), as.vector(time_ortho)))
) %>% ggplot(aes(x = n, y = log(Value, base=10), fill = Method)) +
  geom_boxplot() +
  labs(x = "n", y = NULL) + #title = "Side-by-side Boxplot",
  theme_bw() + theme(legend.position = "none")




## 3matrices timing

TOL <- 1e-8
N <- 30
nseq <- c(200, 500, 1000, 2000, 3000)

# data simulation with upper triangular factor of the Cholesky decomposition of K
simData <- function(h21, h22, h23, sigma2=1, R1, R2, R3, V) { # R are sqrt(eigen of K)
  n <- length(R1)
  e = rnorm(n, mean = 0, sd = sqrt(sigma2 * (1-h21-h22-h23)))
  u1 = rnorm(n, mean = 0, sd = sqrt(sigma2 * h21)) # mvrnorm(1, mu = rep(0, n), Sigma = sigma2 * h12 * K1)
  u2 = rnorm(n, mean = 0, sd = sqrt(sigma2 * h22)) # mvrnorm(1, mu = rep(0, n), Sigma = sigma2 * h22 * K2)
  u3 = rnorm(n, mean = 0, sd = sqrt(sigma2 * h23))
  y = R1 * u1 + R2 * u2 + R3 * u3 + e              #crossprod(R2, u2) + e
  return(V%*%y)
}

alpha <- 0.05
h22 <- 0
sigma2 <- 1
h21 <- 0
h23 <- 0.2
M = 3



time_naive_2 = matrix(0, nrow = N, ncol = length(nseq))
time_test0_2 = matrix(0, nrow = N, ncol = length(nseq))
time_ortho_2 = matrix(0, nrow = N, ncol = length(nseq))

test_naive_2 = matrix(0, nrow = N, ncol = length(nseq))
test_test0_2 = matrix(0, nrow = N, ncol = length(nseq))
test_ortho_2 = matrix(0, nrow = N, ncol = length(nseq))

for (i in 1:length(nseq)) {
  n <- nseq[i]
  K <- 0.5^abs(outer(1:n, 1:n, '-'))
  eoK <- eigen(K)
  V <- eoK$vectors
  lambda <- eoK$val

  n_each <- floor(n / (M + 1))
  idx_mat <- matrix(sample(n, size = M * n_each), ncol = M)
  Klist <- list()
  for(ii in 1:M){
    Klist[[ii]] <- V[, idx_mat[, ii]] %*% (lambda[idx_mat[, ii]] * t(V[, idx_mat[, ii]]))
  }
  #Rlist <- lapply(Klist, function(K) { sqrt(pmax(diag(t(V) %*% K %*% V), 0))})
  R1 <- sqrt(pmax(diag(t(V) %*% Klist[[1]] %*% V), 0))
  R2 <- sqrt(pmax(diag(t(V) %*% Klist[[2]] %*% V), 0))
  R3 <- sqrt(pmax(diag(t(V) %*% Klist[[3]] %*% V), 0))

  for (k in 1:N) {
    Y <- simData(h21, h22, h23, sigma2, R1, R2, R3, V)
    index0 <- sample(1:n, floor(n/2), replace = F)
    index1 <- (1:n)[-index0]

    # naive method
    start <- Sys.time()
    test_naive_2[k, i]<- varcomp::slrt_naive(y=Y, K1=Klist[1:2], K2=Klist[3],
                                             i1=index1, i0=index0, rho = c(h21,h22), parameter.set = "h2")
    time_naive_2[k, i] <- difftime(Sys.time(), start, units = "secs")
    print(paste("n =", nseq[i], "; k =", k, "; Naive time used:", time_naive[k, i]))


    # test0
    start <- Sys.time()
    test_test0_2[k, i] <- varcomp::slrt_test0(y=Y, K1=Klist[1:2], K2=Klist[3],
                                              i1=index1, i0=index0, parameter.set = "h2")
    time_test0_2[k, i] <- difftime(Sys.time(), start, units = "secs")
    print(paste("n =", nseq[i], "; k =", k, "; Test0 time used:", time_test0[k, i]))

    # orthogonal K
    start <- Sys.time()
    eoK1 <- eigen(Klist[[1]])
    Vnew <- eoK1$vectors
    idx1 <- which(abs(eoK1$values) < TOL)  # find indices of 0
    Vnew_good <- Vnew[, -idx1]
    Vnew_null <- Vnew[, idx1]

    K2_sub <- t(Vnew_null) %*% Klist[[2]] %*% Vnew_null # rotate K2 with null space of K1
    eoK2_sub <- eigen(K2_sub)
    idx2 <- which(abs(eoK2_sub$values) < TOL)
    V_K2 <- Vnew_null %*% eoK2_sub$vectors
    Vnew_good <- cbind(Vnew_good, V_K2[, -idx2])
    Vnew_null <- V_K2[, idx2]

    K3_sub <- t(Vnew_null) %*% Klist[[3]] %*% Vnew_null # rotate K3 with null space of K1 and K2
    eoK3_sub <- eigen(K3_sub)
    Vnew_null <- Vnew_null %*% eoK3_sub$vectors
    Vnew <- cbind(Vnew_good, Vnew_null)

    K1new <- eoK1$values
    K2new <- c(rep(0, length(eoK1$values)-length(idx1)), eoK2_sub$values)
    K3new <- c(rep(0, length(eoK1$values)-length(idx1)+length(eoK2_sub$values)-length(idx2)), eoK3_sub$values)
    Ynew <- crossprod(Vnew, Y)

    test_ortho_2[k, i] <- varcomp::slrt_ortho(y=Ynew, K1_list = list(K1new, K2new), K2_list = list(K3new),
                                              rho=c(h21,h22), i1=index1, i0=index0,
                                              parameter.set = 'h2')
    time_ortho_2[k, i] <- difftime(Sys.time(), start, units = "secs")
    print(paste("n =", nseq[i], "; k =", k, "; Orthogonal time used:", time_ortho[k, i]))
  }
}
print("#")
#saveRDS(list(time_naive_2, time_test0_2, time_ortho_2), "time_M=3")


# Create the boxplot
p_3 <- data.frame(
  n = as.factor(c(rep(nseq, each = N),rep(nseq, each = N),rep(nseq, each = N))),
  Method = c(rep("Naive", times = N*length(nseq)),
             rep("Diagonalize under null", times = N*length(nseq)),
             rep("Diagonalize", times = N*length(nseq))),
  Value = (c(as.vector(time_naive_2), as.vector(time_test0_2), as.vector(time_ortho_2)))
) %>% ggplot(df5, aes(x = n, y = log(Value, base=10), fill = Method)) +
  geom_boxplot() +
  labs(x = "n", y = NULL) + #title = "Side-by-side Boxplot",
  theme_bw() #+ theme(legend.position = "none")
ggsave("timeplot_chp3.pdf", plot = p_2+p_3, width = 12, height = 4, dpi = 300)

