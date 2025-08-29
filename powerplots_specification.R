### simulation method
set.seed(145)
library(Matrix)
library(MASS)
library(varcomp) # devtools::install_github("yqzhang5972/varcomp")
library(dplyr)
library(ggplot2)
set.seed(145)
num_sims <- 1e2 # Publication figures created with 1e4
num_params <- 5 # Publication figures created with 10

# data simulation with upper triangular factor of the Cholesky decomposition of K
simData.pow <- function(h21, h22, sigma2=1, R1, R2, V2) { # R are sqrt(eigen of K), V is eigenvalue of K2
  n <- length(R1)
  e = rnorm(n, mean = 0, sd = sqrt(sigma2 * (1-h21-h22)))
  u1 = rnorm(n, mean = 0, sd = sqrt(sigma2 * h21)) # mvrnorm(1, mu = rep(0, n), Sigma = sigma2 * h12 * K1)
  u2 = rnorm(n, mean = 0, sd = sqrt(sigma2 * h22)) # mvrnorm(1, mu = rep(0, n), Sigma = sigma2 * h22 * K2)
  y = R1 * u1 + V2 %*% (R2 * u2) + e              #crossprod(R2, u2) + e
  return(y)
}


### formating K

TOL <- 1e-8
alpha <- 0.05
n <- 300
N <- num_sims

h22 <- 0.2
sigma2 <- 1
h21true <- 0 #.5  # test for different h21, each with k to see the sizes
h21 <- seq(0.001, 0.9999-h22, length.out = num_params)
const <- 100
L1true = L2true = numeric(n)
L1true[1:floor(n/3)] = seq(10, 5, length.out = floor(n/3))
L2true[1:floor(n/2.5)] = seq(10, 5, length.out = floor(n/2.5))
L1true[(floor(n/3)+1):n] = seq(5, 0, length.out = n-floor(n/3)) / const
L2true[(floor(n/2.5)+1):n] = seq(5, 0, length.out = n-floor(n/2.5)) / const
R1 <- sqrt(L1true)
R2 <- sqrt(L2true)    # used to simulate Y

K1true <- diag(L1true)
V <- as.matrix(bdiag(diag(floor(n/2.5)),
                     svd(matrix(rnorm((n-floor(n/2.5))*(n-floor(n/2.5))), nrow = n-floor(n/2.5)))$u))
K2true <- V %*% diag(L2true) %*% t(V)

L1simp = L1true
L2simp = L2true
L1simp[(floor(n/3)+1):n] = 0
L2simp[(floor(n/2.5)+1):n] = 0
# all.equal(diag(300), crossprod(V)) V is orthogonal
# all.equal(t(V)%*%diag(L1simp)%*%V, diag(L1simp)) V is one set of eigenvectors for L1simp
# all.equal(t(V)%*%diag(L2simp)%*%V, diag(L2simp)) V is one set of eigenvectors for L2simp



##################### ordinary power plot ####################

testfun <- function(k) {
  test <- rep(0, length(h21))
  for (i in 1:length(h21)) {
    index0 <- sample(1:n, floor(n/2), replace = F)
    index1 <- (1:n)[-index0]
    y <- simData.pow(h21true, h22, sigma2, R1, R2, V2=V) # let true value to be 0
    # ynew <- crossprod(V, y)

    # test[i]<- varcomp::slrt_naive(y=ynew, K1_list = list(K1true), K2_list=list(K2true),
    #                                   i1=index1, i0=index0, rho = h21[i])
    test[i]<- varcomp::slrt_naive(y=y, K1_list = list(K1true), K2_list=list(K2true),
                         i1=index1, i0=index0, rho = h21[i])
  }
  if(k %% 10 == 0) {system(paste("echo 'now processing:",k,"'", ", system time:", Sys.time()))}
  return(test)
}
start <- Sys.time()
ord_test <- parallel::mclapply(1:N, testfun, mc.cores = getOption("mc.cores", 6L))
difftime(start, Sys.time(), units = "secs")
ord_test <- matrix(unlist(ord_test), nrow = N, byrow = T)
#saveRDS(ord_test, "power_ordinary")

### find power
p_ord <- ord_test# readRDS(file = "power_ordinary")
for (k in 1:nrow(p_ord)) {
  for (i in 1:length(h21)) {
    p_ord[k,i] <- (p_ord[k,i] > (runif(1)/alpha)) # runif(1)
  }
}
(cov_ord <- apply(p_ord, 2, mean))
er_ord <- sqrt(cov_ord*(1-cov_ord)/nrow(p_ord)) * 1.96

### making plots


plot_ord <- data.frame(h21seq = rep(1:length(h21), 1),
                       ub = c(cov_ord+er_ord),
                       lb = c(cov_ord-er_ord),
                       mp = c(cov_ord)) %>%                   # n = 1000, h22 = 0.2
  ggplot(aes(x = h21seq, y = mp)) + #ggtitle(expression(h["1*"]^2 == 0)) +
  xlab(expression(paste("Tested ", h[1]^2))) + ylab(expression(paste("Power with true ", h["1*"]^2, " = 0"))) + coord_cartesian(ylim = c(0, 1)) +
  geom_ribbon(aes(x = h21seq, ymin = lb, ymax =ub), alpha= 0.3, colour = NA) +
  theme_bw() + geom_hline(yintercept = 1) + geom_hline(yintercept = 0.05, linetype = "dashed") +
  scale_x_continuous(breaks = 1:length(h21), labels = round(h21,3)) +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = c(0.15, 0.15),legend.key.size = unit(0.2, 'cm'))
plot_ord





########################## ignore small eigenvalues ###################

testfun <- function(k) {
  test <- rep(0, length(h21))
  for (i in 1:length(h21)) {
    index0 <- sample(1:n, floor(n/2), replace = F)
    index1 <- (1:n)[-index0]
    y <- simData.pow(h21true, h22, sigma2, R1, R2, V2=V) # let true value to be 0
    ynew <- crossprod(V, y)
    #sink(file = nullfile())
    test[i]<- slrt_ortho(y=ynew, K1_list = list(L1simp), K2_list=list(L2simp),
                         i1=index1, i0=index0, rho = h21[i])
    #sink()
  }
  if(k %% 100 == 0) {system(paste("echo 'now processing:",k,"'", ", system time:", Sys.time()))}
  return(test)
}
start <- Sys.time()
simp_test <- parallel::mclapply(1:N, testfun, mc.cores = getOption("mc.cores", 6L))
difftime(start, Sys.time(), units = "secs")

simp_test <- matrix(unlist(simp_test), nrow = N, byrow = T)
# saveRDS(simp_test, "power_simplified")

### find power
p_simp <- simp_test# readRDS(file = "power_simplified")
for (k in 1:nrow(p_simp)) {
  for (i in 1:length(h21)) {
    p_simp[k,i] <- (p_simp[k,i] > (runif(1)/alpha)) # runif(1)
  }
}
(cov_simp <- apply(p_simp, 2, mean))
er_simp <- sqrt(cov_simp*(1-cov_simp)/nrow(p_simp)) * 1.96

### making plots
plot_simp <- data.frame(h21seq = rep(1:length(h21), 1),
                        ub = c(cov_simp+er_simp),
                        lb = c(cov_simp-er_simp),
                        mp = c(cov_simp)) %>%                   # n = 1000, h22 = 0.2
  ggplot(aes(x = h21seq, y = mp)) + # ggtitle(expression(h["1*"]^2 == 0)) +
  xlab(expression(paste("Tested ", h[1]^2))) + coord_cartesian(ylim = c(0, 1)) + # ylab("Power")  +
  geom_ribbon(aes(x = h21seq, ymin = lb, ymax =ub), alpha= 0.3, colour = NA) +
  theme_bw() + geom_hline(yintercept = 1) + geom_hline(yintercept = 0.05, linetype = "dashed") +
  scale_x_continuous(breaks = 1:length(h21), labels = round(h21,3)) +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = c(0.15, 0.15),legend.key.size = unit(0.2, 'cm'))
plot_simp



######################## ignore constrains #####################

testfun <- function(k) {
  test <- rep(0, length(h21))
  for (i in 1:length(h21)) {
    index0 <- sample(1:n, floor(n/2), replace = F)
    index1 <- (1:n)[-index0]
    y <- simData.pow(h21true, h22, sigma2, R1, R2, V2=V) # let true value to be 0
    # ynew <- crossprod(V, y)
    test[i]<- varcomp::slrt_naive_uncon(y=y, K1_list = list(K1true), K2_list=list(K2true),
                               i1=index1, i0=index0, rho = h21[i])
  }
  if(k %% 100 == 0) {system(paste("echo 'now processing:",k,"'", ", system time:", Sys.time()))}
  return(test)
}
start <- Sys.time()
uncon_test <- parallel::mclapply(1:N, testfun, mc.cores = getOption("mc.cores", 6L))
difftime(start, Sys.time(), units = "secs")

uncon_test <- matrix(unlist(uncon_test), nrow = N, byrow = T)
# saveRDS(uncon_test, "power_unconstrained")

### find power
p_uncon <- uncon_test # readRDS(file = "power_unconstrained")
for (k in 1:nrow(p_uncon)) {
  for (i in 1:length(h21)) {
    p_uncon[k,i] <- (p_uncon[k,i] > (runif(1)/alpha)) # runif(1)
  }
}
(cov_uncon <- apply(p_uncon, 2, mean))
er_uncon <- sqrt(cov_uncon*(1-cov_uncon)/nrow(p_uncon)) * 1.96

### making plots
plot_uncon <- data.frame(h21seq = rep(1:length(h21), 1),
                         ub = c(cov_uncon+er_uncon),
                         lb = c(cov_uncon-er_uncon),
                         mp = c(cov_uncon)) %>%                   # n = 1000, h22 = 0.2
  ggplot(aes(x = h21seq, y = mp)) + # ggtitle(expression(h["1*"]^2 == 0)) +
  xlab(expression(paste("Tested ", h[1]^2))) + coord_cartesian(ylim = c(0, 1)) + #ylab("Power") +
  geom_ribbon(aes(x = h21seq, ymin = lb, ymax =ub), alpha= 0.3, colour = NA) +
  theme_bw() + geom_hline(yintercept = 1) + geom_hline(yintercept = 0.05, linetype = "dashed") +
  scale_x_continuous(breaks = 1:length(h21), labels = round(h21,3)) +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = c(0.15, 0.15),legend.key.size = unit(0.2, 'cm'))
plot_uncon


#plot_ord+plot_simp+plot_uncon
#ggsave("power.pdf", plot = plot_ord+plot_simp+plot_uncon, width = 12, height = 4, dpi = 300)
