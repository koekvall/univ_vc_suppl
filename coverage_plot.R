library(tidyverse)
library(patchwork)
library(ggplot2)
library(varcomp) # devtools::install_github("yqzhang5972/varcomp")
set.seed(145)
num_sims <- 1e3 # Publication figures created with 1e4
num_params <- 20 # Publication figures created with 40

# ---------------------------------------
# Make K1, K2
# ---------------------------------------
n <- 300
K1 <- matrix(0, nrow=n, ncol=n)
for(jj in 1:n){
  for(kk in 1:n){
    K1[jj,kk] <- 0.95^abs(jj - kk)
  }
}
K1 <- eigen(K1)$values

# K1 <- seq(1, 5, length.out = 300)

K2 <- 0.5^abs(outer(1:n, 1:n, '-'))
K2 <- eigen(K2)$values

X <- matrix(rep(1,n), ncol = 1) # matrix(rnorm(n), ncol = 1) #






# -------------------------------------------------
# log-likelihood function for h21, h22, s2p, par = h22
# -------------------------------------------------

loglik_score1 <- function(par, rhoh21, X, y, K1, K2, return.s2phat= FALSE) {
  n <- nrow(X)
  p <- ncol(X)
  Sigma0 <- rhoh21*K1 + par*K2 + rep(1-rhoh21-par,n)
  Sigma0_inv <- 1 / Sigma0 # Sigma_inv = Sigma0_inv / s2p
  XS0iX <- crossprod(X, Sigma0_inv * X) # XSX = XS0iX / s2p
  betahat <- solve(XS0iX, crossprod(X, Sigma0_inv * y))

  # yQ0y <- sum((y-X%*%betahat) * Sigma0_inv * (y-X%*%betahat)) # y'Qy = yQ0y / s2p
  s2phat <- sum((y-X%*%betahat) * Sigma0_inv * (y-X%*%betahat)) / (n-p)
  # print(paste("betahat=", betahat,"s2phat=", s2phat))
  if(return.s2phat) {
    return(s2phat)
  }
  l <- -(sum(log(Sigma0)) + (n-p)*log(s2phat) + determinant(XS0iX)$modulus[1] + n-p) / 2   #yQ0y / s2phat)
  return(l)
}


# U = -1/2tr(QV1) + 1/2y^TQ V1 Qy
# Q = Sigma^{-1} - Sigma^{-1}X(X'Sigma^{-1}X)^{-1}X'Sigma^{-1}
U_h21 <- function(h21, h22, s2p, X, y, K1, K2) {
  V1 <- s2p * (K1 - rep(1, n))
  Sigma <- (h21*K1 + h22*K2 + rep(1-h21-h22,n)) * s2p
  Sigma_inv <- 1 / Sigma
  XtSi <- t(Sigma_inv * X)            # = t(X) %*% diag(Sigma_inv) = X'Sigma^{-1}
  Q <- diag(Sigma_inv) - crossprod(XtSi, solve(XtSi%*%X, XtSi))    # matrix
  Qy2 <- crossprod(Q, y) ^ 2          #
  u1 <- (-sum(diag(Q)*V1) + sum(Qy2 * V1)) / 2
  return(u1)
}

# I_{ij} <- 1/2tr(QViQVj)
i11_inv <- function(h21, h22, s2p, X, y, K1, K2) {
  n <- nrow(X)
  p <- ncol(X)
  V1 <- s2p * (K1 - rep(1, n))
  V2 <- s2p * (K2 - rep(1, n))
  Sigma <- (h21*K1 + h22*K2 + rep(1-h21-h22,n)) * s2p
  Sigma_inv <- 1 / Sigma
  XtSi <- t(Sigma_inv * X)
  Q <- diag(Sigma_inv) - crossprod(XtSi, solve(XtSi%*%X, XtSi))


  i11 <- (sum((t(Q)*Q) * outer(V1,V1))) / 2 # sum((t(Q)*Q) * outer(V1,V2))
  i12 <- (sum((t(Q)*Q) * outer(V1,V2))) / 2
  i22 <- (sum((t(Q)*Q) * outer(V2,V2))) / 2
  i13 <- sum(diag(Q)*V1) / s2p / 2
  i23 <- sum(diag(Q)*V2) / s2p / 2
  i33 <- (n-p) / s2p^2 / 2

  detI <- i11 * (i22*i33 - i23^2) - i12 * (i12*i33 - i23*i13) + i13 * (i12*i23 - i22*i13)

  return((i22*i33 - i23^2) / detI)
}


#################### Introduction: h22 varies for score based test ####################
h22true <- seq(0,1, length.out = num_params)
h21true <- 0
s2ptrue <- 1
N <- num_sims
#N0 <- 10000

teststat <- rep(0, N)
coverage <- rep(0, length(h22true))


for (c in 1:length(h22true)) {
  for (k in 1:N) {
    # simulate y
    y <- rnorm(n, sd = sqrt((1-h21true-h22true[c])*s2ptrue)) + sqrt(K1) * rnorm(n, sd = sqrt(h21true*s2ptrue)) + sqrt(K2) * rnorm(n ,sd = sqrt(h22true[c]*s2ptrue)) #  true h22=beta=0 //
    # y <- crossprod(V, y)

    # method 1
    parhat <- optim(1e-8, loglik_score1, control = list(fnscale=-1), rhoh21 = h21true,
                    X=X, y=y, K1=K1, K2=K2, method = "L-BFGS-B", lower = c(0), upper = c(1-1e-8-h21true))
    h22hat <- parhat$par
    s2phat <- loglik_score1(par = h22hat, rhoh21 = h21true, X=X, y=y, K1=K1, K2=K2, return.s2phat= TRUE)

    u1 <- U_h21(h21=h21true, h22 = h22hat, s2p = s2phat,
                X=X, y=y, K1=K1, K2=K2)
    i11i <- i11_inv(h21=h21true, h22 = h22hat, s2p = s2phat,
                    X=X, y=y, K1=K1, K2=K2)

    teststat[k] <- u1^2 * i11i
    if(k%%1000 == 0) print(paste(c,"th coverage:", k, "th finished."))
  }
  coverage[c] <- sum(teststat < qchisq(0.95, 1)) / N
  print(paste("*****", c, "th coverage finished.******"))
}
coverage_score <- coverage
# save(coverage_score, file = "cov_moti_chp3_2.RData")

p <- data.frame(tested = h22true, value = coverage_score) %>%
  mutate(
    # Calculate the standard error of the proportion
    se = 2*sqrt(value * (1 - value) / N),
    # Calculate the lower bound for the ribbon
    lower_bound = value - se,
    # Calculate the upper bound for the ribbon
    upper_bound = value + se
  ) %>%
  ggplot(aes(x = tested, y = value)) +
  geom_line() +
  geom_ribbon(aes(ymin = lower_bound, ymax = upper_bound), alpha = 0.4) +
  # geom_point(alpha = 0.7) + geom_line(alpha = 0.7) +# Add scatter points, set transparency and color
  labs(x=expression(paste(h[2]^2)), y="Coverage") + ylim(c(0.85,1))+
  geom_hline(yintercept = 0.95, linetype = "dashed", size = 0.7) + theme_bw()
# ggsave("motivationforchp3_h22various.pdf", plot = p, width = 12, height = 4, dpi = 300)





###################### Simulation: h21 varies for both score test and split LRT ###################
TOL <- 1e-8
h22true <- 0
h21true <- seq(0,0.9999, length.out = num_params)
s2ptrue <- 1
N <- num_sims

# score test
teststat <- rep(0, N)
coverage <- rep(0, length(h21true))


for (c in 1:length(h21true)) {
  for (k in 1:N) {
    # simulate y
    y <- rnorm(n, sd = sqrt((1-h21true[c]-h22true)*s2ptrue)) + sqrt(K1) * rnorm(n, sd = sqrt(h21true[c]*s2ptrue)) + sqrt(K2) * rnorm(n ,sd = sqrt(h22true*s2ptrue)) #  true h22=beta=0 //
    # y <- crossprod(V, y)

    # method 1
    parhat <- optim(0.1, loglik_score1, control = list(fnscale=-1), rhoh21 = h21true[c],
                    X=X, y=y, K1=K1, K2=K2, method = "L-BFGS-B", lower = c(0), upper = c(1-1e-8-h21true[c]))
    h22hat <- parhat$par
    s2phat <- loglik_score1(par = h22hat, rhoh21 = h21true[c], X=X, y=y, K1=K1, K2=K2, return.s2phat= TRUE)

    u1 <- U_h21(h21=h21true[c], h22 = h22hat, s2p = s2phat,
                X=X, y=y, K1=K1, K2=K2)
    i11i <- i11_inv(h21=h21true[c], h22 = h22hat, s2p = s2phat,
                    X=X, y=y, K1=K1, K2=K2)

    teststat[k] <- u1^2 * i11i
    if(k%%1000 == 0) print(paste(c,"th coverage:", k, "th finished."))
  }
  coverage[c] <- sum(teststat < qchisq(0.95, 1)) / N
  print(paste("*****", c, "th coverage finished.******"))
}
coverage_score <- coverage
# save(coverage_score, file = "cov_moti_chp3.RData")

p1<-data.frame(tested = h21true, value = coverage_score) %>%
  mutate(
    # Calculate the standard error of the proportion
    se = 2*sqrt(value * (1 - value) / N),
    # Calculate the lower bound for the ribbon
    lower_bound = value - se,
    # Calculate the upper bound for the ribbon
    upper_bound = value + se
  ) %>%
  ggplot(aes(x = tested, y = value)) +
  geom_line() +
  geom_ribbon(aes(ymin = lower_bound, ymax = upper_bound), alpha = 0.4) +
  # geom_point(alpha = 0.7) + geom_line(alpha = 0.7) +# Add scatter points, set transparency and color
  labs(x=expression(paste(h[1]^2)), y="Coverage") + #ylim(c(0.85,1))+
  geom_hline(yintercept = 0.95, linetype = "dashed", size = 0.7) + theme_bw()
# ggsave("motivationforchp3_h22various.pdf", plot = p, width = 12, height = 4, dpi = 300)
# save(coverage, file = "moti_10000.RData")

p1

# split LRT
teststat <- rep(0, N)
coverage <- rep(0, length(h21true))


for (c in 1:length(h21true)) {
  for (k in 1:N) {
    index0 <- sample(1:n, floor(n/2), replace = F)
    index1 <- (1:n)[-index0]
    y <- rnorm(n, sd = sqrt((1-h21true[c]-h22true)*s2ptrue)) +
      sqrt(K1) * rnorm(n, sd = sqrt(h21true[c]*s2ptrue)) +
      sqrt(K2) * rnorm(n ,sd = sqrt(h22true*s2ptrue))

    teststat[k] <- varcomp::slrt_ortho(y=y, rho=h21true[c], K1_list = list(K1), K2_list = list(K2), i1 = index1, i0=index0)
    if(k%%200 == 0) print(paste(c,"th coverage:", k, "th finished."))
  }
  coverage[c] <- sum(teststat < runif(N) / 0.05) / N
  print(paste("*****", c, "th coverage finished.******"))
}
coverage_slrt <- coverage
# save(coverage_slrt, file = "cov_moti_chp3_2_slrt.RData")

p2<-data.frame(tested = h21true, value = coverage_slrt) %>%
  mutate(
    # Calculate the standard error of the proportion
    se = 2*sqrt(value * (1 - value) / N),
    # Calculate the lower bound for the ribbon
    lower_bound = value - se,
    # Calculate the upper bound for the ribbon
    upper_bound = value + se
  ) %>%
  ggplot(aes(x = tested, y = value)) +
  geom_line() +
  geom_ribbon(aes(ymin = lower_bound, ymax = upper_bound), alpha = 0.4) +
  # geom_point(alpha = 0.5) + geom_line(alpha = 0.5) +# Add scatter points, set transparency and color
  labs(x=expression(paste(h[1]^2))) + ylim(c(0.8,1)) + # , y="Coverage of split-LRT"
  geom_hline(yintercept = 0.95, linetype = "dashed", size = 0.7) + theme_bw()
p1+p2
#ggsave("coverage_chp3.pdf", plot = p1+p2, width = 12, height = 4, dpi = 300)

