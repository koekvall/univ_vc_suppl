# data example: hicksturner, n1=10, n2=3, replications =2
set.seed(11)
library(RLRsim)
library(lme4)
library(varTestnlme)
library(patchwork)
library(ggplot2)
library(dplyr)
library(merDeriv)

# SET WORKING DIRECTORY
#setwd("~/GitHub/supp_univ/Data example/")

###################### check lmer result #####################
hicksturner <- read.table("hicksturner.txt", row.names = 1)
colnames(hicksturner) <- c("part", "oper", "mohms")

# to factor
hicksturner <- within(hicksturner, {
  oper <- as.factor(oper)
  part <- as.factor(part)})

# center Y
hicksturner$stdmohms <- (hicksturner$mohms - mean(hicksturner$mohms))

# fit lmer model
mod0 <- lmer(stdmohms ~ 0 + (1|part) + (1|oper), data = hicksturner, REML = FALSE)

# Get SEs from observed information
obs_inf <- merDeriv::vcov.lmerMod(mod0, full = TRUE, information = "observed",
                                  ranpar = "var")
sig2_se <- sqrt(diag(obs_inf))

# Compute CIs
ci.profile <- confint(mod0, method = "profile")
ci.boot <- confint(mod0, method = "boot")

sig2_hat <- as.data.frame(VarCorr(mod0))$vcov
ci.wald <- sqrt(pmax(cbind(sig2_hat - 1.96 * sig2_se, sig2_hat + 1.96 * sig2_se), 0))
colnames(ci.wald) <- colnames(ci.boot)
rownames(ci.wald) <- rownames(ci.boot)

# SE for h2 (delta method)

h_fun <- function(v) {
  total <- sum(v)
  c(h2_part = v[1] / total,
    h2_oper = v[2] / total,
    total = total)
}

h2_hat <- h_fun(sig2_hat)
J <- numDeriv::jacobian(h_fun, sig2_hat)
h2_se <- sqrt(diag(J %*% obs_inf %*% t(J)))

# Compare to simulation-based testing of boundary points; only works if
# including intercept and using REML
mod_full <-  lmer(stdmohms ~ 1 + (1|part) + (1|oper), data = hicksturner,
                  REML = TRUE)
mod_oper <- lmer(stdmohms ~ 1 + (1|oper), data = hicksturner, REML = TRUE)
mod_part <- lmer(stdmohms ~ 1 + (1|part), data = hicksturner, REML = TRUE)

# Test whether "part" random effect variance is zero
simtest_part <- RLRsim::exactRLRT(m = mod_part, mA = mod_full, m0 = mod_oper)

# Test whether "operator" random effect variance is zero
simtest_oper <- RLRsim::exactRLRT(m = mod_oper, mA = mod_full, m0 = mod_part)


######################### model setup #######################
K1 <- kronecker(diag(10), matrix(1, 6, 6))
K2 <- kronecker(matrix(1, 20, 20), diag(3))

# The eigenvectors of Sigma are the same for every parameter value, and hence
# the same as those of K1 + K2
O <- eigen(K1 + K2, symmetric = TRUE)$vectors
lam1 <- diag(crossprod(O, K1) %*% O)
lam2 <- diag(crossprod(O, K2) %*% O)

y_O <- crossprod(O, hicksturner$stdmohms)
n <- length(y)


############################# test statistic with split LRT ####################

# Values to test
sigma21_seq <- seq(0, 50, length.out = 500)
sigma22_seq <- c(0, 2^seq(-3, ceiling(log2(10000)), length.out = 499))

# Randomized threshold
ru <- runif(1)

# Random partition of y
kfold <- 4
group <- split(sample(n), 1:kfold)

# Run tests
test_seq <- rep(0, length(sigma21_seq))
test_seq2 <- rep(0, length(sigma22_seq))
for (kk in 1:kfold) {
  index1 <- group[[kk]]
  index0 <- unlist(group[-kk])
  for (ii in 1:length(sigma21_seq)) {
    test_seq[ii] <- test_seq[ii] + varcomp::slrt_ortho(y = y_O,
                                                        # Test sigma_1^2
                                                        K1_list = list(lam1),
                                                        K2_list = list(lam2),
                                                        i1 = index1,
                                                        i0 = index0,
                                                        rho = sigma21_seq[ii],
                                                        parameter.set = "sigma2")

    test_seq2[ii] = test_seq2[ii] + varcomp::slrt_ortho(y = y_O,
                                                        # Test sigma_2^2
                                                        K1_list = list(lam2),
                                                        K2_list = list(lam1),
                                                        i1 = index1,
                                                        i0 = index0,
                                                        rho = sigma22_seq[ii],
                                                        parameter.set = "sigma2")
    if(ii %% 100 == 0) {print(paste(i, "th test finished."))}
  }
  print(paste(kk, "th fold finished."))
}
test_seq = test_seq / kfold
test_seq2 = test_seq2 / kfold

# p-values for testing zero
min(c(ru / test_seq[1], 1))
min(c(ru / test_seq2[1],1))

# CI using randomized SLRT
lower_end1 <- sqrt(sigma21_seq[min(which(test_seq < ru / alpha))])
upper_end1 <- sqrt(sigma21_seq[max(which(test_seq < ru / alpha))])
lower_end2 <- sqrt(sigma22_seq[min(which(test_seq2 < ru / alpha))])
upper_end2 <- sqrt(sigma22_seq[max(which(test_seq2 < ru / alpha))])
ci.slrt <- cbind(c(lower_end1, lower_end2), c(upper_end1, upper_end2))


p1 <- data.frame(tested_value = sqrt(sigma21_seq), test_statistics = test_seq) %>%
  ggplot(aes(x = tested_value, y = test_statistics)) +
  geom_line() +  # Scatter plot
  geom_hline(aes(yintercept = 1/alpha), linetype = "dashed", size = 0.5) +
  geom_hline(aes(yintercept = ru/alpha), linetype = "dotdash", size = 0.5, alpha = 0.7) +
  theme_bw() + ylim(0, 1/alpha+10) + #xlim(0, 1) +
  labs(x = expression(sigma[1]), y = "test statistics")

p2 <- data.frame(tested_value = sqrt(sigma22_seq), test_statistics = test_seq2) %>%
  ggplot(aes(x = tested_value, y = test_statistics)) +
  geom_line() +  # Scatter plot
  geom_hline(aes(yintercept = 1/alpha, linetype = "1/alpha"), size = 0.5) +
  geom_hline(aes(yintercept = ru/alpha, linetype = "ru/alpha"), size = 0.5, alpha = 0.7) +
  scale_linetype_manual(name = "Critical values", values = c("1/alpha" = "dashed", "ru/alpha" = "dotdash"),
                        labels = c("1/alpha" = expression(1 / alpha),"ru/alpha" = expression(Unif(0,1) / alpha))) +
  theme_bw() + ylim(0, 1/alpha+10) + #xlim(0, 18) +
  labs(x = expression(sigma[2]), y = NULL)

p1+p2
ggsave("CI_hicksturner.pdf", plot = p1+p2, width = 12, height = 4, dpi = 300)

###############################################################################
# Simulations for distribution of CI widths
###############################################################################
N <- 1000
ci1 <- matrix(0, nrow = N, ncol = 2)
ci2 <- matrix(0, nrow = N, ncol = 2)
for (iN in 1:N) {
  u <- runif(1)
  # Width for sigma_1
  lower_end <- sqrt(sigma21_seq[min(which(test_seq < u / alpha))])
  upper_end <- sqrt(sigma21_seq[max(which(test_seq < u / alpha))])
  if(is.na(lower_end)){
    len1[iN] <- upper_end
  } else{
    len1[iN] <- upper_end - lower_end
  }
  if(is.na(upper_end)){
    len1[iN] <- 0
  }
  
  # Width for sigma_2
  lower_end <- sqrt(sigma22_seq[min(which(test_seq2 < u / alpha))])
  upper_end <- sqrt(sigma22_seq[max(which(test_seq2 < u/alpha))])
  if(is.na(lower_end)){
    len2[iN] <- upper_end
  } else{
    len2[iN] <- upper_end - lower_end
  }
  if(is.na(upper_end)){
    len2[iN] <- 0
  }
}

# Compare average lengths
mean((len1)) /  (sqrt(sigma21_seq[max(which(test_seq < 1/alpha))]) - sqrt(sigma21_seq[min(which(test_seq < 1/alpha))]))
mean((len2)) /  (sqrt(sigma22_seq[max(which(test_seq2 < 1/alpha))]) - sqrt(sigma22_seq[min(which(test_seq2 < 1/alpha))]))

p01 <- data.frame(len = len1) %>%
  ggplot(aes(x = len)) +
  geom_density(alpha = 0.5, fill = "grey") + theme_bw() + labs(x = expression("Width of CI for " * sigma[1]), y = "Density")

p02 <- data.frame(len = len2) %>%
  ggplot(aes(x = len)) +
  geom_density(alpha = 0.5, fill = "grey") + theme_bw() + labs(x = expression("Width of CI for " * sigma[2]), y = NULL)

p01 + p02

ggsave("hicksturner_CIlength.pdf", plot = p01+p02, width = 12, height = 4, dpi = 300)

Z <- getME(mod0, "Z")
G <- bdiag(sig2_hat[1] * diag(10), sig2_hat[2] * diag(3))
W <- Z %*% G %*%t(Z) + sig2_hat[3] * diag(60)
uhat <- G %*% t(Z ) %*% solve(W, hicksturner$stdmohms)
yhat <- Z %*% uhat
ehat <- hicksturner$stdmohms - yhat


resplot <- data.frame(Fitted = as.vector(yhat), Residuals = as.vector(ehat)) %>%
  ggplot(aes(x = Fitted, y = Residuals)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(title = "", x = "Predicted conditional mean", y = "Predicted error") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))
qqplot <- data.frame(resid = as.vector(ehat)) %>%
  ggplot(aes(sample = resid)) +
  stat_qq(alpha = 0.5) +
  stat_qq_line() +
  labs(title = "", x = "Standard normal quantiles", y = "Empirical quantiles of predicted errors") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))

qqplot + resplot
ggsave("diagnostic_hicksturner.pdf", plot = qqplot + resplot, width = 12, height = 4, dpi = 300)
