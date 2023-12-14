#############
#############
###
### STATS 230: Statistical Computing
### Fall 2023 Final Project
### Author: Yongjun Lee
### Description: Implementation of the AIMH algorithm
###              proposed by Giordani and Kohn (2010)
###
#############
#############

library(latex2exp)
library(kableExtra)

# Experiment 1
##############

# TARGET - AIMH
set.seed(123)
x <- matrix(0, ncol=2)
mu <- rbind(c(-3,-2), c(3,-2), c(0,3))
mu <- rbind(c(-3,-1.5), c(1.5,1.5), c(-3,2))
sigma <- list(diag(2), 1.2*diag(2), 0.7*diag(2))
mix_weight <- c(0.3,0.4,0.3)
target <- list(mu=mu, sigma=sigma, mix_weight=mix_weight)
N = 2

#rslt <- readRDS("experiment1")
rslt <- AIMH(x, k=16, om1=0.05, om2=0.15, L=10, a_L=0.1,
             M=20, a_M=0.02, target, N=N)
# true posterior distribution
x <- c(0,0)
for (i in 1:500) {
  x <- rbind(x, sim_mix_normal(rslt[['lamb_post']]))
}
xlim <- c(min(rslt[['z']][,1]), max(rslt[['z']][,1]))
ylim <- c(min(rslt[['z']][,2]), max(rslt[['z']][,2]))
par(mfrow=c(1,2))
plot(x, main="Target Density", xlab="", ylab="",xlim=xlim, ylim=ylim, pch=4)
plot(rslt[['z']], main="Simulated Density", xlab="", ylab="",xlim=xlim, ylim=ylim, pch=4)


# Q-ITER
# n 0
lamb_init <- rslt[['q_cache']][[1]]
q_initial <-c()
n <- 2
for (i in 1:n) {
  q_initial <- rbind(q_initial, sim_mix_normal(lamb_init))
}
# n 200
lamb_init <- rslt[['q_cache']][[2]][[1]]
lamb_200 <- rslt[['q_cache']][[2]][[2]]
q_200 <- c()
for (i in 1:n) {
  q_200 <- rbind(q_200, q(lamb_init, lamb_200, k=16, om1=0.05, om2=0.15)[[1]])
}
# n 500
lamb_init <- rslt[['q_cache']][[3]][[1]]
lamb_500 <- rslt[['q_cache']][[3]][[2]]
q_500 <- c()
for (i in 1:n) {
  q_500 <- rbind(q_200, q(lamb_init, lamb_500, k=16, om1=0.05, om2=0.15)[[1]])
}
# n last
lamb_init <- rslt[['q_cache']][[4]][[1]]
lamb_last <- rslt[['q_cache']][[4]][[2]]
q_last <- c()
for (i in 1:n) {
  q_last <- rbind(q_last, q(lamb_init, lamb_last, k=16, om1=0.05, om2=0.15)[[1]])
}
par(mfrow=c(2,2))
plot(q_initial, pch=4, main = TeX(r"($q_0(z;\lambda)$)", bold=TRUE), xlab="", ylab="", ylim=c(-6,6), xlim=c(-6,6))
points(lamb_init[['mu']],col='green', pch=4, lwd=4)

plot(q_200, pch=4, main = TeX(r"($q_{200}(z;\lambda)$)", bold=TRUE), xlab="", ylab="", ylim=c(-6,6), xlim=c(-6,6))
points(lamb_init[['mu']],col='green', pch=4, lwd=4)
points(lamb_200[['mu']],col='red', pch=4, lwd=4)

plot(q_500, pch=4, main = TeX(r"($q_{500}(z;\lambda)$)", bold=TRUE), xlab="", ylab="", ylim=c(-6,6), xlim=c(-6,6))
points(lamb_init[['mu']],col='green', pch=4, lwd=4)
points(lamb_500[['mu']],col='red', pch=4, lwd=4)

plot(q_last, pch=4, main = TeX(r"($q_{n^{last}}(z;\lambda)$)", bold=TRUE), xlab="", ylab="", ylim=c(-6,6), xlim=c(-6,6))
points(lamb_init[['mu']],col='green', pch=4, lwd=4)
points(lamb_last[['mu']],col='red', pch=4, lwd=4)


# K-HARMONIC
km <- khmeans(rslt[['z']],3)
plot(km[['mu']],col='red', pch=4, xlab="", ylab="", lwd=1,ylim=c(-6,6), xlim=c(-6,6))
points(target[['mu']],col='blue', pch=4, lwd=1)
legend("bottomright", pch=4, col = c("blue", "red"), legend = c("True Center", "Estimated Center"), inset = 0.02)

saveRDS(rslt, "experiment1")

km[['mu']] <- km[['mu']][c(2,3,1),]
table1 <- cbind(target[['mu']], km[['mu']])
rownames(table1) <- c("mu1", "mu2", "mu3")
colnames(table1) <- c("True x_1", "true x2","est x_1", "est x2")

kable(table1,
      booktabs = T,
      escape = F,
      format = "latex",
      digits = 4,
      caption = "Correct Rank %") %>%
  kable_classic(full_width = F, html_font = "Computer Modern typeface", latex_options = c("hold_position"))
##############


# Experiment 2
##############
# Example 2 from https://arxiv.org/PS_cache/arxiv/pdf/0801/0801.1864v1.pdf
# for posterior ratio
set.seed(4833)
p=20
x <- matrix(0, ncol=p)
mu <- rbind(rep(0,p), c(3, rep(0,(p-1))))
sigma <- list(diag(p), 2*diag(p))
mix_weight <- c(0.7, 0.3)
target <- list(mu=mu, sigma=sigma, mix_weight=mix_weight)

N = 2

rslt <- AIMH(x, k=2, om1=0.05, om2=0.15, L=10, a_L=0.1,
             M=20, a_M=0.02, target, N=N, k_init = 3)

x <- rep(0,p)
for (i in 1:N) {
  x <- rbind(x, sim_mix_normal(target))
}
filtered <- rslt[['z']][,1][rslt[['z']][,1]>=-4 & rslt[['z']][,1]<=8]
par(mfrow=c(1,2))
hist(x[,1], breaks=30, main="Target Density Histogram", xlab="")
hist(filtered, breaks=30, main="Simulated Density Histogram", xlab="")
rslt[['A_n']]
mean(rslt[['accept_list']])
##############


# Sanity check functions
########################
par(mfrow=c(1,1))
# sanity check: simulate mixed normal from 2 mixture
x <- sim_2_mix_normal(2000, 0.4, c(-2,0), diag(2), c(2,0), 0.5*diag(2))
plot(x)
x <- sim_2_mix_normal(2000, 0.4, -2, 1, 2, 0.5)
hist(x, breaks=20)

# sanity checking simulate mix normal function
mu <- rbind(c(0,2*sqrt(2)), c(-4,0), c(4,0))
sigma <- list(diag(2), diag(2), diag(2))
mix_weight <- c(0.3,0.3,0.3)
lambda <- list(mu=mu, sigma=sigma, mix_weight=mix_weight)
x <-c()
for (i in 1:1000) {
  x <- rbind(x, sim_mix_normal(lambda))
}
plot(x, pch=4)

# sanity checking k harmonic means function
km <- khmeans(x,3)
points(km[['mu']],col='red', pch=21, lwd=4)

# sanity check bic criterion
for (i in 1:5) {
  km <- khmeans(x,i)
  print(paste("# cluster: ",i, "BIC: ",km[['bic']]))
}

# sanity check q function
mu <- rbind(c(0,-10), c(0,-10))
sigma <- list(diag(2), 25*diag(2))
mix_weight <- c(0.6, 0.4)
lamb_init <- list(mu=mu, sigma=sigma, mix_weight=mix_weight)
mu <- rbind(c(-20,10), c(0,10), c(20,10))
sigma <- list(diag(2), 1.5*diag(2), 2*diag(2))
mix_weight <- c(0.3,0.4,0.3)
lamb_g <- list(mu=mu, sigma=sigma, mix_weight=mix_weight)
x <-c()
for (i in 1:1000) {
  x <- rbind(x, q(lamb_init, lamb_g, k, om1=0.05, om2=0.15)[[1]])
}
plot(x, pch=4, main = TeX(r"(1000 Draws from $q_n(z;\lambda)$)", bold=TRUE), xlab="", ylab="")




