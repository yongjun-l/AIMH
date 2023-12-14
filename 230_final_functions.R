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

library(matrixcalc)
library(MASS)
library(matrixStats)

#' Multivariate normal likelihood
#'
#' @param X pxn matrix
#' @param mu px1 mean vector
#' @param sigma positive definite symmetric matrix
#' @param log boolean for returning log likelihood vs likelihood
#'
#' @return likelihood
eval_mvn_chol <- function(X, mu, sigma, log=FALSE) {
  d <- length(mu)
  n <- length(X)
  L <- t(chol(sigma))
  u <- forwardsolve(L, t(X-mu))
  if (log==FALSE) {
    return( (2*pi)^(-1/2) * 1/prod(diag(L)) * exp( (-1/2)*t(u)%*%u ) )
  } else {
    return( (-1/2)*log(2*pi) - sum(log(diag(L))) - (1/2)*t(u)%*%u )
  }
}



#' K harmonic means clustering
#'
#' @param x nxp matrix.
#' @param c number of clusters
#'
#' @return centers: cxp matrix. each row is mu for mixture gaussian,
#'         V: estimated variance-covariance matrix for each mixture
#'         bic: for choosing the number of clusters.
khmeans <- function(x, c) {
  n = nrow(x)
  p <- ncol(x)
  centers <- matrix(runif(c*p, -1,1), nrow = c) # number of clusters
  V <- list()
  iter <- 0
  while(1) {
    iter <- iter + 1
    d <- as.matrix(dist(rbind(centers, x)))
    d <- d[(c+1):(c+n), 1:c] # euclidean distance from each data point to each center
    d <- d+0.001 # prevent degeneracy
    d1 <- d^(-c-2) # dim: (n,c)

    if (c==1) {
      w <- d1 / d^(-c*2) # weight for each data point. dim: (n, 1)
      m <- as.matrix(d1 / d1) # membership function for each cluster. dim: (n, c)
    } else {
      w <- rowSums(d1) / rowSums(d^(-c))^2
      m <- d1 / rowSums(d1)
    }

    old_centers <- centers
    mw <- as.matrix(m*w)
    if (c==1) {
      colsum_mw <- sum(mw) # dim: (1,c)
    } else {
      colsum_mw <- colSums(mw) # dim: (1,c)
    }

    # update center
    for (j in 1:c) {
      if (c==1) {
        centers[j,] <- colSums(mw[j]*x) / colsum_mw[j]
      } else {
        centers[j,] <- colSums(mw[,j]*x) / colsum_mw[j]
      }
    }

    # stopping criteria
    if (max(abs(old_centers-centers))<1e-6 | iter > 20) {
      break
    }
  } # end of while loop: convergence

  # variance estimation: once after convergence
  for (j in 1:c) {
    resid <- x-matrix(rep(centers[1,],each=n),nrow=n)
    v <- array(0, dim = c(p, p, n))
    for (k in 1:n) {
      v[,,k] <- mw[k,j] * as.matrix(resid[k,]) %*% t(as.matrix(resid[k,]))
    }
    V[[j]] <- rowSums(v, dims = 2) / colsum_mw[j]

    # replace variance matrix if not pos. def.
    if(is.positive.definite(V[[j]], tol=1e-8) == FALSE) {
      V[[j]] <- 0.25*var(x)
    }
  }

  # calculate BIC
  mix_weight = matrix(rep(colMeans(m), times = n), nrow = n, byrow = 1)

  lkhd <- matrix(nrow = n, ncol = c)
  for (j in 1:c) {
    #lkhd[, j] <- lapply(t(x), eval_mvn_chol, centers[j, ], V[[j]], log = 0)
    for (k in 1:n) {
      lkhd[k, j] <- eval_mvn_chol(t(x[k, ]), centers[j, ], V[[j]], log = 0)
    }
  }
  log_likelihood = sum(log(rowMaxs(lkhd )))
  # c vectors of dimension p. weights sum to one so dof is c - 1
  bic <- (p*c + c - 1) * log(n) - 2 * log_likelihood

  return(list("mu"=centers, "sigma"=V, "mix_weight"=colMeans(m), "bic"=bic))
}

#' Pick number of clusters
#'
#' @param z vector of draws
#'
#' @return return the estimated lambda with the lowest BIC criterion
pick_c <- function(z) {
  bic <- c()
  rslt <- list()
  for (i in 1:5) {
    rslt[[i]] <- khmeans(z, i)
    bic <- c(bic, rslt[[i]][["bic"]])
  }
  return(rslt[[which.min(bic)]])
}

#' simulate n observations from mixed multivariate normal with only 2 mixtures
#'
#' @param n no. observations
#' @param a1 mixture weight
#' @param mu1 mean
#' @param sig1 p.d. symm. mat
#' @param mu2 mean
#' @param sig2 p.d. symm. mat
#'
#' @return simulated data
sim_2_mix_normal <- function(n, a1, mu1, sig1, mu2, sig2) {
  x <- c()
  for (i in 1:n) {
    if (runif(1) <a1) {
      x <- rbind(x, mvrnorm(1, mu1, sig1))
    } else {
      x <- rbind(x, mvrnorm(1,mu2,sig2))
    }
  }
  return(x)
}


#' Simulate mixed multivariate normal with c number of mixtures.
#'
#' @param lambda list
#'               vector of mu, dim: (c, p)
#'                 c - number of clusters
#'                 p - dimention of parameter space
#'               sigma - list of p.d. symm var-cov matrices
#'               mix_weight - estimated weight
#'
#' @return 1 simulated observation from the multivar mixture normal
sim_mix_normal <- function(lambda) {
  mu <- lambda[['mu']]
  sigma <- lambda[['sigma']]
  mix_weight <- lambda[['mix_weight']]
  n_mix <- length(mix_weight)
  c <- sample(1:n_mix, size = 1, prob = mix_weight)
  return(mvrnorm(1, mu[c,], sigma[[c]]))
}


#' Evaluate mixture multvar normal
#'
#' @param z one data point in the mixture
#' @param lambda list
#'               vector of mu, dim: (c, p)
#'                 c - number of clusters
#'                 p - dimention of parameter space
#'               sigma - list of p.d. symm var-cov matrices
#'               mix_weight - estimated weight
#'
#' @return evaluated probability density of the mixture.
eval_mix_mvnormal <- function(z, lambda) {
  mu <- lambda[['mu']]
  sigma <- lambda[['sigma']]
  mix_weight <- lambda[['mix_weight']]
  n_mix <- length(mix_weight)
  rslt <- 0
  for (i in 1:n_mix) {
    p <- eval_mvn_chol(z, t(mu[i,]), sigma[[i]], log=FALSE)
    rslt <- rslt + p * mix_weight[i]
  }
  return(rslt)
}


#' simulate 1 observation from proposal density
#'
#' @param lamb_init initial parameters
#' @param lamb_g n_th iteration paramters
#' @param k g_n variance increase factor (scalar)
#' @param om1 weight for initial g_0
#' @param om2 weight for g_n_tilde (with increased variance)
#' @param curr current draw for ratio calculation
#'
#' @return one observation from proposal density
q <- function(lamb_init, lamb_g, k, om1, om2, curr=0) {
  lamb_g_til <- lamb_g
  lamb_g_til[['sigma']] <- lapply(lamb_g_til[['sigma']], function(x) x * k)

  lamb_list <- list(lamb_init, lamb_g_til, lamb_g)
  g_weight <- c(om1, om2, (1-om1-om2))
  selected_g <- sample(1:3, size=1, prob=g_weight)
  lamb <- lamb_list[[selected_g]]

  prop <- sim_mix_normal(lamb)

  # evaluate curr, prop
  eval_curr <- 0
  eval_prop <- 0
  for (i in 1:3) {
    eval_curr <- eval_curr + g_weight[i] * eval_mix_mvnormal(curr, lamb_list[[i]])
    eval_prop <- eval_prop + g_weight[i] * eval_mix_mvnormal(prop, lamb_list[[i]])
  }

  ratio <- eval_curr/eval_prop

  return(list(prop, ratio))
}


#' Adaptive Independent Metropolis Hastings
#'
#' @param x data dim: (n,p)
#' @param k g_n variance multiplication factor
#' @param om1 g_0 weight
#' @param om2 g_n_tilde weight
#' @param L
#' @param a_L
#' @param M
#' @param a_M
#' @param N number of draws
#' @param k_init g_0 variance multiplication factor
#'
#' @return simulated draws
AIMH <- function(x, k, om1, om2, L, a_L, M, a_M, target, N=1000, k_init=16) {
  n_update <- c(seq(50, 400, by=50), seq(500, 1000, by=100), seq(1000,10000, by=1000))
  A_n <- 0 # number of accepted draws
  n_star <- 0 # iteration at the first
  n_last <- 0 # iteration at the end of preliminary stage
  p <- ncol(x) # dimension of the parameter space
  prelim <- TRUE

  # initial lambda
  mu_init <- rbind(rep(0,p), rep(0,p))
  sigma_init <- list(diag(p), k_init*diag(p))
  mix_weight <- c(0.6, 0.4)
  lamb_init <- list(mu=mu_init, sigma=sigma_init, mix_weight=mix_weight)
  z <- t(as.matrix(sim_mix_normal(lamb_init)))
  n <- dim(z)[1]

  # cache
  pi_ratio_list <- c(0,0)
  prop_list <- z
  accept_list <- c()
  q_cache <- list(lamb_init)

  while(1) {
    if (n %% 100 ==0) { print(paste("n", n, "A_n", A_n)) }
    # draw one from q.
    curr <- z[n]
    if (n_star == 0) {
      prop <- sim_mix_normal(lamb_init)
      eval_curr <- eval_mix_mvnormal(curr, lamb_init)
      eval_prop <- eval_mix_mvnormal(prop, lamb_init)
      q_ratio <- eval_curr/eval_prop
    } else {
      q <- q(lamb_init, lamb_g, k, om1, om2, curr)
      prop <- q[[1]]
      q_ratio <- q[[2]]

      if (n == 200 | n==500 | (n_last>500 && n == n_last) | n == (N-1)) {
        q_cache <- append(q_cache, list(list(lamb_init, lamb_g, k, om1, om2)))
      }
    }

    ####################
    # for additive gaussian mixture model
    # lkhd(prop)*prior*prop_den(curr) / lkhd(curr)*prior*prop_den(prop)
    #
    ####################

    ####################
    # mixture mvnormal posterior
    pi_curr <- eval_mix_mvnormal(curr, target)
    pi_prop <- eval_mix_mvnormal(prop, target)
    pi_ratio <- pi_prop / pi_curr
    ####################

    # accept or reject proposal
    accept <- min(pi_ratio*q_ratio, 1)
    accept_list <- c(accept_list, accept)
    if (runif(1)<accept) {
      A_n <- A_n + 1
      z <- rbind(z, prop)
      n <- dim(z)[1]
    } else {
      z <- rbind(z, z[n,])
      n <- dim(z)[1]
    }

    # PRELIMINARY PHASE
    if (prelim == TRUE && A_n >= 5*p) {
      if (n_star ==0) { # INITIAL update
        n_star <- n
        lamb_g <- pick_c(z)
        #lamb_g <- khmeans(z,2)
      } else {
        if ((n-n_star) %in% n_update) { # specified update iteration
          lamb_g <- pick_c(z)
          #lamb_g <- khmeans(z,2)
        }
        # if mean tail accpet_list smaller than a_L
        if (mean(tail(accept_list, n=L))<a_L) {
          lamb_g <- pick_c(z)
          #lamb_g <- khmeans(z,2)
        }
      }
      if (n > M) { # end prelim criteria
        if (min(accept_list[(n-M):(n-1)])>a_M) {
          print("Entering the Second Stage")
          lamb_init <- lamb_g
          prelim <- FALSE
          n_last <- n
        }
      }
    } else { # SECOND PHASE: UPDATE EVERY 1000
      if ( (n-n_star) %% 1000 == 0 ) {
        lamb_g <- pick_c(z)
        #lamb_g <- khmeans(z,2)
      }
    }
    if (n == (N)) {
      break
    }
  } # end of draws

  return(list(z=z,
              target=target,
              accept_list=accept_list,
              q_cache=q_cache,
              A_n=A_n))
}






