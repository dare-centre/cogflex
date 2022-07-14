library(dplyr)
library(LaplacesDemon)
library(BayesLogit)
library(truncnorm)
transition_prob = function(q, state){
  
  # Return a matrix of Markov transition probabilities between states
  
  # Arguments:
  # q     - vector containing transition probabilities between regimes
  # state - list containing indices representing each regime
  
  # Probabilities of remaining in the same state
  num_regime = length(state)
  omega      = unlist(lapply(state, length))
  val        = c()
  for(i in 1:num_regime){
    val = c(val, rep(0, omega[i]-1), 1-q[i])
  }
  Q = diag(val)
  
  # Probabilities of transitioning into a different state
  val = c()
  for(i in 1:(num_regime-1)){
    val = c(val, rep(1, omega[i]-1), q[i])
  }
  val = c(val, rep(1, omega[num_regime]-1))
  
  # Filling up entries for probabilities of transitioning into a different state
  Q[row(Q)-col(Q)==-1]                     = val
  Q[nrow(Q), min(state[[(num_regime-1)]])] = q[num_regime]
  
  return(Q)
}

likelihood_regime = function(y_star, mu, sigma2){
  
  # Return a matrix of likelihood function under each regime
  
  # Arguments:
  # y_star - vector containing posterior propensity scores
  # mu     - vector containing posterior mean of random effects distribution
  # sigma2 - vector containing posterior variance of random effects distribution
  
  J          = length(y_star)
  num_regime = length(mu)
  y_star_rep = matrix(y_star, J, num_regime)
  
  # Unconditional mean and variance for each regime
  mean     = matrix(mu, J, num_regime, T)
  variance = matrix(sigma2, J, num_regime, T) + 1
  
  # Compute likelihood function
  likelihood = dnorm(y_star_rep, mean, sqrt(variance))
  
  return(likelihood)
}

normalise = function(vector) return(vector/sum(vector))

gibbs_state = function(likelihood, Q, xi, state){
  
  # Return a vector of posterior state for an individual
  
  # Arguments:
  # likelihood - matrix containing likelihood function under each state
  # Q          - matrix containing Markov transition probabilities between states
  # xi         - scalar for the probability of starting from regime 1 in the first trial
  # state      - list containing indices representing each regime
  
  J         = nrow(likelihood)
  num_state = ncol(Q)
  states    = vector("numeric", J)
  
  # Compute predictive probabilities
  pred_prob                     = matrix(0, J, num_state)
  pred_prob[1,1]                = xi
  pred_prob[1, min(state[[2]])] = 1-xi
  pred_prob[1,]                 = normalise(pred_prob[1,] * likelihood[1,])
  for(i in 2:J){
    pred_prob[i,] = normalise((pred_prob[i-1,] %*% Q) * likelihood[i,])
  }

  # Posterior sampling of state in the final trial
  states[J] = sample(1:num_state, 1, F, pred_prob[J,])
  
  # Backwards sampling
  for(i in (J-1):1){
    prob      = normalise(pred_prob[i,] * Q[,states[i+1]])
    states[i] = sample(1:num_state, 1, F, prob)
  }
  
  return(states)
}

regime_transform = function(states_sample, state){
  
  # Return a vector of posterior regime from a vector of posterior state
  
  # Arguments:
  # states_sample - vector containing posterior samples of states
  # state         - list containing indices representing each regime
  
  num_regime = length(state)
  regime     = vector("numeric", length(states_sample))
  
  # Transforming states into regimes
  for(i in 1:num_regime){
    regime[states_sample %in% state[[i]]] = i
  }
  
  return(regime)
}

count_transition = function(regime, state){
  
  # Return a data frame counting number of successful transitions into different regimes
  
  # Arguments:
  # regime - vector containing posterior samples of regimes
  # state  - list containing indices representing each regime
  
  num_regime = length(state)
  J          = length(regime)
  
  # Minimum number of trials in each regime
  omega = unlist(lapply(state, length))
  
  # Consecutive pairs of regimes
  regime_pair = cbind(regime[1:(J-1)], regime[2:J])
  
  # Data frame storing number of successful transitions (freq) from regime x1 to x2 
  # out of N attempts
  n_trans = data.frame(x1=1:num_regime, x2=c(2:num_regime, num_regime-1), freq=NA, N=NA)
  
  for(i in 1:num_regime){
    index = which(regime_pair[,1]==i)
    if(length(index)>0){
      # Create a list containing consecutive trial numbers in state i
      g = cumsum(c(1, abs(index[-length(index)]-index[-1])>1))
      g = by(index, g, identity)
      
      # Counting effective number of attempts by removing (minimum number of trial-1)
      g = unlist(lapply(g, length))
      g = g-omega[i]+1
      
      # Removing negative effective number of attempts (this could happen towards the end of block)
      g = sum(g[g>0])
      
      n_trans$N[i]    = g
      n_trans$freq[i] = sum(regime_pair[index,2]!=i)
    }else{
      n_trans$N[i]    = 0
      n_trans$freq[i] = 0
    }
  }
  
  return(n_trans)
}

gibbs_q = function(count_trans, psi0, psi1){
  
  # Return a vector of posterior transition probabilities between regimes
  
  # Arguments:
  # count_trans - data frame containing number of successful transitions and effective attempts
  # psi0        - vector containing first shape parameter of Beta prior
  # psi1        - vector containing second shape parameter of Beta prior
  
  # Prior distribution - Beta(psi0, psi1)
  
  num_trans  = length(unique(count_trans$x1))
  total_freq = apply(matrix(count_trans$freq, byrow=T, ncol=num_trans), 2, sum)
  total_N    = apply(matrix(count_trans$N, byrow=T, ncol=num_trans), 2, sum)
  q          = rbeta(num_trans, psi0+total_freq, psi1+total_N-total_freq)
  
  return(q)
}

gibbs_mu = function(beta_i, inv_sigma2, m, tau2){
  
  # Return a vector of posterior mean of random effects distribution
  
  # Arguments:
  # beta_i     - matrix containing posterior random effects for all individuals and regimes
  # inv_sigma2 - vector containing posterior precision of random effects distribution
  # m          - vector containing mean of Gaussian prior
  # tau2       - vector containing variance of Gaussian prior
  
  # Prior distribution - Gaussian(m, tau2)
  
  num_regime = length(m)
  mu         = vector("numeric", num_regime)
  
  for(i in 1:num_regime){
    index = which(!is.na(beta_i[i,]))
    var   = 1/(length(index)*inv_sigma2[i] + 1/tau2[i])
    mean  = var * (sum(beta_i[i,index])*inv_sigma2[i] + m[i]/tau2[i])
    mu[i] = rnorm(1, mean, sqrt(var))
  }
  
  return(mu)
}

gibbs_beta_i = function(y_star, regime, mu, inv_sigma2, beta_i){
  
  # Return a vector of posterior random effects for an individual
  
  # Arguments:
  # y_star     - vector containing posterior propensity scores
  # regime     - vector containing posterior samples of regimes
  # mu         - vector containing posterior mean of random effects distribution
  # inv_sigma2 - vector containing posterior precision of random effects distribution
  # beta_i     - vector containing current sampled values of posterior random effects
  
  num_regime = length(inv_sigma2)
  
  for(i in 1:num_regime){
    group_i = regime==i
    if(sum(group_i)==0){
      # No value is sampled if no trial is spent on regime i
      beta_i[i] = NA
    }else{
      var  = 1/(sum(group_i) + inv_sigma2[i])
      mean = var * (sum(y_star[group_i]) + mu[i]*inv_sigma2[i])
      if(i %in% c(1,3)){
        if(is.na(beta_i[2])){
          # Case: regime 1 and beta_2 is null
          # Note: not possible for regime 3 and beta_2 is null because transition to regime 3 can
          #       only occur from regime 2
          beta_i[i] = rnorm(1, mean, sqrt(var))
        }else{
          # Case: regime 1 or 3 and beta_2 is not null
          # beta_i[i] = rtruncnorm(1, beta_i[2], Inf, mean, sqrt(var))
          beta_i[i] = rtruncnorm(1, -Inf, Inf, mean, sqrt(var))
        }
      }else{
        if(all(is.na(c(beta_i[1], beta_i[3])))){
          # Case: regime 2 and both beta_1 and beta_3 are null => stay in regime 2 for entire block
          beta_i[i] = rnorm(1, mean, sqrt(var))
        }else{
          # Case: regime 2 and at least one of beta_1 and beta_3 are not null
          # beta_i[i] = rtruncnorm(1, -Inf, min(c(beta_i[1], beta_i[3]), na.rm=T), mean, sqrt(var))
          beta_i[i] = rtruncnorm(1, -Inf, Inf, mean, sqrt(var))
        }
      }
    }
  }
  
  return(beta_i)
}

gibbs_sigma2 = function(beta_i, mu, inv_sigma2, prior_scale=5){
  
  # Return a vector of posterior variance of random effects distribution
  
  # Arguments:
  # beta_i      - matrix containing posterior random effects for all individuals and regimes
  # mu          - vector containing posterior mean of random effects distribution
  # inv_sigma2  - vector containing posterior precision of random effects distribution
  # prior_scale - scalar for scale parameter of a Cauchy distribution
  
  # Prior distribution - half-Cauchy
  
  num_regime = length(inv_sigma2)
  new_sigma2 = vector("numeric", num_regime)
  
  for(i in 1:num_regime){
    a             = rinvgamma(1, 1, inv_sigma2[i]+1/prior_scale^2)
    beta_i_k      = beta_i[i,!is.na(beta_i[i,])]
    sq_diff       = sum((beta_i_k-mu[i])^2)
    new_sigma2[i] = rinvgamma(1, (length(beta_i_k)+1)/2, 1/a+1/2*sq_diff)
  }
  
  return(new_sigma2)
}

gibbs_y_star = function(beta_i_all, lb, ub){
  
  # Return a vector of posterior propensity scores
  
  # Arguments:
  # beta_i_all - vector containing posterior random effects for all observations
  # lb         - vector containing lower bounds of truncated normal distribution
  # ub         - vector containing upper bounds of truncated normal distribution
    
  y_star = rtruncnorm(1, lb, ub, beta_i_all, 1)
  
  return(y_star)
}

gibbs_sampler = function(y, ID, iter, state, prior=NULL, burnin=NULL, n_update=NULL,
  seed = 1){
  
  # Arguments:
  # y        - vector containing binary outcomes
  # ID       - vector containing ID of individuals for each observation
  # iter     - scalar for number of MCMC iterations
  # state    - list containing indices representing each regime
  # prior    - list containing hyperparameters for prior distributions
  # burnin   - scalar for number of MCMC burnin (for trace plot only)
  # n_update - scalar for number of iterations for each update of MCMC progress
  
  set.seed(seed)
  
  if(is.null(burnin)){
    burnin = ceiling(iter/5)
  }
  if(is.null(n_update)){
    n_update = ceiling(iter/100)
  }
  
  # Number of individuals
  P = max(ID)
  
  # Number of trials in a block
  J = unique(tabulate(ID))
  
  # Total number of observations
  N = length(y)
  
  # Number of regime
  num_regime = length(state)
  
  # Minimum number of trials in each regime
  omega = unlist(lapply(state, length))
  
  # Default values of hyperparameters for prior
  if(is.null(prior)){
    prior = list(psi0 = rep(3, num_regime),
      psi1 = rep(3, num_regime),
      m    = rep(0, num_regime),
      tau2 = rep(1, num_regime),
      phi0 = 3,
      phi1 = 3)
  }
  
  # Initialising values for parameters
  
  # Regime
  regime = array(NA, dim=c(J, P, iter))
  for(i in 1:P){
    regime[,i,1] = sort(sample(1:(sample(1:num_regime, 1)), J, replace=T))
  }
  
  # Transition probability
  q     = matrix(NA, num_regime, iter)
  q[,1] = rbeta(num_regime, prior$psi0, prior$psi1)
  
  # Mean of random effects distribution
  mu     = matrix(NA, num_regime, iter)
  mu[,1] = rnorm(num_regime, prior$m, sqrt(prior$tau2))
  
  # Variance of random effects distribution
  sigma2     = matrix(NA, num_regime, iter)
  sigma2[,1] = abs(rcauchy(num_regime, 0, 5))
  
  # Random effects
  beta_i = array(NA, dim=c(num_regime, P, iter))
  for(i in 1:P){
    beta_i[1:max(regime[,i,1]),i,1] = rnorm(max(regime[,i,1]), mu[1:max(regime[,i,1]),1],
      sqrt(sigma2[1:max(regime[,i,1]),1]))
  }
  regime_vec = as.vector(regime[,,1])
  beta_i_all = sapply(1:N, function(j) beta_i[regime_vec[j], ID[j], 1])
  
  # Propensity scores
  y_star     = matrix(NA, N, iter)
  y_star[,1] = abs(rnorm(N, beta_i_all, 1))*(2*y-1)
  
  # Probability of starting from regime 1 in the first trial
  xi    = vector("numeric", iter)
  xi[1] = rbeta(1, prior$phi0, prior$phi1)
  
  # Lower and upper bounds of truncated normal distributions
  lb = ub = vector("numeric", N)
  lb[y==0] = -Inf
  ub[y==1] = Inf
  
  for(i in 2:iter){
    
    states_val  = sapply(1:P, function(j){
      likelihood = likelihood_regime(y_star[((j-1)*J+1):(j*J),i-1], mu[,i-1], 
        sigma2[,i-1])
      Q = transition_prob(q[,i-1], state)
      gibbs_state(likelihood[,rep(1:num_regime, omega)], Q, xi[i-1], state)})
    regime[,,i] = apply(states_val, 2, regime_transform, state=state)
    
    n1    = sum(regime[1,,i]==1)
    xi[i] = rbeta(1, prior$phi0+n1, prior$phi1+P-n1)
    
    count_trans = do.call("rbind", apply(regime[,,i], 2, count_transition, state=state))
    q[,i]       = gibbs_q(count_trans, prior$psi0, prior$psi1)
    
    inv_sigma2   = 1/sigma2[,i-1]
    beta_i[,,i] = sapply(1:P, function(j){
      gibbs_beta_i(y_star[((j-1)*J+1):(j*J),i-1], regime[,j,i], mu[,i-1], inv_sigma2, 
        beta_i[,j,i-1])})
    
    mu[,i] = gibbs_mu(beta_i[,,i], inv_sigma2, prior$m, prior$tau2)
    
    sigma2[,i] = gibbs_sigma2(beta_i[,,i], mu[,i], inv_sigma2)
    
    regime_vec = as.vector(regime[,,i])
    beta_i_all = sapply(1:N, function(j) beta_i[regime_vec[j], ID[j], i])
    y_star[,i] = gibbs_y_star(beta_i_all, lb, ub)
    
    if(i%%n_update==0){
      print(paste('Iteration', i, 'of', iter))
      if(i>burnin){
        par(mfrow=c(num_regime, num_regime))
        for(k in 1:num_regime){
          plot(1:(i-burnin), mu[k,(burnin+1):i], "l", xlab="Iteration", ylab="Value", 
            main=bquote(mu[.(k)]))
        }
        
        for(k in 1:num_regime){
          plot(1:(i-burnin), q[k,(burnin+1):i], "l", xlab="Iteration", ylab="Value", 
            main=bquote(q[.(k)]))
        }
        
        for(k in 1:num_regime){
          plot(1:(i-burnin), sigma2[k,(burnin+1):i], "l", xlab="Iteration", ylab="Value", 
            main=bquote(sigma[.(k)]^2))
        }
      }
    }
  }
  
  return(list(regime=regime, q=q, mu=mu, beta_i=beta_i, sigma2=sigma2, y_star=y_star, 
    xi=xi))
}
