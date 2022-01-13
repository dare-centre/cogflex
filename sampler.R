library(dplyr)
library(LaplacesDemon)
library(BayesLogit)
library(truncnorm)
transition_prob = function(up_prob, state){
  
  # Return a matrix of Markov transition probabilities between states
  
  # Arguments:
  # up_prob - vector containing transition probabilities between regimes
  # state   - list containing indices representing each regime
  
  # Probabilities of remaining in the same state
  num_regime = length(state)
  min_trial  = unlist(lapply(state, length))
  val        = c()
  for(i in 1:num_regime){
    val = c(val, rep(0, min_trial[i]-1), 1-up_prob[i])
  }
  prob_mat = diag(val)
  
  # Probabilities of transitioning into a different state
  val = c()
  for(i in 1:(num_regime-1)){
    val = c(val, rep(1, min_trial[i]-1), up_prob[i])
  }
  val = c(val, rep(1, min_trial[num_regime]-1))
  
  # Filling up entries for probabilities of transitioning into a different state
  prob_mat[row(prob_mat)-col(prob_mat)==-1]              = val
  prob_mat[nrow(prob_mat), min(state[[(num_regime-1)]])] = up_prob[num_regime]
  
  return(prob_mat)
}

likelihood_regime = function(y_star, mu, Sigma){
  
  # Return a matrix of likelihood function under each regime
  
  # Arguments:
  # y_star - vector containing posterior propensity scores
  # mu     - vector containing posterior mean of intercepts
  # Sigma  - vector containing posterior variance of zero-mean random effects
  
  N          = length(y_star)
  num_regime = length(mu)
  y_star_rep = matrix(y_star, nrow=N, ncol=num_regime)
  
  # Mean and variance of each regime
  mean     = matrix(mu, byrow=T, nrow=N, ncol=num_regime)
  variance = matrix(Sigma, byrow=T, nrow=N, ncol=num_regime) + 1
  
  # Compute likelihood function
  likelihood = dnorm(y_star_rep, mean=mean, sd=sqrt(variance))
  
  return(likelihood)
}

normalise = function(vector) return(vector/sum(vector))

gibbs_regime = function(likelihood, prob_mat, init_prob, state){
  
  # Return a vector of posterior state for an individual
  
  # Arguments:
  # likelihood - matrix containing likelihood function under each regime
  # prob_mat   - matrix containing Markov transition probabilities between states
  # init_prob  - scalar for the probability of starting from regime 1 in the first trial
  # state      - list containing indices representing each regime
  
  N         = nrow(likelihood)
  num_state = ncol(prob_mat)
  states    = vector("numeric", length=N)
  
  # Compute predictive probabilities
  pred_prob                     = matrix(0, nrow=N, ncol=num_state)
  pred_prob[1,1]                = init_prob
  pred_prob[1, min(state[[2]])] = 1-init_prob
  for(i in 2:N){
    pred_prob[i,] = normalise((pred_prob[i-1,] %*% prob_mat) * likelihood[i,])
  }
  
  # Posterior sampling of state in the final trial
  states[N] = sample(1:num_state, 1, prob=pred_prob[N,])
  
  # Backwards sampling
  for(i in (N-1):1){
    prob      = normalise(pred_prob[i,] * prob_mat[,states[i+1]])
    states[i] = sample(1:num_state, 1, F, prob)
  }
  
  return(states)
}

regime_transform = function(states_sample, state){
  
  # Return a vector of posterior regime from a vector of posterior state
  
  # Arguments:
  # states_sample - vector containing posterior samples of states
  # state         - list containing indices representing each regime
  
  regime = rep(NA, length=length(states_sample))
  
  # Transforming states into regimes
  for(i in 1:length(state)){
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
  N          = length(regime)
  
  # Minimum number of trials in each regime
  min_trial = unlist(lapply(state, length))
  
  # Consecutive pairs of regimes
  regime_pair = cbind(regime[1:(N-1)], regime[2:N])
  
  # Data frame storing number of successful transitions (freq) from regime x1 to x2 
  # out of N attempts
  n_trans = data.frame(x1=1:num_regime, x2=c(2:num_regime, num_regime-1), freq=NA, N=NA)
  
  for(i in 1:num_regime){
    index = which(regime_pair[,1]==i)
    if(length(index)>0){
      g = cumsum(c(1, abs(index[-length(index)]-index[-1])>1))
      g = by(index, g, identity)
      g = unlist(lapply(g, length))
      g = g-min_trial[i]+1
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

gibbs_gamma = function(count_trans, gamma, prior_prec){
  
  # Return a vector of posterior regression coefficients in the Markov transition
  # probability
  
  # Arguments:
  # count_trans - data frame containing number of successful transitions
  # gamma       - vector containing previous value of posterior regression coefficients
  # prior_prec  - vector containing precision of prior distribution
  
  # Prior distribution - Gaussian(0, 1/prior_prec)
  
  # Reference: Polson et al. (2013). Bayesian Inference for Logistic Models Using 
  #            Pólya–Gamma Latent Variables
  
  num_regime = length(gamma)
  N          = nrow(count_trans)/num_regime
  new_gamma  = matrix(NA, nrow=1, ncol=num_regime)
  freq_i     = matrix(count_trans$freq, byrow=T, ncol=num_regime)
  N_i        = matrix(count_trans$N, byrow=T, ncol=num_regime)
  psi_gamma  = matrix(gamma, byrow=T, nrow=N, ncol=num_regime)
  
  for(i in 1:num_regime){
    jump          = N_i[,i]>0
    omega         = rpg(sum(jump), N_i[jump,i], psi_gamma[jump,i])
    var           = 1/(sum(omega) + prior_prec[i])
    mean          = var * sum(freq_i[jump,i] - N_i[jump,i]/2)
    new_gamma[,i] = rnorm(1, mean, sqrt(var))
  }
  
  return(new_gamma)
}

gibbs_mu = function(y_star, regime_vec, beta_i_all, prior_prec){
  
  # Return a vector of posterior mean of intercepts
  
  # Arguments:
  # y_star     - vector containing posterior propensity scores
  # regime_vec - vector containing posterior regime of all observations
  # beta_i_all - vector containing posterior zero-mean random effects of all observations
  # prior_prec - vector containing precision of prior distribution
  
  # Prior distribution - Gaussian(0, 1/prior_prec)
  
  num_regime = length(prior_prec)
  mu         = matrix(NA, nrow=1, ncol=num_regime)
  fe         = y_star - beta_i_all
  
  for(i in 1:num_regime){
    index  = regime_vec==i
    var    = 1/(sum(index) + prior_prec[i])
    mean   = var * sum(fe[index])
    mu[,i] = rnorm(1, mean, sqrt(var))
  }
  
  return(mu)
}

gibbs_beta_i = function(y_star, regime, mu, inv_Sigma){
  
  # Return a vector of posterior zero-mean random effects for an individual
  
  # Arguments:
  # y_star    - vector containing posterior propensity scores
  # regime    - vector containing posterior samples of regimes
  # mu        - vector containing posterior mean of intercepts
  # inv_Sigma - vector containing posterior precision of zero-mean random effects
  
  num_regime = length(inv_Sigma)
  beta_i     = matrix(NA, nrow=1, ncol=num_regime)
  re         = y_star - mu[regime]
  
  for(i in 1:num_regime){
    group_i = regime==i
    if(sum(group_i)==0){
      beta_i[,i] = NA
    }else{
      var        = 1/(sum(group_i) + inv_Sigma[i])
      mean       = var * sum(re[group_i])
      beta_i[,i] = rnorm(1, mean, sqrt(var))
    }
  }
  
  return(beta_i)
}

gibbs_Sigma = function(beta_i, inv_Sigma, prior_scale=5){
  
  # Return a vector of posterior variance of zero-mean random effects
  
  # Arguments:
  # beta_i      - matrix containing posterior zero-mean random effects for all individuals
  # inv_Sigma   - vector containing previous value of posterior precision of zero-mean random effects
  # prior_scale - scalar for scale parameter of a Cauchy distribution
  
  # Prior distribution - half-Cauchy
  
  num_regime = length(inv_Sigma)
  new_Sigma  = array(NA, dim=c(1,1,num_regime))
  
  for(i in 1:num_regime){
    a              = rinvgamma(1, 1, inv_Sigma[i]+1/prior_scale^2)
    beta_i_tmp     = beta_i[i,!is.na(beta_i[i,])]
    outer_prod     = sum(beta_i_tmp^2)
    new_Sigma[,,i] = rinvgamma(1, (length(beta_i_tmp)+1)/2, 1/a+1/2*outer_prod)
  }
  
  return(new_Sigma)
}

gibbs_y_star = function(regime_vec, beta_i_all, mu, lb, ub){
  
  # Return a vector of posterior propensity scores
  
  # Arguments:
  # regime_vec - vector containing posterior regime of all observations
  # beta_i_all - vector containing posterior zero-mean random effects of all observations
  # mu         - vector containing posterior mean of intercepts
  # lb         - vector containing lower bounds of truncated normal distribution
  # ub         - vector containing upper bounds of truncated normal distribution
  
  y_star = rtruncnorm(1, lb, ub, mu[regime_vec]+beta_i_all, 1)
  
  return(y_star)
}

gibbs_sampler = function(y, ID, iter, burnin, n_update, state){
  
  # Arguments:
  # y        - vector containing binary outcomes
  # ID       - vector containing ID of individuals for each observation
  # iter     - scalar for number of MCMC iterations
  # burnin   - scalar for number of MCMC burnin (for trace plot only)
  # n_update - scalar for number of iterations for each update of MCMC progress
  # state    - list containing indices representing each regime
  
  set.seed(54321)
  
  # Number of individuals
  P = max(ID)
  
  # Number of trials in a block
  J = unique(tabulate(ID))
  
  # Total number of observations
  N = length(y)
  
  # Number of regime
  num_regime = length(state)
  
  # Minimum number of trials in each regime
  min_trial = unlist(lapply(state, length))
  
  # Number of predictor in regression and Markov transition probability (only include intercept here)
  n_psi = 1
  n_X   = 1
  
  # Initialising values for parameters
  
  # Regime
  regime = array(NA, dim=c(J, P, iter))
  for(i in 1:P){
    regime[,i,1] = sort(sample(1:(sample(1:num_regime, 1)), J, replace=T))
  }
  
  # Regression coefficient in Markov transition probability
  gamma          = array(NA, dim=c(n_psi, num_regime, iter))
  gamma[,,1]     = rnorm(n_psi*num_regime)
  prior_prec_psi = rep(1/(pi^2/(3*dim(gamma)[1])), num_regime) # Precision of prior distribution
  
  # Intercepts
  mu           = array(NA, dim=c(n_X, num_regime, iter))
  mu[,,1]      = rnorm(n_X*num_regime)
  prior_prec_X = rep(n_X, num_regime) # Precision of prior distribution
  
  # Zero-mean random effects
  beta_i = array(NA, dim=c(n_X, num_regime, P, iter))
  for(i in 1:P){
    beta_i[,1:max(regime[,i,1]),i,1] = rnorm(max(regime[,i,1]))
  }
  
  # Variance of zero-mean random effects
  Sigma       = array(NA, dim=c(n_X, n_X, num_regime, iter))
  Sigma[,,,1] = diag(1, n_X)
  
  # Propensity scores
  y_star     = array(NA, dim=c(N, iter))
  y_star[,1] = abs(rnorm(N))*(2*y-1)
  
  # Probability of starting from regime 1 in the first trial
  init_prob    = vector("numeric", length=iter)
  init_prob[1] = runif(1)
  
  # Lower and upper bounds of truncated normal distributions
  lb = ub = vector("numeric", length=N)
  lb[y==0] = -Inf
  ub[y==1] = Inf
  
  for(i in 2:iter){
    
    
    probability = exp(matrix(gamma[,,i-1], byrow=T, nrow=P, ncol=num_regime))
    probability = probability/(1+probability)
    tmp_regime  = sapply(1:P, function(j){
                  likelihood = likelihood_regime(y_star[((j-1)*J+1):(j*J),i-1], mu[,,i-1], 
                                                 Sigma[,,,i-1])
                  prob_mat = transition_prob(probability[j,], state)
                  gibbs_regime(likelihood[,rep(1:num_regime, min_trial)], prob_mat, init_prob[i-1], 
                               state)})
    regime[,,i] = apply(tmp_regime, 2, regime_transform, state=state)
    
    init_prob[i] = rbeta(1, 1+sum(regime[1,,i]==1), 1+sum(regime[1,,i]==2))
    
    count_trans = do.call("rbind", apply(regime[,,i], 2, count_transition, state=state))
    gamma[,,i]  = gibbs_gamma(count_trans, gamma[,,i-1], prior_prec_psi)
    
    inv_Sigma    = 1/Sigma[,,,i-1]
    beta_i[,,,i] = sapply(1:P, function(j){
                   gibbs_beta_i(y_star[((j-1)*J+1):(j*J),i-1], regime[,j,i], mu[,,i-1], inv_Sigma)})
    
    regime_vec = as.vector(regime[,,i])
    beta_i_all = sapply(1:N, function(j) beta_i[, regime_vec[j], ID[j], i])
    mu[,,i]    = gibbs_mu(y_star[,i-1], regime_vec, beta_i_all, prior_prec_X)

    Sigma[,,,i] = gibbs_Sigma(beta_i[,,,i], inv_Sigma)
    
    y_star[,i] = gibbs_y_star(regime_vec, beta_i_all, mu[,,i], lb, ub)
    
    if(i%%n_update==0){
      print(paste('Iteration', i, 'of', iter))
      if(i>burnin){
        par(mfrow=c(num_regime, num_regime))
        for(k in 1:num_regime){
          plot(1:(i-burnin), mu[1,k,(burnin+1):i], "l", xlab="Iteration", ylab="Value", 
               main=bquote(mu[.(k)]))
        }
        
        for(k in 1:num_regime){
          plot(1:(i-burnin), gamma[1,k,(burnin+1):i], "l", xlab="Iteration", ylab="Value", 
               main=bquote(gamma[.(k)]))
        }
        
        for(k in 1:num_regime){
          plot(1:(i-burnin), Sigma[1,1,k,(burnin+1):i], "l", xlab="Iteration", ylab="Value", 
               main=bquote(sigma[.(k)]^2))
        }
      }
    }
  }
  
  return(list(regime=regime, gamma=gamma, mu=mu, beta_i=beta_i, Sigma=Sigma, y_star=y_star, 
              init_prob=init_prob))
}
