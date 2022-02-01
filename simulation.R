# Number of individual
P = 80

# Number of observations per individual
J = 57

# ID of individual for each observation
ID = rep(1:P, each=J)

# True value of parameters
q      = c(0.05, 0.35, 0.5)
mu     = c(1.15, 0.1, 0.85)
sigma2 = c(0.02, 0.025, 0.015)
xi     = 0.65

# Minimum number of trials in each regime
omega = c(4,4,4)

# Possible transitions
transition = matrix(c(1,2,2,3,3,2), byrow=T, ncol=2)

# Simulate regime for each individual based on transition probabilities
regime = matrix(NA, nrow=J, ncol=P)
set.seed(12345)
for(i in 1:P){
  # Sample regime for the first trial and set the following trials to be the same using omega
  regime_i = sample(transition[1,], 1, prob=c(xi, 1-xi))
  regime_i = rep(regime_i, omega[regime_i])
  n        = length(regime_i)
  while(n<J){
    # Sample regime for the next trial
    regime_next_trial = sample(transition[regime_i[n],], 1, prob=c(1-q[regime_i[n]], q[regime_i[n]]))
    if(regime_next_trial==regime_i[n]){
      regime_i = c(regime_i, regime_next_trial)
    }else{
      regime_i = c(regime_i, rep(regime_next_trial, omega[regime_next_trial]))
    }
    n = length(regime_i)
  }
  regime[,i] = regime_i[1:J]
}

# Simulate random effects
beta_i = matrix(NA, nrow=3, ncol=P)
for(i in 1:P){
  beta_i[,i] = rnorm(3, mu, sqrt(sigma2))
}
# Check constraints for random effects
all(beta_i[1,]>beta_i[2,] & beta_i[3,]>beta_i[2,])

# Simulate propensity scores
y_star = matrix(NA, nrow=J, ncol=P)
for(i in 1:P){
  y_star[,i] = beta_i[regime[,i],i] + rnorm(J)
}

# Obtain binary outcome from propensity scores
y = as.numeric(y_star > 0)
