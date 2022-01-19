library(readr)
library(dplyr)
library(ggplot2)
library(grid)
library(gridExtra)

# Path to data file
data_path = "~/Code/cognitive/HMM/real data interference/nback_interference_data.csv"

# Directory where Gibbs sampler is stored
code_dir = "~/Code/cognitive/HMM"

# .R file name for Gibbs sampler
code_file = "sampler.R"

# Save MCMC samples?
save = T
if(save){
  # Directory to save MCMC samples
  result_dir = "~/Code/cognitive/HMM/real data interference"
  
  # File name for MCMC samples
  result_file = "cogflex_result.Rdata"
}

# Directory to save plots
plot_dir = "~/Code/cognitive/HMM/real data interference"

# File name for plots
plot_file = "cogflex_plot.pdf"

# Read data and filter out the first 3 trials
data = read_csv(data_path)
data = data %>% filter(trialnum>3)

# Binary outcomes
y = data$correct

# Number of observations
N = length(y)

# Number of observations per individual
J = max(data$trialnum)-3

# Number of individual
P = N/J

# ID of individual for each observation
ID = rep(1:P, each=J)

# Run Gibbs sampler
setwd(code_dir)
source(code_file)
iter       = 60000
state      = list("1"=1:4, "2"=5:8, "3"=9:12)
num_regime = length(state)
prior      = list(psi0 = rep(3, num_regime),
                  psi1 = rep(3, num_regime),
                  m    = rep(0, num_regime),
                  tau2 = rep(1, num_regime),
                  phi0 = 3,
                  phi1 = 3)
burnin     = 10000
n_update   = iter/100
result     = gibbs_sampler(y, ID, iter, state, prior, burnin, n_update)
if(save){
  setwd(result_dir)
  save.image(file=result_file)
}

# Plotting
iter_range = (burnin+1):iter

# Function for plotting trajectory
trajectory = function(beta_i, regime, indv, ID, y, cond, n_cp=NULL, top3=NULL){
  
  # Arguments:
  # beta_i - Posterior samples of beta_i
  # regime - Posterior samples of regime
  # indv   - Individual to plot
  # ID     - Vector containing ID of individuals for each observation
  # y      - Binary outcome
  # cond   - T for plotting conditional trajectory, F for plotting unconditional trajectory
  # n_cp   - Number of change points for starting from different regimes (only when cond=T)
  # top3   - Number of change points with the top 3 highest conditional posterior probability (only when cond=T)
  
  J = nrow(regime)
  
  # Posterior probability of a correct outcome
  prob_correct = sapply(1:dim(beta_i)[2], function(i) pnorm(beta_i[regime[,i],i]))
  
  if(cond==T){ # Conditional trajectory
    
    # Iterations where number of change points has the top 3 highest conditional posterior probability
    index = which(n_cp$n_cp %in% top3)
    
    # Data frame for: n_cp - number of change points, prob - probability of a correct outcome, trial - trial number
    df = data.frame(n_cp=rep(n_cp$n_cp[index], each=J), 
                    prob=as.vector(prob_correct[,index]), 
                    trial=rep(1:J, length(index)))
    
    # Group by number of change points and trial to compute posterior mean and 95% credible interval
    mean_prob_correct = df %>% 
                        group_by(n_cp, trial) %>% 
                        dplyr::summarise(prob=mean(prob))
    ci_prob_correct   = df %>% 
                        group_by(n_cp, trial) %>% 
                        dplyr::summarise(prob_lb=quantile(prob, probs=0.025), prob_ub=quantile(prob, probs=0.975))
    
    # Data frame for: n_cp - number of change points, Trial - trial number, Outcome - binary outcome,
    #                 Mean - posterior mean of probability of a correct outcome,
    #                 lb   - lower bound of 95% credible interval, ub - upper bound of 95% credible interval
    traj = data.frame()
    
    for(i in 1:length(top3)){
      data1 = mean_prob_correct %>% filter(n_cp==top3[i])
      data2 = ci_prob_correct %>% filter(n_cp==top3[i])
      traj  = rbind(traj, data.frame(n_cp=top3[i], 
                                     Trial=1:J, 
                                     Outcome=y[ID==indv], 
                                     Mean=data1$prob, 
                                     lb=data2$prob_lb, 
                                     ub=data2$prob_ub))
    }
  }else{ # Unonditional trajectory
    
    # Posterior mean and 95% credible interval of probability of a correct outcome
    mean_prob_correct = apply(prob_correct, 1, mean)
    ci_prob_correct   = apply(prob_correct, 1, function(x) quantile(x, probs=c(0.025,0.975)))
    
    traj = data.frame(Trial=1:J, 
                      Outcome=y[ID==indv], 
                      Mean=mean_prob_correct, 
                      lb=ci_prob_correct[1,], 
                      ub=ci_prob_correct[2,])
  }

  return(traj)
}

# Function to compute number of change points
num_change_point = function(regime){
  d = diff(regime)
  return(sum(d!=0))
}

# Save plots
setwd(plot_dir)
pdf(plot_file, height=8, width=15)

for(i in 1:P){
  
  # Index of iterations starting from regime 1 and 2
  index = list(which(result$regime[1,i,]==1), which(result$regime[1,i,]==2))
  
  # Posterior samples of regime conditioned on first regime (only include samples after burnin)
  regime = list(as.matrix(result$regime[,i,intersect(index[[1]], iter_range)]),
                as.matrix(result$regime[,i,intersect(index[[2]], iter_range)]))
  
  # Data frame for: init_regime - regime in the first trial, 
  #                 prop - posterior probability of starting in either regime
  init = data.frame(init_regime=as.factor(c(1,2)), 
                    prop=c(ncol(regime[[1]])/(iter-burnin), ncol(regime[[2]])/(iter-burnin)))
  
  p1 = ggplot(init, aes(x=init_regime, y=prop)) + 
       geom_bar(stat="identity") +
       scale_y_continuous(limits=c(0,1)) +
       xlab("Initial Regime") +
       ylab("Probability") +
       theme_bw() +
       theme(axis.title = element_text(size=16),
             axis.text = element_text(size=14))
  
  # Compute number of change points in each iteration 
  if(ncol(regime[[1]])>0 & ncol(regime[[2]])>0){
    n_cp = rbind(data.frame(init_regime="1", 
                            n_cp = apply(regime[[1]], 2, num_change_point)), 
                 data.frame(init_regime="2", 
                            n_cp = apply(regime[[2]], 2, num_change_point)))
  }else if(ncol(regime[[1]])>0 & ncol(regime[[2]])==0){
    n_cp = data.frame(init_regime="1", 
                      n_cp = apply(regime[[1]], 2, num_change_point))
  }else if(ncol(regime[[1]])==0 & ncol(regime[[2]])>0){
    n_cp = data.frame(init_regime="2", 
                      n_cp = apply(regime[[2]], 2, num_change_point))
  }
    
  p2=ggplot(n_cp, aes(x=n_cp)) + 
     stat_count(aes(y=..prop.., group=init_regime)) +
     scale_x_continuous(breaks=0:max(as.numeric(names(table(n_cp$n_cp))))) +
     scale_y_continuous(limits=c(0,1)) +
     facet_wrap(~init_regime) +
     xlab("Number of Change Points") +
     ylab("Conditional Distribution") +
     theme_bw() +
     theme(axis.title = element_text(size=16),
           strip.text = element_text(size=14),
           axis.text = element_text(size=14))
  
  # Top 3 number of change points conditioned on the most probable initial regime
  top3 = table(n_cp$n_cp[n_cp$init_regime==as.character(which.max(init$prop))])
  top3 = as.numeric(names(top3)[order(top3, decreasing=T)[1:3]])
  
  most_probable = which.max(init$prop)
  beta_i        = result$beta_i[,i,intersect(index[[most_probable]], iter_range)]
  n_cp          = n_cp[n_cp$init_regime==as.character(most_probable),]
  traj          = trajectory(beta_i, regime[[most_probable]], i, ID, y, T, n_cp, top3)
  traj$n_cp     = factor(traj$n_cp, levels=top3)
  
  p3=ggplot(traj, aes(x=Trial, y=Outcome)) +
     geom_point() +
     geom_errorbar(aes(ymin=lb, ymax=ub), size=0.25, col="red") +
     facet_grid(~n_cp) +
     geom_point(data=traj, aes(x=Trial, y=Mean), col="red", pch=18) +
     coord_cartesian(ylim=c(0,1)) +
     ylab("Conditional Success Probability") +
     theme_bw()+
     theme(axis.title = element_text(size=16),
           axis.text = element_text(size=14),
           strip.text = element_text(size=14),
           legend.position = "right")
  
  # Unconditional trajectory
  regime = result$regime[,i,iter_range]
  beta_i = result$beta_i[,i,iter_range]
  traj   = trajectory(beta_i, regime, i, ID, y, F)
  
  p4=ggplot(traj, aes(x=Trial, y=Outcome)) +
     geom_point() +
     geom_errorbar(aes(ymin=lb, ymax=ub), size=0.25, col="red") +
     geom_point(data=traj, aes(x=Trial, y=Mean), col="red", pch=18) +
     coord_cartesian(ylim=c(0,1)) +
     ylab("Success Probability") +
     theme_bw()+
     theme(axis.title = element_text(size=16),
           axis.text = element_text(size=14),
           strip.text = element_text(size=14))
  
  grid.arrange(p1, p2, p3, p4, layout_matrix = rbind(c(4,1,2,2), c(3,3,3,3)), top = textGrob(paste0("Individual: ", i),gp=gpar(fontsize=20,font=3)))
}

dev.off()
