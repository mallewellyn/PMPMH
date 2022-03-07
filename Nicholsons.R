source("run_MH.R")
source("gen_blocks.R")


set.seed(66)
n=100
N=numeric(n)
R=numeric(n)
S=numeric(n)
y=numeric(n)
tau=5
P=50
N0=50
delta=0.9
phi=1
sigma2.e=1
sigma2.epsilon=1
e=rgamma(n, 1/sigma2.e, rate=1/sigma2.e)
epsilon=rgamma(n, 1/sigma2.epsilon, rate=1/sigma2.epsilon)
#e=100
#epsilon=10

S[1]=rbinom(1, N0, exp(-delta*epsilon[1]))
N[1]=S[1]
for(t in 2:(tau)){
  S[t]=rbinom(1, N[t-1], exp(-delta*epsilon[t]))
  N[t]=S[t]
}
S[tau+1]=rbinom(1, N[tau], exp(-delta*epsilon[tau+1]))
R[tau+1]=rpois(1, P*N0*e[tau+1]*exp(-1))
N[tau+1]=S[tau+1]+R[tau+1]
for(i in (tau+2):n){
  R[i]=rpois(1, P*N[i-tau-1]*e[i]*exp(-N[i-tau-1]/N0))
  S[i]=rbinom(1, N[i-1], exp(-delta*epsilon[i]))
  N[i]=R[i]+S[i]
}
for(i in 1:n){
  y[i]=rpois(1, phi*N[i])
}


Nits=5e4
l=3
N=25
q1=0.1
qN=1-0.1
sigma=10
delta.e=2
thresh=1/100
end.prop.var=10
l=4
states.init=rep(1, length(y))
theta.init<-rep(0.5, 4)

################ user-defined quantile function
quantile_func<-function(perc, t, states_curr, y, theta, sigma){
  quants=qnorm(perc, states_curr[t], sqrt(sigma))
  unique(quants)
}
midpoint_func<-function(q, delta.e, t){
  N.adapt<-length(q)+1

  mp<-numeric(N.adapt)
  len.cell<-numeric(N.adapt)


  mp[1]<-q[1]-delta.e/2
  mp[N.adapt]<-q[length(q)] + delta.e/2
  len.cell[c(1, N.adapt)]<-delta.e


  mp[-c(1, N.adapt)]=(q[-1]+q[-(length(q))])/2
  len.cell=c(delta.e, q[-1]-q[-length(q)], delta.e)

  return(list("N.adapt"=N.adapt, "mp"=mp, "len.cell"=len.cell))
}

midpoint_int_func_system<-function(mpoints, bin.len, theta, t){
  if(t==1){
    dens=sapply(mpoints[[1]], log_state_dens, t=1, theta=theta)
  } else {
    dens=matrix(rep(0, length(mpoints[[t-1]])*length(mpoints[[t]])), ncol=length(mpoints[[t]]))
    for(j in 1:length(mpoints[[t-1]])){
      for(i in 1:length(mpoints[[t]])){
        temp<-numeric(t)
        temp[t-1]=mpoints[[t-1]][j]
        temp[t]=mpoints[[t]][i]
        dens[j,i]=log_state_dens(temp, t, theta)
      }
    }
  }

  mass=dens + log(bin.len[[t]])
  mass

}
midpoint_int_func_obs<-function(y, mpoints, bin.len, theta, t){
  mass=numeric(length(mpoints[[t]]))
  for(i in 1:length(mpoints[[t]])){
    mass[i]=log_obs_dens(y, c(numeric(t-1), mpoints[[t]][i]), t, theta)
  }
  mass + log(bin.len[[t]])
}
########################## scheme for updating theta
prop_lim<-c(1, 1, 1, 1)
beta.sigma=2
gamma.sigma=2
beta.small=2
gamma.small=2
beta.large=2
gamma.large=700
log_prior_density_sigmasmall<-function(val, beta.small, gamma.small){
  dinvgamma(val, shape=beta.small, rate=gamma.small, log=TRUE)
}
log_prior_density_sigmalarge<-function(val, beta.large, gamma.large){
  dinvgamma(val, shape=beta.large, rate=gamma.large, log=TRUE)
}
theta_update<-function(states, theta){
  theta_curr<-theta
  theta_prop<-theta
  pacc<-numeric(length(theta))

  ## first p
  theta_prop[3]<-runif(1, theta_curr[3]-prop_lim[3], theta_curr[3]+prop_lim[3])
  if(theta_prop[3]>0 && theta_prop[3]<1){
    pacc[3]<-min(1, exp( sum(sapply(seq(1, len.y, by=1), log_state_dens, states=states, theta=theta_curr)) -
                           sum(sapply(seq(1, len.y, by=1), log_state_dens, states=states, theta=theta_prop))))

  }
  if(runif(1)<pacc[3]){
    theta_curr[3]<-theta_prop[3]
  } else {
    theta_prop[3]<-theta_curr[3]
  }

  ## sigma eta small
  theta_prop[1]<-runif(1, theta_curr[1]-prop_lim[1], theta_curr[1]+prop_lim[1])
  if(theta_prop[1]>0 && theta_prop[1]<theta_curr[2]){
    pacc[1]<-min(1, exp(log_prior_density_sigmasmall(theta_prop[1], beta.small, gamma.small) +
                          sum(sapply(seq(1, len.y, by=1), log_state_dens, states=states, theta=theta_prop)) -
                          log_prior_density_sigmasmall(theta_curr[1], beta.small, gamma.small) -
                          sum(sapply(seq(1, len.y, by=1), log_state_dens, states=states, theta=theta_curr))))
  }
  if(runif(1)<pacc[1]){
    theta_curr[1]<-theta_prop[1]
  } else {
    theta_prop[1]<-theta_curr[1]
  }

  ## sigma eta large
  theta_prop[2]<-runif(1, theta_curr[2]-prop_lim[2], theta_curr[2]+prop_lim[2])
  if(theta_prop[2]>theta_curr[2]){
    pacc[2]<-min(1, exp(log_prior_density_sigmalarge(theta_prop[2], beta.large, gamma.large) +
                          sum(sapply(seq(1, len.y, by=1), log_state_dens, states=states, theta=theta_prop)) -
                          log_prior_density_sigmalarge(theta_curr[2], beta.large, gamma.large) -
                          sum(sapply(seq(1, len.y, by=1), log_state_dens, states=states, theta=theta_curr))))
  }
  if(runif(1)<pacc[2]){
    theta_curr[2]<-theta_prop[2]
  } else {
    theta_prop[2]<-theta_curr[2]
  }

  ## Gibbs step for observation variance
  theta_curr[4]<-rinvgamma(1, shape=beta.sigma+len.y/2, rate=gamma.sigma + 0.5*sum((y-states)^2))
  pacc[4]<-1

  theta_curr
} ## fixed input and output

##################### log observation and state density functions
log_obs_dens<-function(y, states, t, theta){

  dnorm(y[t], states[t], sqrt(theta[4]), log=TRUE)

} ## fixed input and output
log_state_dens<-function(states, t, theta){
  if(t==1){
    log(theta[3]*dnorm(states[t], 1, sqrt(theta[1])) + (1-theta[3])*dnorm(states[t], 1, sqrt(theta[2])))
  } else {
    log(theta[3]*dnorm(states[t], states[t-1], sqrt(theta[1])) + (1-theta[3])*dnorm(states[t], states[t-1], sqrt(theta[2])))
  }
} ## fixed input and output

###################### proposal distribution for the states
end_grid_cell_prop<-function(const, end.prop.var){
  const + rnorm(1, 0, sqrt(end.prop.var))
}
start_grid_cell_prop<-function(const, end.prop.var){
  const - rnorm(1, 0, sqrt(end.prop.var))
}
mid_grid_cell_prop<-function(lower_quant, upper_quant){
  runif(1, lower_quant, upper_quant)
}

end_grid_cell_dens<-function(prop.state, const, end.prop.var){
  #log density
  dnorm(prop.state-const, 0, sqrt(end.prop.var), log=TRUE)
}
start_grid_cell_dens<-function(prop.state, const, end.prop.var){
  dnorm(const-prop.state, 0, sqrt(end.prop.var), log=TRUE)
}
mid_grid_cell_dens<-function(prop.state, lower_quant, upper_quant){
  dunif(prop.state, lower_quant, upper_quant, log=TRUE)
}


################## run algorithm
len.y<-length(y)

bl<-gen_blocks(l, len.y)
blocks=bl$blocks
blocksize=bl$blocksize


states_curr<-matrix(rep(0, len.y*Nits), nrow=Nits)
theta_curr<-matrix(rep(0, length(theta.init)*Nits), nrow=Nits)
states_prop<-matrix(rep(0, len.y*Nits), nrow=Nits)
theta_prop<-matrix(rep(0, length(theta.init)*Nits), nrow=Nits)
pacc_states<-matrix(rep(0, dim(blocks)[1]*Nits), nrow=Nits)
states_curr[1,]<-states.init
theta_curr[1,]<-theta.init

for(i in 2:Nits){
  #update states
  run_it<-PMPMH(states_curr[i-1,], y, blocks, N, q1, qN, sigma, theta_curr[i-1,], delta.e, thresh)
  states_curr[i,]<-run_it$states
  states_prop[i,]<-run_it$`proposed states`
  pacc_states[i,]<-run_it$states_pacc

  #update theta
  theta_curr[i,]<-theta_update(states_curr[i,], theta_curr[i-1,])
}



