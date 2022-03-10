source("run_MH.R")
source("gen_blocks.R")
library(invgamma)
library(truncnorm)


##################### user defined ##########################
set.seed(66)
n=600
x=numeric(n)
y=numeric(n)
mu.prior=1
sigma.epsilon=1
jump.prop=0.90
sigma.eta.large=700
sigma.eta.small=1

u=runif(n)
if(u[1]<jump.prop){
  x[1]=rnorm(1, mu.prior, sqrt(sigma.eta.small))
} else {
  x[1]=rnorm(1, mu.prior, sqrt(sigma.eta.large))
}
y[1]=rnorm(1, x[1], sqrt(sigma.epsilon))
for(i in 2:n){
  if(u[i]<jump.prop){
    x[i]=rnorm(1, x[i-1], sqrt(sigma.eta.small))
  } else {
    x[i]=rnorm(1, x[i-1], sqrt(sigma.eta.large))
  }
  y[i]=rnorm(1, x[i], sqrt(sigma.epsilon))
}


Nits=1e5
l=4
N=10
q1=0.3
qN=1-q1
var.infl=1
delta.e=0.25
thresh=0
end.prop.var=5
states.init=rep(1, length(y))
theta.init<-c(5, 100, 0.5, 5)
approach=3
state_lag=1
##################### user defined log observation and state density functions
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

################ user-defined quantile function
#calculate quantiles
quantile_func<-function(perc, t, states_curr, y, theta, var.infl){
  quants=qnorm(perc, states_curr[t], sqrt(var.infl))
  unique(quants)
}
#calculate midpoints
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

#calculate midpoint integration
midpoint_int_func_system<-function(mpoints, bin.len, theta, t, state_lag){
  reformat<-function(mpminus, mp) log_state_dens(c(numeric(t-2), mpminus, mp), t, theta)
  vreformat<-Vectorize(reformat, vectorize.args = c("mpminus", "mp"))

  if(t==1){
    dens=sapply(mpoints[[1]], log_state_dens, t=1, theta=theta) + log(bin.len[[1]])
  } else {
    dens=outer(mpoints[[t-state_lag]], mpoints[[t]], vreformat)
    dens=dens+bin.len[[t-state_lag]]
    dens=t(t(dens) + bin.len[[t]])
  }


  dens

}
midpoint_int_func_obs<-function(y, mpoints, bin.len, theta, t){
  reformat<-function(mp, t) log_obs_dens(y, c(numeric(t-1), mp), t, theta)

  sapply(mpoints[[t]], reformat, t=t) + log(bin.len[[t]])

}

########################## scheme for updating theta
prop_lim<-c(1.6, 80, 0.3, 1)
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
    pacc[3]<-min(1, exp(sum(sapply(seq(1, len.y, by=1), log_state_dens, states=states, theta=theta_prop)) -
                           sum(sapply(seq(1, len.y, by=1), log_state_dens, states=states, theta=theta_curr))))

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
  if(theta_prop[2]>theta_curr[1]){
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


###################### proposal distribution for the states
end_grid_cell_prop<-function(const, end.prop.var){
  rtruncnorm(1, a=const, b=Inf, mean=const, sqrt(end.prop.var))
}
start_grid_cell_prop<-function(const, end.prop.var){
  rtruncnorm(1, a=-Inf, b=const, mean=const, sqrt(end.prop.var))
}
mid_grid_cell_prop<-function(lower_quant, upper_quant){
  runif(1, lower_quant, upper_quant)
}

end_grid_cell_dens<-function(prop.state, const, end.prop.var){
  log(dtruncnorm(prop.state, a=const, b=Inf, mean=const, sqrt(end.prop.var)))
}
start_grid_cell_dens<-function(prop.state, const, end.prop.var){
  log(dtruncnorm(prop.state, a=-Inf, b=const, mean=const, sqrt(end.prop.var)))
}
mid_grid_cell_dens<-function(prop.state, lower_quant, upper_quant){
  dunif(prop.state, lower_quant, upper_quant, log=TRUE)
}


################## run algorithm
len.y<-length(y)

bl<-gen_blocks(l, len.y)
blocks=bl$blocks


states_curr<-matrix(rep(0, len.y*Nits), nrow=Nits)
theta_curr<-matrix(rep(0, length(theta.init)*Nits), nrow=Nits)
states_prop<-matrix(rep(0, len.y*Nits), nrow=Nits)
theta_prop<-matrix(rep(0, length(theta.init)*Nits), nrow=Nits)
pacc_states<-matrix(rep(0, dim(blocks)[1]*Nits), nrow=Nits)
states_curr[1,]<-states.init
theta_curr[1,]<-theta.init

for(i in 2:Nits){
    #update states
    run_it<-PMPMH(states_curr[i-1,], y, blocks, N, q1, qN, var.infl, theta_curr[i-1,], delta.e, thresh, approach, state_lag)

    states_curr[i,]<-run_it$states
    states_prop[i,]<-run_it$`proposed states`
    pacc_states[i,]<-run_it$states_pacc

    #update parameters
    #theta_curr[i,]<-theta_update(states_curr[i,], theta_curr[i-1,])
    theta_curr[i,]<-c(sigma.eta.small, sigma.eta.large, jump.prop, sigma.epsilon)
    print(i)
  }

t=1
plot.ts(states_curr[,t][1:i])
lines(states10_1[,t][1:i], col="blue")
lines(states_prop[,t][1:i], col="red")
lines(rep(x[t], i), col="orange")
mean(pacc_states[1:i,])

var(states10_1[,t][4000:7000])
var(states_curr[,t][4000:7000])
var(states7_1[,t][4000:7000])

t=10
plot.ts(states_10_1[,t+1][1:7750])
lines(states_10_2[,t][1:7750], col="red")
lines(rep(x[t], 7750), col="orange")
lines(rep(states_curr[1,t], 7750), col="red")


write.csv(states_curr, file="/home/s1521656/OneDrive/PMPMH paper/checks/states_N10_q103_var1_l4_fixedparams.csv")
write.csv(states_prop, file="/home/s1521656/OneDrive/PMPMH paper/checks/statesprop_N10_q103_var1_l4_fixedparams.csv")

states_10_1<-read.csv(file="/home/s1521656/OneDrive/PMPMH paper/checks/states_N10_q10001_var5_l4_fixedparams.csv")
