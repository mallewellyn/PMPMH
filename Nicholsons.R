library(invgamma)

source("run_MH.R")
source("gen_blocks.R")


##################### user defined ##########################
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
beta.e=1/sigma2.e
beta.epsilon=1/sigma2.epsilon
e=rgamma(n, beta.e, rate=beta.e)
epsilon=rgamma(n, beta.epsilon, rate=beta.epsilon)
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

set.seed(NULL)


Nits=5e4
l=3
N=50
q1=0.001
qN=1-q1
var.infl=2
delta.e=2
thresh=0.001
end.prop.var=10
approach=3
state_lag=tau+1
states.init=c(rep(0, state_lag-1), y[-(1:(state_lag-1))])
theta.init<-c(0.5, 62.5, 1/3, 1/3, 1.5) #(delta, P, beta_epsilon, beta_e, phi)
e.init<-rep(3, length(y))
epsilon.init<-rep(3, length(y))

##################### user defined log observation and state density functions
log_obs_dens<-function(y, states, t, theta){
  density=dpois(y[t], theta[5]*(S[t]+states[t]), log=TRUE)

} ## fixed input and output

log_state_dens<-function(states, t, theta){
  if(t < state_lag){
    dens=log(1)
    } else {
      if(t==state_lag){
        dens=dpois(states[t], theta[2]*N0*e_curr[t]*exp(-1), log=TRUE)
      } else {
        N=S[t-state_lag]+states[t-state_lag]
        dens=dpois(states[t], theta[2]*N*e_curr[t]*exp(-N/N0), log=TRUE)
      }
    }
  dens

} ## fixed input and output

log_S_dens<-function(S, states, t, theta){
  if(t==1){
    dens=dbinom(S[t], N0, exp(-theta[1]*epsilon_curr[t]), log=TRUE)
  } else {
  if(t<=state_lag-1){
    dens=dbinom(S[t], S[t-1], exp(-theta[1]*epsilon_curr[t]), log=TRUE)
  } else {
    dens=dbinom(S[t], S[t-1] + states[t-1], exp(-theta[1]*epsilon_curr[t]), log=TRUE)
  }
  }
  dens
}

################ user-defined quantile function
#calculate quantiles
quantile_func<-function(perc, t, states_curr, y, theta, var.infl){
  if(t<state_lag){
    quants=0
  } else {
   if(states_curr[t]==0){
    quants=qpois(perc, var.infl) - var.infl + 1
  } else {
    quants=pmax(qpois(perc, states_curr[t]*var.infl)-states_curr[t]*var.infl + states_curr[t], 0)
  }
  }
  unique(round(quants/2)*2)
}
#calculate midpoints
midpoint_func<-function(q, delta.e, t){
  if(t<state_lag){
    N.adapt=1
    mp=0
    len.cell=1
  } else {
    N.adapt<-length(q)+1

    mp<-numeric(N.adapt)
    len.cell<-numeric(N.adapt)


    mp[1]<-pmax(as.integer(q[1]-delta.e/2), 0)
    mp[N.adapt]<-as.integer(q[length(q)] + delta.e/2)
    len.cell[c(1, N.adapt)]<-rep(delta.e, 2)


    mp[-c(1, N.adapt)]=(q[-1]+q[-(length(q))])/2
    len.cell[-c(1, N.adapt)]=c(q[-1]-q[-length(q)])
  }
  return(list("N.adapt"=N.adapt, "mp"=mp, "len.cell"=len.cell))
}

#calculate midpoint integration
midpoint_int_func_system<-function(mpoints, bin.len, theta, t, state_lag){
  reformat<-function(mpminus, mp) log_state_dens(c(numeric(t-2), mpminus, mp), t, theta)
  vreformat<-Vectorize(reformat, vectorize.args = c("mpminus", "mp"))

  if(t<state_lag){
    dens=1
    } else {
      if(t==state_lag){
        dens=outer(0, mpoints[[t]], vreformat) + log(bin.len[[t]])
      } else {
     dens=outer(mpoints[[t-state_lag]], mpoints[[t]], vreformat)
      dens=dens+bin.len[[t-state_lag]]
      dens=t(t(dens) + bin.len[[t]])
      }
    }


  dens

}
midpoint_int_func_obs<-function(y, mpoints, bin.len, theta, t){
  reformat<-function(mp, t) log_obs_dens(y, c(numeric(t-1), mp), t, theta)
  sapply(mpoints[[t]], reformat, t=t) + log(bin.len[[t]])

}

########################## scheme for updating theta
delta.a=0.0081
delta.b=0.01
P.a=50
P.b=1
beta.epsilon.a=100
beta.epsilon.b=100
beta.e.a=100
beta.e.b=100
phi.a=0.01
phi.b=0.01


prop.lim=c(0.1, 10, 1, 1/3, 1)
prop.epsilon=2
prop.e=4

log_prior_density_delta<-function(val, delta.a, delta.b){
  dgamma(val, delta.a, rate=delta.b, log=TRUE)
}
log_prior_density_beta_epsilon<-function(val, beta.epsilon.a, beta.epsilon.b){
  dinvgamma(val, shape=beta.epsilon.a, rate=beta.epsilon.b, log=TRUE)
}
log_prior_density_beta_e<-function(val, beta.e.a, beta.e.b){
  dinvgamma(val, shape=beta.e.a, rate=beta.e.b, log=TRUE)
}

theta_update<-function(states, theta){
  N.curr=S+states
  theta_curr<-theta
  theta_prop<-theta
  pacc<-numeric(length(theta))
  e_prop=e_curr
  e_pacc=numeric(length(y))
  epsilon_prop=epsilon_curr
  epsilon_pacc=numeric(length(y))

  ## first delta
  theta_prop[1]=runif(1, theta_curr[1]-prop.lim[1], theta_curr[1]+prop.lim[1])
  if(theta_prop[1]>0){
    pacc[1]=min(1, exp(log_prior_density_delta(theta_prop[1], delta.a, delta.b) +
                          sum(sapply(seq(1, length(y), by=1), log_S_dens, S=S, states=states, theta=theta_prop)) -
                                (log_prior_density_delta(theta_curr[1], delta.a, delta.b) +
                                   sum(sapply(seq(1, n, by=1), log_S_dens, S=S, states=states, theta=theta_curr)))))
  }
  if(is.na(pacc[1])){
    pacc[1]=0
  }
  if(runif(1)<pacc[1]){
      theta_curr[1]=theta_prop[1]
   } else {
      theta_prop[1]=theta_curr[1]
  }



  ## P
  N=states+S
  theta_prop[2]=rgamma(1, P.a + sum(states[state_lag:length(y)]),
                        rate=P.b+ N0*exp(-1)*e_curr[state_lag] + sum(N.curr[1:(length(y)-state_lag)]*exp(-N.curr[1:(length(y)-state_lag)]/N0)*e_curr[(state_lag+1):length(y)]))
  theta_curr[2]=theta_prop[2]
  pacc[2]=1

  #beta_epsilon
  theta_prop[3]=runif(1, theta_curr[3]-prop.lim[3], theta_curr[3]+prop.lim[3])
  if(theta_prop[3]>0){
    pacc[3]=min(1, exp(log_prior_density_beta_epsilon(theta_prop[3], beta.epsilon.a, beta.epsilon.b) +
                                 sum(sapply(epsilon_curr, dgamma, shape=theta_prop[3], rate=theta_prop[3], log=TRUE)) -
                                (log_prior_density_beta_epsilon(theta_curr[3], beta.epsilon.a, beta.epsilon.b) +
                                   sum(sapply(epsilon_curr, dgamma, shape=theta_curr[3], rate=theta_curr[3], log=TRUE)))))
  }
  if(is.na(pacc[3])){
    pacc[3]=0
  }
    if(runif(1)<pacc[3]){
      theta_curr[3]=theta_prop[3]
    } else {
      theta_prop[3]=theta_curr[3]
    }

  #beta_e
  theta_prop[4]=runif(1, theta_curr[4]-prop.lim[4], theta_curr[4] + prop.lim[4])
  if(theta_prop[4]>0){
    pacc[4]=min(1, exp(log_prior_density_beta_e(theta_prop[4], beta.e.a, beta.e.b) +
                         sum(sapply(e_curr, dgamma, shape=theta_prop[4], rate=theta_prop[4], log=TRUE)) - (
                           log_prior_density_beta_e(theta_curr[4], beta.e.a, beta.e.b) +
                             sum(sapply(e_curr, dgamma, shape=theta_curr[4], rate=theta_curr[4], log=TRUE)))))
  }
  if(is.na(pacc[4])){
    pacc[4]=0
  }

  if(runif(1)<pacc[4]){
    theta_curr[4]=theta_prop[4]
  } else {
    theta_prop[4]=theta_curr[4]
  }

  #epsilon
  epsilon_prop[1]=runif(1, epsilon_curr[1]-prop.epsilon, epsilon_curr[1]+prop.epsilon)
  e=epsilon_curr
  epsilon_curr=epsilon_prop
  num=dgamma(epsilon_prop[1], theta_curr[3], rate=theta_curr[3], log=TRUE) +
    log_S_dens(S, states, 1, theta_curr)
  epsilon_curr=e
  if(epsilon_prop[1]>0){
    epsilon_pacc[1]=min(1, exp(num - (dgamma(epsilon_curr[1], theta_curr[3], rate=theta_curr[3], log=TRUE) +
                                    log_S_dens(S, states, 1, theta_curr))))
  }
  if(is.na(epsilon_pacc[1])){
    epsilon_pacc[1]=0
  }
    if(runif(1)<epsilon_pacc[1]){
      epsilon_curr[1]=epsilon_prop[1]
    } else {
      epsilon_prop[1]=epsilon_curr[1]
    }

  for(t in 2:length(y)){

    epsilon_prop[t]=runif(1, epsilon_curr[t]-prop.epsilon, epsilon_curr[t]+prop.epsilon)
    e=epsilon_curr
    epsilon_curr=epsilon_prop
    num=dgamma(epsilon_prop[t], theta_curr[3], rate=theta_curr[3], log=TRUE) +
      log_S_dens(S, states, t, theta_curr)
    epsilon_curr=e
    if(epsilon_prop[t]>0){
      epsilon_pacc[t]=min(1, exp( num - (
                                      dgamma(epsilon_curr[t], theta_curr[3], rate=theta_curr[3], log=TRUE) +
                                        log_S_dens(S, states, t, theta_curr))))

    }
    if(is.na(epsilon_pacc[t])){
      epsilon_pacc[t]=0
    }
    if(runif(1)<epsilon_pacc[t]){
        epsilon_curr[t]=epsilon_prop[t]
      } else {
        epsilon_prop[t]=epsilon_curr[t]
      }
  }

  ##e
  e_prop[state_lag]=rgamma(1, theta_curr[4] + states[state_lag], rate=theta_curr[4] + theta_curr[2]*N0*exp(-1))
  e_pacc[state_lag]=1
  e_curr[state_lag]=e_prop[state_lag]

  for(t in (state_lag+1):length(y)){
    e_prop[t]=rgamma(1, theta_curr[4]+states[t], rate=theta_curr[4]+theta_curr[2]*N.curr[t-state_lag]*exp(-N.curr[t-state_lag]/N0))
    e_pacc[t]=1
    e_curr[t]=e_prop[t]
  }


  theta_prop[5]=rgamma(1, phi.a + sum(y), rate=phi.b + sum(states+S))
  theta_curr[5]=theta_prop[5]
  pacc[5]=1

  return(list(theta_curr, epsilon_curr, e_curr))
}


###################### proposal distribution for the states
end_grid_cell_prop<-function(const, end.prop.var){
  const+rpois(1, end.prop.var)
}
mid_grid_cell_prop<-function(lower_quant, upper_quant){
  sample(seq(lower_quant, upper_quant, by=1), size=1)
}
start_grid_cell_prop<-function(const, end.prop.var){
  sample(seq(0, const, by=1), size=1)
}

end_grid_cell_dens<-function(prop.state, const, end.prop.var){
  dpois(prop.state-const, end.prop.var, log=TRUE)
}
mid_grid_cell_dens<-function(prop.state, lower_quant, upper_quant){
  log(1/(upper_quant-lower_quant+1))
}
start_grid_cell_dens<-function(prop.state, const, end.prop.var){
  log(1/(const+1))
}

################## run algorithm
len.y<-length(y)
bl<-gen_blocks(l, len.y, state_lag)
blocks=bl$blocks


states_curr<-matrix(rep(0, len.y*Nits), nrow=Nits)
theta_curr<-matrix(rep(0, length(theta.init)*Nits), nrow=Nits)
e_samps<-matrix(rep(0, len.y*Nits), nrow=Nits)
epsilon_samps<-matrix(rep(0, len.y*Nits), nrow=Nits)
states_prop<-matrix(rep(0, len.y*Nits), nrow=Nits)
pacc_states<-matrix(rep(0, dim(blocks)[1]*Nits), nrow=Nits)
states_curr[1,]<-states.init
theta_curr[1,]<-theta.init
e_samps[1,]<-e.init
epsilon_samps[1,]<-epsilon.init

for(i in 2:Nits){
  epsilon_curr=epsilon_samps[i-1,]
  e_curr=e_samps[i-1,]

  #update states
  run_it<-PMPMH(states_curr[i-1,], y, blocks, N, q1, qN, var.infl, theta_curr[i-1,], delta.e, thresh, approach, state_lag)
  states_curr[i,]<-run_it$states
  states_prop[i,]<-run_it$`proposed states`
  pacc_states[i,]<-run_it$states_pacc

  #update parameters
  param_temp<-theta_update(states_curr[i,], theta_curr[i-1,])
  theta_curr[i,]<-param_temp[[1]]
  epsilon_samps[i,]=param_temp[[2]]
  e_samps[i,]=param_temp[[3]]
  print(i)
}
