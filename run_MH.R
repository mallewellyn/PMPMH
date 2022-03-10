source("gen_HMM.R")
source("gen_proposal.R")


gen_pacc<-function(block, prop, states_curr, theta_curr, state_lag){
  states_prop<-states_curr
  states_prop[block]<-prop$states_prop[block]

  if(block[1]==1){
    t.eval.states<-c(block, max(block)+state_lag)
  }
  if(max(block)==len.y){
    t.eval.states<-c(min(block)-state_lag, block)
  }
  if(min(block)!=1 && max(block)!=len.y){
    t.eval.states<-c(min(block)-state_lag, block, max(block)+state_lag)
  }

  pacc<-min(1, exp(prop$log.prob.curr + sum(sapply(t.eval.states, log_state_dens, states=states_prop, theta=theta_curr)) +
                     sum(sapply(block, log_obs_dens, y=y, states=states_prop, theta=theta_curr)) -
                     (prop$log.prob.new + sum(sapply(t.eval.states, log_state_dens, states=states_curr, theta=theta_curr)) +
                        sum(sapply(block, log_obs_dens, y=y, states=states_curr, theta=theta_curr)))))

  if(is.na(pacc)){
    pacc=0
  }
  pacc
}

PMPMH<-function(states_curr, y, blocks, N, q1, qN, var.infl, theta_curr, delta.e, thresh, approach, state_lag){

  len.y=length(y)

  pacc_states<-numeric(dim(blocks)[1])
  states_prop<-numeric(len.y)

  if(approach!=3){
  grid=gen_HMM(seq(1, len.y, by=1), y, N, q1, qN, var.infl, states_curr, theta_curr, delta.e, thresh, state_lag)
  }

  for(l in 1:dim(blocks)[1]){
    if(approach==3){
      grid=gen_HMM(na.omit(blocks[l,]), y, N, q1, qN, var.infl, states_curr, theta_curr, delta.e, thresh, state_lag)

    }
    proposal<-sample_HMM(states_curr, y, na.omit(blocks[l,]), grid, end.prop.var, state_lag)
    pacc_states[l]<-gen_pacc(na.omit(blocks[l,]), proposal, states_curr, theta_curr, state_lag)
    states_prop[na.omit(blocks[l,])]=proposal$states_prop[na.omit(blocks[l,])]

    if(runif(1)<pacc_states[l]){
      states_curr[na.omit(blocks[l,])]<-proposal$states_prop[na.omit(blocks[l,])]
    }
  }

  return(list("states"=states_curr, "proposed states"=states_prop, "states_pacc"=pacc_states))
}

