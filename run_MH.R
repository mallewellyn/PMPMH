source("/home/s1521656/OneDrive/PMPMH paper/Code/PMPMH/R/Github_code/gen_HMM.R")
source("/home/s1521656/OneDrive/PMPMH paper/Code/PMPMH/R/Github_code/gen_proposal.R")

find_curr_index<-function(states, Quants){
  inters<-matrix(c(-Inf, rep(Quants, each=2), Inf), nrow=2, byrow=FALSE)
  index=which(inters[2,which(inters[1,]<states)]>=states)
  index
}

gen_pacc<-function(block, prop, states_curr, theta_curr){
  block<-na.omit(block)
  states_prop<-states_curr
  states_prop[block]<-prop$states.prop

  if(min(block)==1){
    t.eval.states<-c(block, max(block)+1)
  }
  if(max(block)==len.y){
    t.eval.states<-c(min(block)-1, block)
  }
  if(min(block)!=1 && max(block)!=len.y){
    t.eval.states<-c(min(block)-1, block, max(block)+1)
  }

  pacc<-min(1, exp(prop$log.prob.curr + sum(sapply(t.eval.states, log_state_dens, states=states_prop, theta=theta_curr)) +
                     sum(sapply(block, log_obs_dens, y=y, states=states_prop, theta=theta_curr)) -
                     (prop$log.prob.new + sum(sapply(t.eval.states, log_state_dens, states=states_curr, theta=theta_curr)) +
                        sum(sapply(block, log_obs_dens, y=y, states=states_curr, theta=theta_curr)))))

  pacc
}

PMPMH<-function(states_curr, y, blocks, N, q1, qN, sigma, theta_curr, delta.e, thresh){


  pacc_states<-numeric(dim(blocks)[1])
  states_prop<-numeric(length(y))

  grid=gen_HMM(y, N, q1, qN, sigma, states_curr, theta_curr, delta.e, thresh)

  for(l in 1:dim(blocks)[1]){
    proposal<-gen_proposal(blocks[l,], grid, states_curr)
    pacc_states[l]<-gen_pacc(blocks[l,], proposal, states_curr, theta_curr)
    states_prop[na.omit(blocks[l,])]<-proposal$states.prop

    if(runif(1)<pacc_states[l]){
      states_curr[na.omit(blocks[l,])]<-proposal$states.prop
    }
  }

  return(list("states"=states_curr, "proposed states"=states_prop, "states_pacc"=pacc_states))
}

