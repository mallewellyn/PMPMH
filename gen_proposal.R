
sample_HMM<-function(states.curr, blocksize, index.curr, trProbs, obProbs, N.adapt, quants, end.prop.var){
  index.samp<-numeric(blocksize)
  log.prob.curr<-0
  log.prob.new<-0
  grid.traj.curr<-numeric(blocksize)
  grid.traj.new<-numeric(blocksize)
  states.prop<-numeric(blocksize)

  filt_temp<-trProbs[[1]]*obProbs[[1]]
  filt_temp<-norm_rows(filt_temp)

  filt=list(filt_temp)

  for(i in 2:blocksize){
    filt_temp<-apply(filt[[i-1]]*trProbs[[i]], MARGIN=2, sum)*obProbs[[i]]
    filt_temp<-norm_rows(filt_temp)
    filt[[i]]=filt_temp
  }

  ## sampling
  if(length(trProbs)==blocksize){
    smooth=filt[[blocksize]]
  } else {
    smooth<-filt[[blocksize]]*trProbs[[blocksize+1]]
  }
  smooth<-norm_rows(smooth)
  index.samp[blocksize]<-which((runif(1)-(cumsum(smooth)/cumsum(smooth)[length(smooth)])<0))[1]

  for(i in (blocksize-1):1){
    smooth<-filt[[i]]*trProbs[[i+1]][,index.samp[i+1]]
    smooth<-norm_rows(smooth)
    index.samp[i]<-which((runif(1)-(cumsum(smooth)/cumsum(smooth)[length(smooth)])<0))[1]
  }

  #proposal and current probability calculations
  log.prob.curr<-log.prob.curr+log(trProbs[[1]][index.curr[1]]) + log(obProbs[[1]][index.curr[1]])
  for(i in 2:blocksize){
    log.prob.curr<-log.prob.curr + log(trProbs[[i]][index.curr[i-1], index.curr[i]]) +
      log(obProbs[[i]][index.curr[i]])
  }
  if(length(trProbs)!=blocksize){
    log.prob.curr<-log.prob.curr+log(trProbs[[blocksize+1]][index.curr[blocksize]])
  }
  log.prob.new<-log.prob.new+log(trProbs[[1]][index.samp[1]]) + log(obProbs[[1]][index.samp[1]])
  for(i in 2:blocksize){
    log.prob.new<-log.prob.new + log(trProbs[[i]][index.samp[i-1], index.samp[i]]) +
      log(obProbs[[i]][index.samp[i]])
  }
  if(length(trProbs)!=blocksize){
  log.prob.new<-log.prob.new+log(trProbs[[blocksize+1]][index.samp[blocksize]])
  }


  for(k in 1:blocksize){
    if(index.samp[k]==N.adapt[k]){
      states.prop[k]<-end_grid_cell_prop(quants[[k]][N.adapt[k]-1], end.prop.var)
      log.prob.new<-log.prob.new + end_grid_cell_dens(states.prop[k], quants[[k]][N.adapt[k]-1], end.prop.var)
    }
    if(index.curr[k]==N.adapt[k]){
      log.prob.curr<-log.prob.curr + end_grid_cell_dens(states.curr[k], quants[[k]][N.adapt[k]-1], end.prop.var)
    }

    if(index.samp[k]==1){
      states.prop[k]<-start_grid_cell_prop(quants[[k]][1], end.prop.var)
      log.prob.new<-log.prob.new + start_grid_cell_dens(states.prop[k], quants[[k]][1], end.prop.var)
    }
    if(index.curr[k]==1){
      log.prob.curr<-log.prob.curr + start_grid_cell_dens(states.curr[k], quants[[k]][1], end.prop.var)
    }

    if(index.samp[k]!=1 & index.samp[k]!=N.adapt[k]){
      if(index.samp[k]!=N.adapt[k]){
        states.prop[k]<-mid_grid_cell_prop(quants[[k]][index.samp[k]-1], quants[[k]][index.samp[k]])
      log.prob.new<-log.prob.new+mid_grid_cell_dens(states.prop[k], quants[[k]][index.samp[k]-1], quants[[k]][index.samp[k]])
      }
    }
    if(index.curr[k]!=1 & index.curr[k]!=N.adapt[k]){
      if(index.samp[k]!=N.adapt[k]){
        log.prob.curr<-log.prob.curr + mid_grid_cell_dens(states.curr[k], quants[[k]][index.curr[k]-1], quants[[k]][index.curr[k]])
      }
    }

  }
  return(list("index.samp"=index.samp, "states.prop"=states.prop, "log.prob.new"=log.prob.new, "log.prob.curr"=log.prob.curr))
}

gen_proposal<-function(block, grid, states_curr){
  Midpoints=grid$Midpoints
  Quants=grid$Quants
  N.adapt=grid$`Number bins`
  Transition=grid$Transition
  Observation=grid$Observation

  index.curr<-numeric(len.y)
  for(k in 1:len.y){
    index.curr[k]<-find_curr_index(states_curr[k], Quants[[k]])
  }


  #first block
  if(min(na.omit(block))==1){
  trEnd<-Transition[[block[length(na.omit(block))]+1]][,index.curr[block[length(na.omit(block))]+1]]
  trAll<-Transition[block]
  trAll<-c(trAll, list(trEnd))
  }
  if(max(na.omit(block))==len.y){
    trStart<-Transition[[block[1]]][index.curr[block[1]-1],]
    trAll<-c(list(trStart), Transition[as.numeric(na.omit(block))[-1]])
  }

  if(min(na.omit(block))!=1 && max(na.omit(block))!=len.y){
    trStart<-Transition[[block[1]]][index.curr[block[1]-1],]
    trAll<-c(list(trStart), Transition[block[-1]])
    trEnd<-Transition[[block[length(na.omit(block))]+1]][,index.curr[block[length(na.omit(block))]+1]]
    trAll<-c(trAll, list(trEnd))
  }

  propose<-sample_HMM(states_curr[block], length(na.omit(block)), index.curr[block], trAll,
                      Observation[block], N.adapt[block], Quants[blocks], end.prop.var)
  propose
}

