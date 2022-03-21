gen_blocks<-function(l, len.y, state_lag){
  x<-matrix(seq(state_lag, state_lag-1+state_lag*l, by=1), ncol=l, nrow=state_lag)
  i=1
  while(max(x)<len.y){
    x<-rbind(x, matrix(seq(x[((i-1)*state_lag+1),l], x[((i-1)*state_lag+1),l]+state_lag*l-1, by=1), ncol=l, nrow=state_lag))
    i=i+1
  }

  for(i in dim(x)[1]:1){
    if(length(which(x[i,]>len.y))==l){
      x=x[-i,]
    } else {
      if(length(which(x[i,]>len.y))>0){
        x[i,which(x[i,]>len.y)]=NA
      }
    }

  }

  blocks=x
  find.block.size<-function(row){
    length(which(!is.na(row)))
  }
  blocksize<-as.numeric(apply(blocks, MARGIN=1, find.block.size))

  return(list("blocks"=blocks, "blocksize"=blocksize))

}
