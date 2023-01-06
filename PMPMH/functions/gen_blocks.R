genBlocks<-function(l, leny, stateLag){
  x<-matrix(seq(stateLag, stateLag-1+stateLag*l, by=1), ncol=l, nrow=stateLag)
  i=1
  while(max(x)<leny){
    x<-rbind(x, matrix(seq(x[((i-1)*stateLag+1),l], x[((i-1)*stateLag+1),l]+stateLag*l-1, by=1), ncol=l, nrow=stateLag))
    i=i+1
  }

  for(i in dim(x)[1]:1){
    if(length(which(x[i,]>leny))==l){
      x=x[-i,]
    } else {
      if(length(which(x[i,]>leny))>0){
        x[i,which(x[i,]>leny)]=NA
      }
    }

  }

  blocks=x
  findBlockSize<-function(row){
    length(which(!is.na(row)))
  }
  blockSize<-as.numeric(apply(blocks, MARGIN=1, findBlockSize))

  return(list("blocks"=blocks, "blockSize"=blockSize))

}
