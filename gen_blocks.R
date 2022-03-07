gen_blocks<-function(l, len.y){
  split=seq(1, len.y, by=l-1)
  while(split[length(split)]<len.y){
    split=c(split, len.y)
  }

  blocks<-matrix(seq(split[1], split[2], by=1), ncol=split[2]-split[1]+1, byrow=FALSE)
  for(k in 2:(length(split)-2)){
    blocks<-rbind(blocks, seq(split[k], split[k+1]))
  }

  if(max(blocks)<len.y){
    temp=seq(split[length(split)-1], len.y, by=1)
    while(length(temp)<l){
      temp=c(temp, NA)
    }
    blocks=rbind(blocks,temp)
  }

  find.block.size<-function(row){
    length(which(!is.na(row)))
  }
  blocksize<-as.numeric(apply(blocks, MARGIN=1, find.block.size))

  return(list("blocks"=blocks, "blocksize"=blocksize))

}
