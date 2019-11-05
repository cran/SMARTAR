#getag() gives tag of all pairwise comparison, which presented as "ats.I vs. ats.J" (I,J defined in atsscan)

getag=function(nats){
      tag=NULL
      for (i in 1:nats){
           for (j in 1:nats){
                if (i!=j) {tag=rbind(tag,c(i,j))}
                }}
      colnames(tag)=c("I","J")
      return(tag)
}

