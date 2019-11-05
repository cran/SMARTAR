#getcmat() is the function to build the contrast design matrix

getcmat=function(control,nats){
        CON=control; G=nats                            #total number of stragegies
        Cmat=matrix(rep(0,G*(G-1)),ncol=G)             #contrast matrix
        for (j in 1:G){
             if (j==CON) {Cmat[ ,j]= 1} else
             if (j<CON)  {Cmat[j,j]=-1} else
             if (j>CON)  {Cmat[j-1,j]=-1}
             }
        return(Cmat)
}

