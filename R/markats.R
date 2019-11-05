#markats() assign index ind=1 for those subjects following strategy of interest, otherwise ind=0

markats=function(data,ats){
        D=as.data.frame(data); ATS=ats
        if (is.null(D$O1)) {Base=0} else {Base=1}
        Nstage=nstage(data=D)
        N=nrow(D)
        IND=rep(0,N)
        D=data.frame(D,IND)
        if (Nstage==1 & Base==0) {d1=ATS[1]                                                   #Base=0 => d1 is a scalar
                                  D$IND[which(D$A1==d1)]=1
                                  }
        if (Nstage==1 & Base==1) {Hvec=sort(unique(D$O1))
                                  for (h in Hvec) {
                                       dp1=which(Hvec==h)
                                       d1=ATS[0+dp1]                                            #DM for Hvec[h]
                                       D$IND[which(D$O1==h & D$A1==d1)]=1
                                       }}
        if (Nstage==2 & Base==0) {x=0
                                  d1=ATS[1]                                                   #Base=0 => d1 is a scalar
                                  Jvec=sort(unique(D$O2[which(D$A1==d1)]))
                                  for (j in Jvec){
                                       x=x+1
                                       d2=ATS[1+x]
                                       #dp2=which(Jvec==j)
                                       #d2=ATS[1+dp2]                                                       #1st is d1, 2:3 is d2
                                       D$IND[which(D$A1==d1 & D$O2==j & D$A2==d2)]=1
                                       }}
        if (Nstage==2 & Base==1) {x=y=0
                                  Hvec=sort(unique(D$O1))
                                  for (h in Hvec){
                                       x=x+1
                                       d1=ATS[x]
                                       #dp1=which(Hvec==h)
                                       #d1=ATS[dp1]                                                       #DM for Hvec[h]
                                       Jvec=sort(unique(D$O2[which(D$O1==h & D$A1==d1)]))
                                       for (j in Jvec){
                                            y=y+1
                                            d2=ATS[2+y]
                                            #dp2=which(Jvec==j)
                                            #d2=ATS[2+dp2]
                                            D$IND[which(D$O1==h & D$A1==d1 & D$O2==j & D$A2==d2)]=1
                                            }}}
        if (Nstage==3 & Base==0) {x=y=0
                                  d1=ATS[1]                                                               #Base=0 => d1 is a scalar
                                  Jvec=sort(unique(D$O2[which(D$A1==d1)]))
                                  for (j in Jvec){
                                       x=x+1
                                       d2=ATS[1+x]
                                       #dp2=which(Jvec==j)
                                       #d2=ATS[1+dp2]                                                       #1st is d1, 2:3 is d2
                                       Lvec=sort(unique(D$O3[which(D$A1==d1 & D$O2==j & D$A2==d2)]))
                                       for (l in Lvec){
                                            y=y+1
                                            d3=ATS[3+y]
                                            #dp3=which(Lvec==l)
                                            #d3=ATS[3+dp3]
                                            D$IND[which(D$A1==d1 & D$O2==j & D$A2==d2 & D$O3==l & D$A3==d3)]=1
                                            }}}
        if (Nstage==3 & Base==1) {x=y=z=0
                                  Hvec=sort(unique(D$O1))
                                  for (h in Hvec){
                                       x=x+1
                                       d1=ATS[x]
                                       #dp1=which(Hvec==h)
                                       #d1=ATS[dp1]                                                       #DM for Hvec[h]
                                       Jvec=sort(unique(D$O2[which(D$O1==h & D$A1==d1)]))
                                       for (j in Jvec){
                                            y=y+1
                                            d2=ATS[2+y]
                                            #dp2=which(Jvec==j)
                                            #d2=ATS[2+dp2]
                                            Lvec=sort(unique(D$O3[which(D$O1==h & D$A1==d1 & D$O2==j & D$A2==d2)]))
                                            for (l in Lvec){
                                                 z=z+1
                                                 d3=ATS[6+z]
                                                 #dp3=which(Lvec==l)
                                                 #d3=ATS[6+dp3]
                                                 D$IND[which(D$O1==h & D$A1==d1 & D$O2==j & D$A2==d2 & D$O3==l & D$A3==d3)]=1
                                                 }}}}
        return(D$IND)
        #return(D)
}


