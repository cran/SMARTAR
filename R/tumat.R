#tumat() get true vector of strategy values
#inout sequence information matrix

tumat<-function(sim,family="gaussian"){
  Smat<-sim
  FA<-family
  if (is.null(Smat$O1)) {Base<-0} else {Base<-1}

  Nstage<-nstage(data=Smat)
  Dmat<-atsscan(data=Smat)
  G<-nrow(Dmat)

  Uvec<-NULL
  if (Nstage==1 && Base==0) {
      for (g in 1:G) {ATS<-as.numeric(Dmat[g,2])
                      ATS[which(is.na(ATS))]<-0
                      S1<-Smat[which(Smat$A1==ATS[1]),]
                      u<-S1$MEAN
                      Uvec<-c(Uvec,u)
                      }} else
  if (Nstage==1 && Base==1) {
      for (g in 1:G) {ATS<-as.numeric(Dmat[g,2:3])
                      ATS[which(is.na(ATS))]<-0
                      S0<-Smat[which(Smat$O1==0 & Smat$A1==ATS[1]),]
                      S1<-Smat[which(Smat$O1==1 & Smat$A1==ATS[2]),]
                      u0<-S0$MEAN; if (length(u0)==0L) {u0<-0}
                      u1<-S1$MEAN; if (length(u1)==0L) {u1<-0}
                      p0<-S0$P1; if (length(p0)==0L) {p0<-0}
                      p1<-S1$P1; if (length(p1)==0L) {p1<-0}
                      u<-sum(p0*u0,p1*u1)
                      Uvec<-c(Uvec,u)
                      }} else
  if (Nstage==2 && Base==0) {
      for (g in 1:G){
        ATS<-as.numeric(Dmat[g,2:4])
        ATS[which(is.na(ATS))]<-0
        S0<-Smat[which(Smat$A1==ATS[1] &
                        Smat$O2==0 & Smat$A2==ATS[2]),] #seq (d1,O2=0,d20)
        S1<-Smat[which(Smat$A1==ATS[1] &
                        Smat$O2==1 & Smat$A2==ATS[3]),] #seq (d1,O2=1.d21)
        u0<-S0$MEAN; if (length(u0)==0L) {u0<-0}
        u1<-S0$MEAN; if (length(u1)==0L) {u1<-0}
        p0<-S0$P2; if (length(p0)==0L) {p0<-0}
        p1<-S1$P2; if (length(p1)==0L) {p1<-0}
        u<-sum(p0*u0,p1*u1)
        Uvec<-c(Uvec,u)
        }} else
  if (Nstage==2 && Base==1) {
    for (g in 1:G){
      ATS<-as.numeric(Dmat[g,2:8])
      ATS[which(is.na(ATS))]<-0
      S00<-Smat[which(Smat$O1==0 & Smat$A1==ATS[1] &
                       Smat$O2==0 & Smat$A2==ATS[3]),]
      S01<-Smat[which(Smat$O1==0 & Smat$A1==ATS[1] &
                       Smat$O2==1 & Smat$A2==ATS[4]),]
      S10<-Smat[which(Smat$O1==1 & Smat$A1==ATS[2] &
                       Smat$O2==0 & Smat$A2==ATS[5]),]
      S11<-Smat[which(Smat$O1==1 & Smat$A1==ATS[2] &
                       Smat$O2==1 & Smat$A2==ATS[6]),]
      u00<-S00$MEAN; if (length(u00)==0L) {u00<-0}
      u01<-S01$MEAN; if (length(u01)==0L) {u01<-0}
      u10<-S10$MEAN; if (length(u10)==0L) {u10<-0}
      u11<-S11$MEAN; if (length(u11)==0L) {u11<-0}
      p0<-S01$P1;  if (length(p0)==0L)  {p0<-0}
      p1<-S11$P1;  if (length(p1)==0L)  {p1<-0}
      p00<-S00$P2; if (length(p00)==0L) {p00<-0}
      p01<-S01$P2; if (length(p01)==0L) {p01<-0}
      p10<-S00$P2; if (length(p10)==0L) {p10<-0}
      p11<-S01$P2; if (length(p11)==0L) {p11<-0}

      u<-sum(p0*p00*u00,p0*p01*u01,p1*p10*u10,p1*p11*u11)
      Uvec<-c(Uvec,u)
      }} else
  if (Nstage==3 && Base==0) {
    for (g in 1:G){
      ATS<-as.numeric(Dmat[g,2:8])
      ATS[which(is.na(ATS))]<-0
      S00<-Smat[which(Smat$A1==ATS[1] & Smat$O2==0 &
                       Smat$A2==ATS[2] & Smat$O3==0 & Smat$A3==ATS[4]),]
      S01<-Smat[which(Smat$A1==ATS[1] & Smat$O2==0 &
                       Smat$A2==ATS[2] & Smat$O3==1 & Smat$A3==ATS[5]),]
      S10<-Smat[which(Smat$A1==ATS[1] & Smat$O2==1 &
                       Smat$A2==ATS[3] & Smat$O3==0 & Smat$A3==ATS[6]),]
      S11<-Smat[which(Smat$A1==ATS[1] & Smat$O2==1 &
                       Smat$A2==ATS[3] & Smat$O3==1 & Smat$A3==ATS[7]),]
      u00<-S00$MEAN; if (length(u00)==0L) {u00<-0}
      u01<-S01$MEAN; if (length(u01)==0L) {u01<-0}
      u10<-S10$MEAN; if (length(u10)==0L) {u10<-0}
      u11<-S11$MEAN; if (length(u11)==0L) {u11<-0}
      p0<-S01$P2;  if (length(p0)==0L)  {p0<-0}
      p1<-S11$P2;  if (length(p1)==0L)  {p1<-0}
      p00<-S00$P3; if (length(p00)==0L) {p00<-0}
      p01<-S01$P3; if (length(p01)==0L) {p01<-0}
      p10<-S00$P3; if (length(p10)==0L) {p10<-0}
      p11<-S01$P3; if (length(p11)==0L) {p11<-0}
      u<-sum(p0*p00*u00,p0*p01*u01,p1*p10*u10,p1*p11*u11)
      Uvec<-c(Uvec,u)
      }} 
  ##else
  ##if (Nstage==3 && Base==1) {
  ##  for (g in 1:G){
  ##    ATS<-as.numeric(Dmat[g,2:15])
  ##    ATS[which(is.na(ATS))]<-0
  ##    S000<-Smat[which(Smat$O1==0 & Smat$A1==ATS[1] &
  ##                      Smat$O2==0 & Smat$A2==ATS[3] &
  ##                      Smat$O3==0 & Smat$A3==ATS[7]),]
  ##    S001<-Smat[which(Smat$O1==0 & Smat$A1==ATS[1] &
  ##                      Smat$O2==0 & Smat$A2==ATS[3] &
  ##                      Smat$O3==1 & Smat$A3==ATS[8]),]
  ##    S010<-Smat[which(Smat$O1==0 & Smat$A1==ATS[1] &
  ##                      Smat$O2==1 & Smat$A2==ATS[4] &
  ##                      Smat$O3==0 & Smat$A3==ATS[9]),]
  ##    S011<-Smat[which(Smat$O1==0 & Smat$A1==ATS[1] &
  ##                      Smat$O2==1 & Smat$A2==ATS[4] &
  ##                      Smat$O3==1 & Smat$A3==ATS[10]),]
  ##    S100<-Smat[which(Smat$O1==1 & Smat$A1==ATS[2] &
  ##                      Smat$O2==0 & Smat$A2==ATS[5] &
  ##                      Smat$O3==0 & Smat$A3==ATS[11]),]
  ##    S101<-Smat[which(Smat$O1==1 & Smat$A1==ATS[2] &
  ##                      Smat$O2==0 & Smat$A2==ATS[5] &
  ##                      Smat$O3==1 & Smat$A3==ATS[12]),]
  ##    S110<-Smat[which(Smat$O1==1 & Smat$A1==ATS[2] &
  ##                      Smat$O2==1 & Smat$A2==ATS[6] &
  ##                      Smat$O3==0 & Smat$A3==ATS[13]),]
  ##    S111<-Smat[which(Smat$O1==1 & Smat$A1==ATS[2] &
  ##                      Smat$O2==1 & Smat$A2==ATS[6] &
  ##                      Smat$O3==1 & Smat$A3==ATS[14]),]
  ##    u000<-S000$MEAN; if (length(u000)==0L) {u000<-0}
  ##    u001<-S001$MEAN; if (length(u001)==0L) {u001<-0}
  ##    u010<-S010$MEAN; if (length(u010)==0L) {u010<-0}
  ##    u011<-S011$MEAN; if (length(u011)==0L) {u011<-0}
  ##    u100<-S100$MEAN; if (length(u100)==0L) {u100<-0}
  ##    u101<-S101$MEAN; if (length(u101)==0L) {u101<-0}
  ##    u110<-S110$MEAN; if (length(u110)==0L) {u110<-0}
  ##    u111<-S111$MEAN; if (length(u111)==0L) {u111<-0}
  ##    p0<-S001$P1;  if (length(p0)==0L)  {p0<-0}
  ##    p1<-S101$P1;  if (length(p1)==0L)  {p1<-0}
##
  ##    p00<-S001$P2; if (length(p00)==0L) {p00<-0}
  ##    p01<-S011$P2; if (length(p01)==0L) {p01<-0}
  ##    p10<-S101$P2; if (length(p10)==0L) {p10<-0}
  ##    p11<-S111$P2; if (length(p11)==0L) {p11<-0}
##
  ##    p000<-S000$P3; if (length(p000)==0L) {p000<-0}
  ##    p001<-S001$P3; if (length(p001)==0L) {p001<-0}
  ##    p010<-S010$P3; if (length(p010)==0L) {p010<-0}
  ##    p011<-S011$P3; if (length(p011)==0L) {p011<-0}
  ##    p100<-S100$P3; if (length(p100)==0L) {p100<-0}
  ##    p101<-S101$P3; if (length(p101)==0L) {p101<-0}
  ##    p110<-S110$P3; if (length(p110)==0L) {p110<-0}
  ##    p111<-S111$P3; if (length(p111)==0L) {p111<-0}
##
  ##    u<-sum(p0*p00*p000*u000,
  ##          p0*p00*p001*u001,
  ##          p0*p01*p010*u010,
  ##          p0*p01*p011*u011,
  ##          p1*p10*p100*u100,
  ##          p1*p10*p101*u101,
  ##          p1*p11*p110*u110,
  ##          p1*p11*p111*u111)
  ##    Uvec<-c(Uvec,u)
  ##    }}
      Umat<-matrix(Uvec,ncol=1)
      return(Umat)
}
