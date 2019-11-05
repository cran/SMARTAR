#ecov() estimted covariance between two DTRs (variance when two DTR identical)
#' @importFrom stats var

ecov=function(data,ats1=c(0),ats2=c(0),family="normal",common=T){
     D=as.data.frame(data); ATS1=ats1; ATS2=ats2; FA=family; C=common
     if (is.null(D$O1)) {Base=0} else {Base=1}
     Nstage=nstage(data=D)
     N=nrow(D)
     ATS1[which(is.na(ATS1))]=0
     ATS2[which(is.na(ATS2))]=0

     I=overlap(ats1=ATS1,ats2=ATS2,nstage=Nstage,baseline=Base)

     if (Nstage==1 & Base==0) {DJ=D[which(D$A1==ATS1[1]),]
                               DK=D[which(D$A1==ATS2[1]),]

                               nJ=length(DK$Y)

                               uJ=mean(DJ$Y); if (is.na(uJ)) {uJ=0}                        # within-sequence true mean (EY)
                               uK=mean(DK$Y); if (is.na(uK)) {uK=0}

                               vJ=var(DJ$Y);  if (is.na(vJ)) {vJ=0}                        #within-sequence variance for cont outcome
                               vK=var(DK$Y);  if (is.na(vK)) {vK=0}
                               bvJ=uJ*(1-uJ); if (is.na(bvJ)) {bvJ=0}                      #within-sequence var est for bin outcome
                               bvK=uK*(1-uK); if (is.na(bvK)) {bvK=0}
                               comvar=cv(data=D);  if (is.na(comvar)) {comvar=0}

                               if (FA=="binomial") {vc=I[1]*bvJ/nJ}
                                       else if (C) {vc=I[1]*comvar/nJ}
                                              else {vc=I[1]*vJ/nJ}
                              } else
     if (Nstage==1 & Base==1) {DJ=D[which( (D$O1==0 & D$A1==ATS1[1])|(D$O1==1 & D$A1==ATS1[2]) ),]
                               DK=D[which( (D$O1==0 & D$A1==ATS2[1])|(D$O1==1 & D$A1==ATS2[2]) ),]

                               n0J=length(DJ$Y[which(DJ$O1==0)]); if (n0J==0|is.na(n0J)) {n0J=1}               #sequence size of R=0 under dJ
                               n1J=length(DJ$Y[which(DJ$O1==1)]); if (n1J==0|is.na(n1J)) {n1J=1}               #sequence size of R=1 under dJ

                               P=length(D$Y[which(D$O1==1)])/length(D$Y)

                               u0J=mean(DJ$Y[which(DJ$O1==0)]); if (is.na(u0J)) {u0J=0}                        # within-sequence true mean (EY)
                               u1J=mean(DJ$Y[which(DJ$O1==1)]); if (is.na(u1J)) {u1J=0}
                               u0K=mean(DK$Y[which(DK$O1==0)]); if (is.na(u0K)) {u0K=0}
                               u1K=mean(DK$Y[which(DK$O1==1)]); if (is.na(u1K)) {u1K=0}
                               v0J=var(DJ$Y[which(DJ$O1==0)]);  if (is.na(v0J)) {v0J=0}                        #within-sequence variance for cont outcome
                               v1J=var(DJ$Y[which(DJ$O1==1)]);  if (is.na(v1J)) {v1J=0}
                               v0K=var(DK$Y[which(DK$O1==0)]);  if (is.na(v0K)) {v0K=0}
                               v1K=var(DK$Y[which(DK$O1==1)]);  if (is.na(v1K)) {v1K=0}
                               bv0J=u0J*(1-u0J);  if (is.na(bv0J)) {bv0J=0}                                    #within-sequence var est for bin outcome
                               bv1J=u1J*(1-u1J);  if (is.na(bv1J)) {bv1J=0}
                               bv0K=u0K*(1-u0K);  if (is.na(bv0K)) {bv0K=0}
                               bv1K=u1K*(1-u1K);  if (is.na(bv1K)) {bv1K=0}
                               comvar=cv(data=D);  if (is.na(comvar)) {comvar=0}

                               if (FA=="binomial") {B1=1*P*(1-P)*(u1J-u0J)*(u1K-u0K)/N          #MLE for binary outcome
                                                    B2=I[1]*(1-P)*(1-P)*bv0J/n0J
                                                    B3=I[2]*P*P*bv1J/n1J
                                                    vc=B1+B2+B3}
                                       else if (C) {B1=1*P*(1-P)*(u1J-u0J)*(u1K-u0K)/N          #MLE for cont outcome with com assumption
                                                    B2=I[1]*(1-P)*(1-P)*comvar/n0J
                                                    B3=I[2]*P*P*comvar/n1J
                                                    vc=B1+B2+B3}
                                              else {B1=1*P*(1-P)*(u1J-u0J)*(u1K-u0K)/N          #MLE for cont outcome without com assumption
                                                    B2=I[1]*(1-P)*(1-P)*v0J/n0J
                                                    B3=I[2]*P*P*v1J/n1J
                                                    vc=B1+B2+B3}
                              } else
     if (Nstage==2 & Base==0) {DJ=D[which((D$A1==ATS1[1] & D$O2==0 & D$A2==ATS1[2])|
                                          (D$A1==ATS1[1] & D$O2==1 & D$A2==ATS1[3])),] #so far no O1
                               DK=D[which((D$A1==ATS2[1] & D$O2==0 & D$A2==ATS2[2])|
                                          (D$A1==ATS2[1] & D$O2==1 & D$A2==ATS2[3])),]

                               nj=length(D$Y[which(D$A1==ATS1[1])])                                            #sample size of A1=d1
                               n0J=length(DJ$Y[which(DJ$O2==0)]); if (n0J==0|is.na(n0J)) {n0J=1}               #sequence size of R=0 under dJ
                               n1J=length(DJ$Y[which(DJ$O2==1)]); if (n1J==0|is.na(n1J)) {n1J=1}               #sequence size of R=1 under dJ

                               P=length(D$Y[which(D$A1==ATS1[1]&D$O2==1)])/length(D$Y[which(D$A1==ATS1[1])])   #estimate respose rate P=P(R=1|A1=d1J)

                               u0J=mean(DJ$Y[which(DJ$O2==0)]); if (is.na(u0J)) {u0J=0}                        # within-sequence true mean (EY)
                               u1J=mean(DJ$Y[which(DJ$O2==1)]); if (is.na(u1J)) {u1J=0}
                               u0K=mean(DK$Y[which(DK$O2==0)]); if (is.na(u0K)) {u0K=0}
                               u1K=mean(DK$Y[which(DK$O2==1)]); if (is.na(u1K)) {u1K=0}
                               v0J=var(DJ$Y[which(DJ$O2==0)]);  if (is.na(v0J)) {v0J=0}                        #within-sequence variance for cont outcome
                               v1J=var(DJ$Y[which(DJ$O2==1)]);  if (is.na(v1J)) {v1J=0}
                               v0K=var(DK$Y[which(DK$O2==0)]);  if (is.na(v0K)) {v0K=0}
                               v1K=var(DK$Y[which(DK$O2==1)]);  if (is.na(v1K)) {v1K=0}
                               bv0J=u0J*(1-u0J);  if (is.na(bv0J)) {bv0J=0}                                    #within-sequence var est for bin outcome
                               bv1J=u1J*(1-u1J);  if (is.na(bv1J)) {bv1J=0}
                               bv0K=u0K*(1-u0K);  if (is.na(bv0K)) {bv0K=0}
                               bv1K=u1K*(1-u1K);  if (is.na(bv1K)) {bv1K=0}
                               comvar=cv(data=D);  if (is.na(comvar)) {comvar=0}                 #common within-sequence variance

                               if (FA=="binomial") {B1=I[1]*P*(1-P)*(u1J-u0J)*(u1K-u0K)/nj          #MLE for binary outcome
                                                    B2=I[2]*(1-P)*(1-P)*bv0J/n0J
                                                    B3=I[3]*P*P*bv1J/n1J
                                                    vc=B1+B2+B3}
                                       else if (C) {B1=I[1]*P*(1-P)*(u1J-u0J)*(u1K-u0K)/nj          #MLE for cont outcome with com assumption
                                                    B2=I[2]*(1-P)*(1-P)*comvar/n0J
                                                    B3=I[3]*P*P*comvar/n1J
                                                    vc=B1+B2+B3}
                                              else {B1=I[1]*P*(1-P)*(u1J-u0J)*(u1K-u0K)/nj          #MLE for cont outcome without com assumption
                                                    B2=I[2]*(1-P)*(1-P)*v0J/n0J
                                                    B3=I[3]*P*P*v1J/n1J
                                                    vc=B1+B2+B3}
                               } else
     if (Nstage==2 & Base==1) {DJ=D[which((D$O1==0 & D$A1==ATS1[1] & D$O2==0 & D$A2==ATS1[3])|
                                          (D$O1==0 & D$A1==ATS1[1] & D$O2==1 & D$A2==ATS1[4])|
                                          (D$O1==1 & D$A1==ATS1[2] & D$O2==0 & D$A2==ATS1[5])|
                                          (D$O1==1 & D$A1==ATS1[2] & D$O2==1 & D$A2==ATS1[6])),]     #so far no O1
                               DK=D[which((D$O1==0 & D$A1==ATS2[1] & D$O2==0 & D$A2==ATS2[3])|
                                          (D$O1==0 & D$A1==ATS2[1] & D$O2==1 & D$A2==ATS2[4])|
                                          (D$O1==1 & D$A1==ATS2[2] & D$O2==0 & D$A2==ATS2[5])|
                                          (D$O1==1 & D$A1==ATS2[2] & D$O2==1 & D$A2==ATS2[6])),]

                               n0J=length(D$Y[which(D$O1==0 & D$A1==ATS1[1])]); if (n0J==0|is.na(n0J)) {n0J=1}                                                           #sample size of A1=d1
                               n1J=length(D$Y[which(D$O1==1 & D$A1==ATS1[2])]); if (n1J==0|is.na(n1J)) {n1J=1}

                               n00J=length(DJ$Y[which(DJ$O1==0 & DJ$O2==0)]); if (n00J==0|is.na(n00J)) {n00J=1}               #sequence size of R=0 under dJ
                               n01J=length(DJ$Y[which(DJ$O1==0 & DJ$O2==1)]); if (n01J==0|is.na(n01J)) {n01J=1}               #sequence size of R=1 under dJ
                               n10J=length(DJ$Y[which(DJ$O1==1 & DJ$O2==0)]); if (n10J==0|is.na(n10J)) {n10J=1}               #sequence size of R=0 under dJ
                               n11J=length(DJ$Y[which(DJ$O1==1 & DJ$O2==1)]); if (n11J==0|is.na(n11J)) {n11J=1}

                               P0=length(D$Y[which(D$O1==0)])/length(D$Y)   #estimate respose rate P=P(R=1|A1=d1J)
                               P1=1-P0

                               P00J=length(D$Y[which(D$O1==0 & D$A1==ATS1[1] & D$O2==0)])/
                                    length(D$Y[which(D$O1==0 & D$A1==ATS1[1])])
                               P01J=1-P00J
                               P10J=length(D$Y[which(D$O1==1 & D$A1==ATS1[2] & D$O2==0)])/
                                    length(D$Y[which(D$O1==1 & D$A1==ATS1[2])])
                               P11J=1-P10J

                               P00K=length(D$Y[which(D$O1==0 & D$A1==ATS2[1] & D$O2==0)])/
                                    length(D$Y[which(D$O1==0 & D$A1==ATS2[1])])
                               P01K=1-P00K
                               P10K=length(D$Y[which(D$O1==1 & D$A1==ATS2[2] & D$O2==0)])/
                                    length(D$Y[which(D$O1==1 & D$A1==ATS2[2])])
                               P11K=1-P10K

                               u00J=mean(DJ$Y[which(DJ$O1==0 & DJ$O2==0)]); if (is.na(u00J)) {u00J=0}        # within-sequence true mean (EY)
                               u01J=mean(DJ$Y[which(DJ$O1==0 & DJ$O2==1)]); if (is.na(u01J)) {u01J=0}
                               u10J=mean(DJ$Y[which(DJ$O1==1 & DJ$O2==0)]); if (is.na(u10J)) {u10J=0}
                               u11J=mean(DJ$Y[which(DJ$O1==1 & DJ$O2==1)]); if (is.na(u11J)) {u11J=0}

                               u00K=mean(DJ$Y[which(DK$O1==0 & DK$O2==0)]); if (is.na(u00K)) {u00K=0}
                               u01K=mean(DJ$Y[which(DK$O1==0 & DK$O2==1)]); if (is.na(u01K)) {u01K=0}
                               u10K=mean(DJ$Y[which(DK$O1==1 & DK$O2==0)]); if (is.na(u10K)) {u10K=0}
                               u11K=mean(DJ$Y[which(DK$O1==1 & DK$O2==1)]); if (is.na(u11K)) {u11K=0}

                               v00J=var(DJ$Y[which(DJ$O1==0 & DJ$O2==0)]);  if (is.na(v00J)) {v00J=0}         #within-sequence variance for cont outcome
                               v01J=var(DJ$Y[which(DJ$O1==0 & DJ$O2==1)]);  if (is.na(v01J)) {v01J=0}
                               v10J=var(DJ$Y[which(DJ$O1==1 & DJ$O2==0)]);  if (is.na(v10J)) {v10J=0}
                               v11J=var(DJ$Y[which(DJ$O1==1 & DJ$O2==1)]);  if (is.na(v11J)) {v11J=0}

                               bv00J=u00J*(1-u00J);  if (is.na(bv00J)) {bv00J=0}                                   #within-sequence var est for bin outcome
                               bv01J=u01J*(1-u01J);  if (is.na(bv01J)) {bv01J=0}
                               bv10J=u10J*(1-u10J);  if (is.na(bv10J)) {bv10J=0}
                               bv11J=u11J*(1-u11J);  if (is.na(bv11J)) {bv11J=0}
                               comvar=cv(data=D);  if (is.na(comvar)) {comvar=0}

                               A0000= (P0*(1-P0)*P00J*P00K/N+I[1]*P0*P0*P00J*(1-P00J)/n0J)*u00J*u00K
                               A0001= (P0*(1-P0)*P00J*P01K/N-I[1]*P0*P0*P00J*(1-P00J)/n0J)*u00J*u01K
                               A0010=-(P0*(1-P0)*P00J*P10K/N)*u00J*u10K
                               A0011=-(P0*(1-P0)*P00J*P11K/N)*u00J*u11K

                               A0100= (P0*(1-P0)*P01J*P00K/N-I[1]*P0*P0*P01J*(1-P01J)/n0J)*u01J*u00K
                               A0101= (P0*(1-P0)*P01J*P01K/N+I[1]*P0*P0*P01J*(1-P01J)/n0J)*u01J*u01K
                               A0110=-(P0*(1-P0)*P01J*P10K/N)*u01J*u10K
                               A0111=-(P0*(1-P0)*P01J*P11K/N)*u01J*u11K

                               A1000=-(P1*(1-P1)*P10J*P00K/N)*u10J*u10K
                               A1001=-(P1*(1-P1)*P10J*P01K/N)*u10J*u11K
                               A1010= (P1*(1-P1)*P10J*P10K/N+I[1]*P1*P1*P10J*(1-P10J)/n1J)*u10J*u10K
                               A1011= (P1*(1-P1)*P10J*P11K/N-I[1]*P1*P1*P10J*(1-P10J)/n1J)*u10J*u11K

                               A1100=-(P1*(1-P1)*P10J*P00K/N)*u10J*u10K
                               A1101=-(P1*(1-P1)*P10J*P01K/N)*u10J*u11K
                               A1110= (P1*(1-P1)*P11J*P10K/N+I[1]*P1*P1*P11J*(1-P11J)/n1J)*u11J*u10K
                               A1111= (P1*(1-P1)*P10J*P11K/N-I[1]*P1*P1*P11J*(1-P11J)/n1J)*u11J*u11K
                               A=sum(A0000,A0001,A0010,A0011,
                                     A0100,A0101,A0110,A0111,
                                     A1000,A1001,A1010,A1011,
                                     A1100,A1101,A1110,A1111)

                               if (FA=="binomial"){if (I[1]==1 & I[3]==1) {B0000=P0*P0*P00J*P00K*bv00J/n00J}  else {B0000=0}
                                                   if (I[1]==1 & I[4]==1) {B0101=P0*P0*P01J*P01K*bv01J/n01J}  else {B0101=0}
                                                   if (I[2]==1 & I[5]==1) {B1010=P1*P1*P10J*P10K*bv10J/n10J}  else {B1010=0}
                                                   if (I[2]==1 & I[6]==1) {B1111=P1*P1*P11J*P11K*bv11J/n11J}  else {B1111=0}
                                                   }
                                      else if (C) {if (I[1]==1 & I[3]==1) {B0000=P0*P0*P00J*P00K*comvar/n00J} else {B0000=0}
                                                   if (I[1]==1 & I[4]==1) {B0101=P0*P0*P01J*P01K*comvar/n01J} else {B0101=0}
                                                   if (I[2]==1 & I[5]==1) {B1010=P1*P1*P10J*P10K*comvar/n10J} else {B1010=0}
                                                   if (I[2]==1 & I[6]==1) {B1111=P1*P1*P11J*P11K*comvar/n11J} else {B1111=0}
                                                   }
                                             else {if (I[1]==1 & I[3]==1) {B0000=P0*P0*P00J*P00K*v00J/n00J}   else {B0000=0}
                                                   if (I[1]==1 & I[4]==1) {B0101=P0*P0*P01J*P01K*v01J/n01J}   else {B0101=0}
                                                   if (I[2]==1 & I[5]==1) {B1010=P1*P1*P10J*P10K*v10J/n10J}   else {B1010=0}
                                                   if (I[2]==1 & I[6]==1) {B1111=P1*P1*P11J*P11K*v11J/n11J}   else {B1111=0}
                                                   }
                               B=sum(B0000,B0101,B1010,B1111)
                               vc=A+B
                               } else
     if (Nstage==3 & Base==0) {DJ=D[which((D$A1==ATS1[1] & D$O2==0 & D$A2==ATS1[2] & D$O3==0 & D$A3==ATS1[4])|
                                          (D$A1==ATS1[1] & D$O2==0 & D$A2==ATS1[2] & D$O3==1 & D$A3==ATS1[5])|
                                          (D$A1==ATS1[1] & D$O2==1 & D$A2==ATS1[3] & D$O3==0 & D$A3==ATS1[6])|
                                          (D$A1==ATS1[1] & D$O2==1 & D$A2==ATS1[3] & D$O3==1 & D$A3==ATS1[7])),]
                               DK=D[which((D$A1==ATS2[1] & D$O2==0 & D$A2==ATS2[2] & D$O3==0 & D$A3==ATS2[4])|
                                          (D$A1==ATS2[1] & D$O2==0 & D$A2==ATS2[2] & D$O3==1 & D$A3==ATS2[5])|
                                          (D$A1==ATS2[1] & D$O2==1 & D$A2==ATS2[3] & D$O3==0 & D$A3==ATS2[6])|
                                          (D$A1==ATS2[1] & D$O2==1 & D$A2==ATS2[3] & D$O3==1 & D$A3==ATS2[7])),]

                               nJ=length(D$Y[which(D$A1==ATS1[1])]); if (nJ==0|is.na(nJ)) {nJ=1}  #sampe size to est Phat(O2)

                               n0J=length(D$Y[which(D$A1==ATS1[1] & D$O2==0 & D$A1==ATS1[2])]); if (n0J==0|is.na(n0J)) {n0J=1} #sample size est Phat(O3=0)                                                         #sample size of A1=d1
                               n1J=length(D$Y[which(D$A1==ATS1[1] & D$O2==1 & D$A1==ATS1[3])]); if (n1J==0|is.na(n1J)) {n1J=1} #sample size est Phat(O3==1)

                               n00J=length(DJ$Y[which(DJ$O2==0 & DJ$O3==0)]); if (n00J==0|is.na(n00J)) {n00J=1}  #seq of (O1,O2)=(0,0) under dJ
                               n01J=length(DJ$Y[which(DJ$O2==0 & DJ$O3==1)]); if (n01J==0|is.na(n01J)) {n01J=1}  #seq of (O1,O2)=(0,1) under dJ
                               n10J=length(DJ$Y[which(DJ$O2==1 & DJ$O3==0)]); if (n10J==0|is.na(n10J)) {n10J=1}  #seq of (O1,O2)=(1,0) under dJ
                               n11J=length(DJ$Y[which(DJ$O2==1 & DJ$O3==1)]); if (n11J==0|is.na(n11J)) {n11J=1}  #seq of (O1,O2)=(1,1) under dJ

                               P0J=length(D$Y[which(D$A1==ATS1[1] & D$O2==0)])/
                                   length(D$Y[which(D$A1==ATS1[1])])
                               P1J=1-P0J

                               P00J=length(D$Y[which(D$A1==ATS1[1] & D$O2==0 & D$A2==ATS1[2] & D$O3==0)])/
                                    length(D$Y[which(D$A1==ATS1[1] & D$O2==0 & D$A2==ATS1[2])])
                               P01J=1-P00J
                               P10J=length(D$Y[which(D$A1==ATS1[1] & D$O2==1 & D$A2==ATS1[3] & D$O3==0)])/
                                    length(D$Y[which(D$A1==ATS1[1] & D$O2==1 & D$A2==ATS1[3])])
                               P11J=1-P10J

                               P00K=length(D$Y[which(D$A1==ATS2[1] & D$O2==0 & D$A2==ATS2[2] & D$O3==0)])/
                                    length(D$Y[which(D$A1==ATS2[1] & D$O2==0 & D$A2==ATS2[2])])
                               P01K=1-P00K
                               P10K=length(D$Y[which(D$A1==ATS2[1] & D$O2==1 & D$A2==ATS2[3] & D$O3==0)])/
                                    length(D$Y[which(D$A1==ATS1[1] & D$O2==1 & D$A2==ATS2[3])])
                               P11K=1-P10K

                               u00J=mean(DJ$Y[which(DJ$O2==0 & DJ$O3==0)]); if (is.na(u00J)) {u00J=0}        # within-sequence true mean (EY)
                               u01J=mean(DJ$Y[which(DJ$O2==0 & DJ$O3==1)]); if (is.na(u01J)) {u01J=0}
                               u10J=mean(DJ$Y[which(DJ$O2==1 & DJ$O3==0)]); if (is.na(u10J)) {u10J=0}
                               u11J=mean(DJ$Y[which(DJ$O2==1 & DJ$O3==1)]); if (is.na(u11J)) {u11J=0}

                               u00K=mean(DJ$Y[which(DK$O2==0 & DK$O3==0)]); if (is.na(u00K)) {u00K=0}
                               u01K=mean(DJ$Y[which(DK$O2==0 & DK$O3==1)]); if (is.na(u01K)) {u01K=0}
                               u10K=mean(DJ$Y[which(DK$O2==1 & DK$O3==0)]); if (is.na(u10K)) {u10K=0}
                               u11K=mean(DJ$Y[which(DK$O2==1 & DK$O3==1)]); if (is.na(u11K)) {u11K=0}

                               v00J=var(DJ$Y[which(DJ$O2==0 & DJ$O3==0)]);  if (is.na(v00J)) {v00J=0}         #within-sequence variance for cont outcome
                               v01J=var(DJ$Y[which(DJ$O2==0 & DJ$O3==1)]);  if (is.na(v01J)) {v01J=0}
                               v10J=var(DJ$Y[which(DJ$O2==1 & DJ$O3==0)]);  if (is.na(v10J)) {v10J=0}
                               v11J=var(DJ$Y[which(DJ$O2==1 & DJ$O3==1)]);  if (is.na(v11J)) {v11J=0}

                               bv00J=u00J*(1-u00J);  if (is.na(bv00J)) {bv00J=0}                              #within-sequence var est for bin outcome
                               bv01J=u01J*(1-u01J);  if (is.na(bv01J)) {bv01J=0}
                               bv10J=u10J*(1-u10J);  if (is.na(bv10J)) {bv10J=0}
                               bv11J=u11J*(1-u11J);  if (is.na(bv11J)) {bv11J=0}
                               comvar=cv(data=D);    if (is.na(comvar)) {comvar=0}                            #pooled variance

                               A0000= (I[1]*P0J*(1-P0J)*P00J*P00K/nJ+I[2]*P0J*P0J*P00J*(1-P00J)/n0J)*u00J*u00K
                               A0001= (I[1]*P0J*(1-P0J)*P00J*P01K/nJ-I[2]*P0J*P0J*P00J*(1-P00J)/n0J)*u00J*u01K
                               A0010=-(I[1]*P0J*(1-P0J)*P00J*P10K/nJ)*u00J*u10K
                               A0011=-(I[1]*P0J*(1-P0J)*P00J*P11K/nJ)*u00J*u11K

                               A0100= (I[1]*P0J*(1-P0J)*P01J*P00K/nJ-I[2]*P0J*P0J*P01J*(1-P01J)/n0J)*u01J*u00K
                               A0101= (I[1]*P0J*(1-P0J)*P01J*P01K/nJ+I[2]*P0J*P0J*P01J*(1-P01J)/n0J)*u01J*u01K
                               A0110=-(I[1]*P0J*(1-P0J)*P01J*P10K/nJ)*u01J*u10K
                               A0111=-(I[1]*P0J*(1-P0J)*P01J*P11K/nJ)*u01J*u11K

                               A1000=-(I[1]*P1J*(1-P1J)*P10J*P00K/nJ)*u10J*u10K
                               A1001=-(I[1]*P1J*(1-P1J)*P10J*P01K/nJ)*u10J*u11K
                               A1010= (I[1]*P1J*(1-P1J)*P10J*P10K/nJ+I[3]*P1J*P1J*P10J*(1-P10J)/n1J)*u10J*u10K
                               A1011= (I[1]*P1J*(1-P1J)*P10J*P11K/nJ-I[3]*P1J*P1J*P10J*(1-P10J)/n1J)*u10J*u11K

                               A1100=-(I[1]*P1J*(1-P1J)*P10J*P00K/nJ)*u10J*u10K
                               A1101=-(I[1]*P1J*(1-P1J)*P10J*P01K/nJ)*u10J*u11K
                               A1110= (I[1]*P1J*(1-P1J)*P11J*P10K/nJ+I[3]*P1J*P1J*P11J*(1-P11J)/n1J)*u11J*u10K
                               A1111= (I[1]*P1J*(1-P1J)*P10J*P11K/nJ-I[3]*P1J*P1J*P11J*(1-P11J)/n1J)*u11J*u11K
                               A=sum(A0000,A0001,A0010,A0011,
                                     A0100,A0101,A0110,A0111,
                                     A1000,A1001,A1010,A1011,
                                     A1100,A1101,A1110,A1111)

                               if (FA=="binomial"){if (I[1]==1 & I[2]==1 & I[4]==1) {B0000=P0J*P0J*P00J*P00K*bv00J/n00J}  else {B0000=0}
                                                   if (I[1]==1 & I[2]==1 & I[5]==1) {B0101=P0J*P0J*P01J*P01K*bv01J/n01J}  else {B0101=0}
                                                   if (I[1]==1 & I[3]==1 & I[6]==1) {B1010=P1J*P1J*P10J*P10K*bv10J/n10J}  else {B1010=0}
                                                   if (I[1]==1 & I[3]==1 & I[7]==1) {B1111=P1J*P1J*P11J*P11K*bv11J/n11J}  else {B1111=0}
                                                   }
                                      else if (C) {if (I[1]==1 & I[2]==1 & I[4]==1) {B0000=P0J*P0J*P00J*P00K*comvar/n00J} else {B0000=0}
                                                   if (I[1]==1 & I[2]==1 & I[5]==1) {B0101=P0J*P0J*P01J*P01K*comvar/n01J} else {B0101=0}
                                                   if (I[1]==1 & I[3]==1 & I[6]==1) {B1010=P1J*P1J*P10J*P10K*comvar/n10J} else {B1010=0}
                                                   if (I[1]==1 & I[3]==1 & I[7]==1) {B1111=P1J*P1J*P11J*P11K*comvar/n11J} else {B1111=0}
                                                   }
                                             else {if (I[1]==1 & I[2]==1 & I[4]==1) {B0000=P0J*P0J*P00J*P00K*v00J/n00J}   else {B0000=0}
                                                   if (I[1]==1 & I[2]==1 & I[5]==1) {B0101=P0J*P0J*P01J*P01K*v01J/n01J}   else {B0101=0}
                                                   if (I[1]==1 & I[3]==1 & I[6]==1) {B1010=P1J*P1J*P10J*P10K*v10J/n10J}   else {B1010=0}
                                                   if (I[1]==1 & I[3]==1 & I[7]==1) {B1111=P1J*P1J*P11J*P11K*v11J/n11J}   else {B1111=0}
                                                   }
                               B=sum(B0000,B0101,B1010,B1111)
                               vc=A+B
                               } else
     if (Nstage==3 & Base==1) {DJ=D[which((D$O1==0 & D$A1==ATS1[1] & D$O2==0 & D$A2==ATS1[3] & D$O3==0 & D$A3==ATS1[7])|
                                          (D$O1==0 & D$A1==ATS1[1] & D$O2==0 & D$A2==ATS1[3] & D$O3==1 & D$A3==ATS1[8])|
                                          (D$O1==0 & D$A1==ATS1[1] & D$O2==1 & D$A2==ATS1[4] & D$O3==0 & D$A3==ATS1[9])|
                                          (D$O1==0 & D$A1==ATS1[1] & D$O2==1 & D$A2==ATS1[4] & D$O3==1 & D$A3==ATS1[10])|
                                          (D$O1==1 & D$A1==ATS1[2] & D$O2==0 & D$A2==ATS1[5] & D$O3==0 & D$A3==ATS1[11])|
                                          (D$O1==1 & D$A1==ATS1[2] & D$O2==0 & D$A2==ATS1[5] & D$O3==1 & D$A3==ATS1[12])|
                                          (D$O1==1 & D$A1==ATS1[2] & D$O2==1 & D$A2==ATS1[6] & D$O3==0 & D$A3==ATS1[13])|
                                          (D$O1==1 & D$A1==ATS1[2] & D$O2==1 & D$A2==ATS1[6] & D$O3==1 & D$A3==ATS1[14])),]
                               DK=D[which((D$O1==0 & D$A1==ATS2[1] & D$O2==0 & D$A2==ATS2[3] & D$O3==0 & D$A3==ATS2[7])|
                                          (D$O1==0 & D$A1==ATS2[1] & D$O2==0 & D$A2==ATS2[3] & D$O3==1 & D$A3==ATS2[8])|
                                          (D$O1==0 & D$A1==ATS2[1] & D$O2==1 & D$A2==ATS2[4] & D$O3==0 & D$A3==ATS2[9])|
                                          (D$O1==0 & D$A1==ATS2[1] & D$O2==1 & D$A2==ATS2[4] & D$O3==1 & D$A3==ATS2[10])|
                                          (D$O1==1 & D$A1==ATS2[2] & D$O2==0 & D$A2==ATS2[5] & D$O3==0 & D$A3==ATS2[11])|
                                          (D$O1==1 & D$A1==ATS2[2] & D$O2==0 & D$A2==ATS2[5] & D$O3==1 & D$A3==ATS2[12])|
                                          (D$O1==1 & D$A1==ATS2[2] & D$O2==1 & D$A2==ATS2[6] & D$O3==0 & D$A3==ATS2[13])|
                                          (D$O1==1 & D$A1==ATS2[2] & D$O2==1 & D$A2==ATS2[6] & D$O3==1 & D$A3==ATS2[14])),]

                               n0J=length(D$Y[which(D$O1==0 & D$A1==ATS1[1])]); if (n0J==0|is.na(n0J)) {n0J=1}
                               n1J=length(D$Y[which(D$O1==1 & D$A1==ATS1[2])]); if (n1J==0|is.na(n1J)) {n1J=1}

                               n00J=length(D$Y[which(D$O1==0 & D$A1==ATS1[1] & D$O2==0 & D$A2==ATS1[3])]); if (n00J==0|is.na(n00J)) {n00J=1}
                               n01J=length(D$Y[which(D$O1==0 & D$A1==ATS1[1] & D$O2==1 & D$A2==ATS1[4])]); if (n01J==0|is.na(n01J)) {n01J=1}
                               n10J=length(D$Y[which(D$O1==1 & D$A1==ATS1[2] & D$O2==0 & D$A2==ATS1[5])]); if (n10J==0|is.na(n10J)) {n10J=1}
                               n11J=length(D$Y[which(D$O1==1 & D$A1==ATS1[2] & D$O2==1 & D$A2==ATS1[6])]); if (n11J==0|is.na(n11J)) {n11J=1}

                               n000J=length(DJ$Y[which(D$O1==0 & D$O2==0 & D$O3==0)]); if (n000J==0|is.na(n000J)) {n000J=1}
                               n001J=length(DJ$Y[which(D$O1==0 & D$O2==0 & D$O3==1)]); if (n001J==0|is.na(n001J)) {n001J=1}
                               n010J=length(DJ$Y[which(D$O1==0 & D$O2==1 & D$O3==0)]); if (n010J==0|is.na(n010J)) {n010J=1}
                               n011J=length(DJ$Y[which(D$O1==0 & D$O2==1 & D$O3==1)]); if (n011J==0|is.na(n011J)) {n011J=1}
                               n100J=length(DJ$Y[which(D$O1==1 & D$O2==0 & D$O3==0)]); if (n100J==0|is.na(n100J)) {n100J=1}
                               n101J=length(DJ$Y[which(D$O1==1 & D$O2==0 & D$O3==1)]); if (n101J==0|is.na(n101J)) {n101J=1}
                               n110J=length(DJ$Y[which(D$O1==1 & D$O2==1 & D$O3==0)]); if (n110J==0|is.na(n110J)) {n110J=1}
                               n111J=length(DJ$Y[which(D$O1==1 & D$O2==1 & D$O3==1)]); if (n111J==0|is.na(n111J)) {n111J=1}

                               P0=length(D$Y[which(D$O1==0)])/length(D$Y)   #estimate respose rate P=P(R=1|A1=d1J)
                               P1=1-P0

                               P00J=length(D$Y[which(D$O1==0 & D$A1==ATS1[1] & D$O2==0)])/
                                    length(D$Y[which(D$O1==0 & D$A1==ATS1[1])])
                               P01J=1-P00J
                               P10J=length(D$Y[which(D$O1==1 & D$A1==ATS1[2] & D$O2==0)])/
                                    length(D$Y[which(D$O1==1 & D$A1==ATS1[2])])
                               P11J=1-P10J

                               P00K=length(D$Y[which(D$O1==0 & D$A1==ATS2[1] & D$O2==0)])/
                                    length(D$Y[which(D$O1==0 & D$A1==ATS2[1])])
                               P01K=1-P00K
                               P10K=length(D$Y[which(D$O1==1 & D$A1==ATS2[2] & D$O2==0)])/
                                    length(D$Y[which(D$O1==1 & D$A1==ATS2[2])])
                               P11K=1-P10K

                               P000J=length(D$Y[which(D$O1==0 & D$A1==ATS1[1] & D$O2==0 & D$A2==ATS1[3] & D$O3==0)])/
                                     length(D$Y[which(D$O1==0 & D$A1==ATS1[1] & D$O2==0 & D$A2==ATS1[3])])
                               P001J=1-P000J
                               P010J=length(D$Y[which(D$O1==0 & D$A1==ATS1[1] & D$O2==1 & D$A2==ATS1[4] & D$O3==0)])/
                                     length(D$Y[which(D$O1==0 & D$A1==ATS1[1] & D$O2==1 & D$A2==ATS1[4])])
                               P011J=1-P010J
                               P100J=length(D$Y[which(D$O1==1 & D$A1==ATS1[2] & D$O2==0 & D$A2==ATS1[5] & D$O3==0)])/
                                     length(D$Y[which(D$O1==1 & D$A1==ATS1[2] & D$O2==0 & D$A2==ATS1[5])])
                               P101J=1-P100J
                               P110J=length(D$Y[which(D$O1==1 & D$A1==ATS1[2] & D$O2==1 & D$A2==ATS1[6] & D$O3==0)])/
                                     length(D$Y[which(D$O1==1 & D$A1==ATS1[2] & D$O2==1 & D$A2==ATS1[6])])
                               P111J=1-P110J

                               P000K=length(D$Y[which(D$O1==0 & D$A1==ATS2[1] & D$O2==0 & D$A2==ATS2[3] & D$O3==0)])/
                                     length(D$Y[which(D$O1==0 & D$A1==ATS2[1] & D$O2==0 & D$A2==ATS2[3])])
                               P001K=1-P000K
                               P010K=length(D$Y[which(D$O1==0 & D$A1==ATS2[1] & D$O2==1 & D$A2==ATS2[4] & D$O3==0)])/
                                     length(D$Y[which(D$O1==0 & D$A1==ATS2[1] & D$O2==1 & D$A2==ATS2[4])])
                               P011K=1-P010K
                               P100K=length(D$Y[which(D$O1==1 & D$A1==ATS2[2] & D$O2==0 & D$A2==ATS2[5] & D$O3==0)])/
                                     length(D$Y[which(D$O1==1 & D$A1==ATS2[2] & D$O2==0 & D$A2==ATS2[5])])
                               P101K=1-P100K
                               P110K=length(D$Y[which(D$O1==1 & D$A1==ATS2[2] & D$O2==1 & D$A2==ATS2[6] & D$O3==0)])/
                                     length(D$Y[which(D$O1==1 & D$A1==ATS2[2] & D$O2==1 & D$A2==ATS2[6])])
                               P111K=1-P110K

                               u000J=mean(DJ$Y[which(DJ$O1==0 & DJ$O2==0 & DJ$O3==0)]); if (is.na(u000J)) {u000J=0}        # within-sequence true mean (EY)
                               u001J=mean(DJ$Y[which(DJ$O1==0 & DJ$O2==0 & DJ$O3==1)]); if (is.na(u001J)) {u001J=0}
                               u010J=mean(DJ$Y[which(DJ$O1==0 & DJ$O2==1 & DJ$O3==0)]); if (is.na(u010J)) {u010J=0}
                               u011J=mean(DJ$Y[which(DJ$O1==0 & DJ$O2==1 & DJ$O3==1)]); if (is.na(u011J)) {u011J=0}
                               u100J=mean(DJ$Y[which(DJ$O1==1 & DJ$O2==0 & DJ$O3==0)]); if (is.na(u100J)) {u100J=0}
                               u101J=mean(DJ$Y[which(DJ$O1==1 & DJ$O2==0 & DJ$O3==1)]); if (is.na(u101J)) {u101J=0}
                               u110J=mean(DJ$Y[which(DJ$O1==1 & DJ$O2==1 & DJ$O3==0)]); if (is.na(u110J)) {u110J=0}
                               u111J=mean(DJ$Y[which(DJ$O1==1 & DJ$O2==1 & DJ$O3==1)]); if (is.na(u111J)) {u111J=0}

                               u000K=mean(DK$Y[which(DK$O1==0 & DK$O2==0 & DK$O3==0)]); if (is.na(u000K)) {u000K=0}        # within-sequence true mean (EY)
                               u001K=mean(DK$Y[which(DK$O1==0 & DK$O2==0 & DK$O3==1)]); if (is.na(u001K)) {u001K=0}
                               u010K=mean(DK$Y[which(DK$O1==0 & DK$O2==1 & DK$O3==0)]); if (is.na(u010K)) {u010K=0}
                               u011K=mean(DK$Y[which(DK$O1==0 & DK$O2==1 & DK$O3==1)]); if (is.na(u011K)) {u011K=0}
                               u100K=mean(DK$Y[which(DK$O1==1 & DK$O2==0 & DK$O3==0)]); if (is.na(u100K)) {u100K=0}
                               u101K=mean(DK$Y[which(DK$O1==1 & DK$O2==0 & DK$O3==1)]); if (is.na(u101K)) {u101K=0}
                               u110K=mean(DK$Y[which(DK$O1==1 & DK$O2==1 & DK$O3==0)]); if (is.na(u110K)) {u110K=0}
                               u111K=mean(DK$Y[which(DK$O1==1 & DK$O2==1 & DK$O3==1)]); if (is.na(u111K)) {u111K=0}

                               v000J=var(DJ$Y[which(DJ$O1==0 & DJ$O2==0 & DJ$O3==0)]);  if (is.na(v000J)) {v000J=0}         #within-sequence variance for cont outcome
                               v001J=var(DJ$Y[which(DJ$O1==0 & DJ$O2==0 & DJ$O3==1)]);  if (is.na(v001J)) {v001J=0}
                               v010J=var(DJ$Y[which(DJ$O1==0 & DJ$O2==1 & DJ$O3==0)]);  if (is.na(v010J)) {v010J=0}
                               v011J=var(DJ$Y[which(DJ$O1==0 & DJ$O2==1 & DJ$O3==1)]);  if (is.na(v011J)) {v011J=0}
                               v100J=var(DJ$Y[which(DJ$O1==1 & DJ$O2==0 & DJ$O3==0)]);  if (is.na(v100J)) {v100J=0}
                               v101J=var(DJ$Y[which(DJ$O1==1 & DJ$O2==0 & DJ$O3==1)]);  if (is.na(v101J)) {v101J=0}
                               v110J=var(DJ$Y[which(DJ$O1==1 & DJ$O2==1 & DJ$O3==0)]);  if (is.na(v110J)) {v110J=0}
                               v111J=var(DJ$Y[which(DJ$O1==1 & DJ$O2==1 & DJ$O3==1)]);  if (is.na(v111J)) {v111J=0}

                               bv000J=u000J*(1-u000J);  if (is.na(bv000J)) {bv000J=0}                                       #within-sequence var est for bin outcome
                               bv001J=u001J*(1-u001J);  if (is.na(bv001J)) {bv001J=0}
                               bv010J=u010J*(1-u010J);  if (is.na(bv010J)) {bv010J=0}
                               bv011J=u011J*(1-u011J);  if (is.na(bv011J)) {bv011J=0}
                               bv100J=u100J*(1-u100J);  if (is.na(bv100J)) {bv100J=0}
                               bv101J=u101J*(1-u101J);  if (is.na(bv101J)) {bv101J=0}
                               bv110J=u110J*(1-u110J);  if (is.na(bv110J)) {bv110J=0}
                               bv111J=u111J*(1-u111J);  if (is.na(bv111J)) {bv111J=0}
                               comvar=cv(data=D);  if (is.na(comvar)) {comvar=0}

                               A000.000=( P0*(1-P0)*P00J*P00K*P000J*P000K/N+I[1]*P0*P0*P00J*(1-P00J)*P000J*P000K/n0J+I[1]*I[3]*P0*P0*P00J*P00K*P000J*(1-P000J)/n00J)*u000J*u000K
                               A000.001=( P0*(1-P0)*P00J*P00K*P000J*P001K/N+I[1]*P0*P0*P00J*(1-P00J)*P000J*P001K/n0J-I[1]*I[3]*P0*P0*P00J*P00K*P000J*(1-P000J)/n00J)*u000J*u001K
                               A000.010=( P0*(1-P0)*P00J*P10K*P000J*P010K/N-I[1]*P0*P0*P00J*(1-P00J)*P000J*P010K/n0J)*u000J*u010K
                               A000.011=( P0*(1-P0)*P00J*P01K*P000J*P011K/N-I[1]*P0*P0*P00J*(1-P00J)*P000J*P011K/n0J)*u000J*u011K
                               A000.100=(-P0*(1-P0)*P00J*P10K*P000J*P100K/N)*u000J*u100K
                               A000.101=(-P0*(1-P0)*P00J*P10K*P000J*P101K/N)*u000J*u101K
                               A000.110=(-P0*(1-P0)*P00J*P11K*P000J*P110K/N)*u000J*u110K
                               A000.111=(-P0*(1-P0)*P00J*P11K*P000J*P111K/N)*u000J*u111K

                               A001.000=( P0*(1-P0)*P00J*P00K*P001J*P000K/N+I[1]*P0*P0*P00J*(1-P00J)*P001J*P000K/n0J-I[1]*I[3]*P0*P0*P00J*P00K*P001J*(1-P001J)/n00J)*u001J*u000K
                               A001.001=( P0*(1-P0)*P00J*P00K*P001J*P001K/N+I[1]*P0*P0*P00J*(1-P00J)*P001J*P001K/n0J+I[1]*I[3]*P0*P0*P00J*P00K*P001J*(1-P001J)/n00J)*u001J*u001K
                               A001.010=( P0*(1-P0)*P00J*P10K*P001J*P010K/N-I[1]*P0*P0*P00J*(1-P00J)*P001J*P010K/n0J)*u001J*u010K
                               A001.011=( P0*(1-P0)*P00J*P10K*P001J*P011K/N-I[1]*P0*P0*P00J*(1-P00J)*P001J*P011K/n0J)*u001J*u011K
                               A001.100=(-P0*(1-P0)*P00J*P10K*P001J*P100K/N)*u001J*u100K
                               A001.101=(-P0*(1-P0)*P00J*P10K*P001J*P101K/N)*u001J*u101K
                               A001.110=(-P0*(1-P0)*P00J*P11K*P001J*P110K/N)*u001J*u110K
                               A001.111=(-P0*(1-P0)*P00J*P11K*P001J*P111K/N)*u001J*u111K

                               A010.000=( P0*(1-P0)*P01J*P00K*P010J*P000K/N-I[1]*P0*P0*P01J*(1-P01J)*P010J*P000K/n0J)*u010J*u000K
                               A010.001=( P0*(1-P0)*P01J*P00K*P010J*P001K/N-I[1]*P0*P0*P01J*(1-P01J)*P010J*P001K/n0J)*u010J*u001K
                               A010.010=( P0*(1-P0)*P01J*P01K*P010J*P010K/N+I[1]*P0*P0*P01J*(1-P01J)*P010J*P010K/n0J+I[1]*I[3]*P0*P0*P01J*P01K*P010J*(1-P010J)/n01J)*u010J*u010K
                               A010.011=( P0*(1-P0)*P01J*P01K*P010J*P011K/N+I[1]*P0*P0*P01J*(1-P01J)*P010J*P011K/n0J-I[1]*I[3]*P0*P0*P01J*P01K*P010J*(1-P010J)/n01J)*u010J*u011K
                               A010.100=(-P0*(1-P0)*P01J*P10K*P010J*P100K/N)*u010J*u100K
                               A010.101=(-P0*(1-P0)*P01J*P10K*P010J*P101K/N)*u010J*u101K
                               A010.110=(-P0*(1-P0)*P01J*P11K*P010J*P110K/N)*u010J*u110K
                               A010.111=(-P0*(1-P0)*P01J*P11K*P010J*P111K/N)*u010J*u111K

                               A011.000=( P0*(1-P0)*P01J*P00K*P011J*P000K/N-I[1]*P0*P0*P01J*(1-P01J)*P011J*P000K/n0J)*u011J*u000K
                               A011.001=( P0*(1-P0)*P01J*P00K*P011J*P001K/N-I[1]*P0*P0*P01J*(1-P01J)*P011J*P001K/n0J)*u011J*u001K
                               A011.010=( P0*(1-P0)*P01J*P01K*P011J*P010K/N+I[1]*P0*P0*P01J*(1-P01J)*P011J*P010K/n0J-I[1]*I[3]*P0*P0*P01J*P01K*P011J*(1-P011J)/n01J)*u011J*u010K
                               A011.011=( P0*(1-P0)*P01J*P01K*P011J*P011K/N+I[1]*P0*P0*P01J*(1-P01J)*P011J*P011K/n0J+I[1]*I[3]*P0*P0*P01J*P01K*P011J*(1-P011J)/n01J)*u011J*u011K
                               A011.100=(-P0*(1-P0)*P01J*P10K*P011J*P100K/N)*u011J*u100K
                               A011.101=(-P0*(1-P0)*P01J*P10K*P011J*P101K/N)*u011J*u101K
                               A011.110=(-P0*(1-P0)*P01J*P11K*P011J*P110K/N)*u011J*u110K
                               A011.111=(-P0*(1-P0)*P01J*P11K*P011J*P111K/N)*u011J*u111K

                               A100.000=(-P1*(1-P1)*P10J*P00K*P100J*P000K/N)*u100J*u000K
                               A100.001=(-P1*(1-P1)*P10J*P00K*P100J*P001K/N)*u100J*u001K
                               A100.010=(-P1*(1-P1)*P10J*P01K*P100J*P010K/N)*u100J*u010K
                               A100.011=(-P1*(1-P1)*P10J*P01K*P100J*P011K/N)*u100J*u011K
                               A100.100=( P1*(1-P1)*P10J*P10K*P100J*P100K/N+I[1]*P1*P1*P10J*(1-P10J)*P100J*P100K/n1J+I[1]*I[3]*P1*P1*P10J*P10K*P100J*(1-P100J)/n10J)*u100J*u100K
                               A100.101=( P1*(1-P1)*P10J*P10K*P100J*P101K/N+I[1]*P1*P1*P10J*(1-P10J)*P100J*P101K/n1J-I[1]*I[3]*P1*P1*P10J*P10K*P100J*(1-P100J)/n10J)*u100J*u101K
                               A100.110=( P1*(1-P1)*P10J*P11K*P100J*P110K/N-I[1]*P1*P1*P10J*(1-P10J)*P100J*P110K/n1J)*u100J*u110K
                               A100.111=( P1*(1-P1)*P10J*P11K*P100J*P111K/N-I[1]*P1*P1*P10J*(1-P10J)*P100J*P111K/n1J)*u100J*u111K

                               A101.000=(-P1*(1-P1)*P10J*P00K*P101J*P000K/N)*u101J*u000K
                               A101.001=(-P1*(1-P1)*P10J*P00K*P101J*P001K/N)*u101J*u001K
                               A101.010=(-P1*(1-P1)*P10J*P01K*P101J*P010K/N)*u101J*u010K
                               A101.011=(-P1*(1-P1)*P10J*P01K*P101J*P011K/N)*u101J*u011K
                               A101.100=( P1*(1-P1)*P10J*P10K*P101J*P100K/N+I[1]*P1*P1*P10J*(1-P10J)*P101J*P100K/n1J-I[1]*I[3]*P1*P1*P10J*P10K*P101J*(1-P101J)/n10J)*u101J*u100K
                               A101.101=( P1*(1-P1)*P10J*P10K*P101J*P101K/N+I[1]*P1*P1*P10J*(1-P10J)*P101J*P101K/n1J+I[1]*I[3]*P1*P1*P10J*P10K*P101J*(1-P101J)/n10J)*u101J*u101K
                               A101.110=( P1*(1-P1)*P10J*P11K*P101J*P110K/N+I[1]*P1*P1*P10J*(1-P10J)*P101J*P110K/n1J)*u101J*u110K
                               A101.111=( P1*(1-P1)*P10J*P11K*P101J*P111K/N+I[1]*P1*P1*P10J*(1-P10J)*P101J*P111K/n1J)*u101J*u111K

                               A110.000=(-P1*(1-P1)*P11J*P00K*P110J*P000K/N)*u110J*u000K
                               A110.001=(-P1*(1-P1)*P11J*P00K*P110J*P001K/N)*u110J*u001K
                               A110.010=(-P1*(1-P1)*P11J*P01K*P110J*P010K/N)*u110J*u010K
                               A110.011=(-P1*(1-P1)*P11J*P01K*P110J*P011K/N)*u110J*u011K
                               A110.100=( P1*(1-P1)*P11J*P10K*P110J*P100K/N-I[1]*P1*P1*P11J*(1-P11J)*P110J*P100K/n1J)*u110J*u100K
                               A110.101=( P1*(1-P1)*P11J*P10K*P110J*P101K/N-I[1]*P1*P1*P11J*(1-P11J)*P110J*P101K/n1J)*u110J*u101K
                               A110.110=( P1*(1-P1)*P11J*P11K*P110J*P110K/N+I[1]*P1*P1*P11J*(1-P11J)*P110J*P110K/n1J-I[1]*I[3]*P1*P1*P11J*P11K*P110J*(1-P110J)/n10J)*u110J*u110K
                               A110.111=( P1*(1-P1)*P11J*P11K*P110J*P111K/N+I[1]*P1*P1*P11J*(1-P11J)*P111J*P111K/n1J-I[1]*I[3]*P1*P1*P11J*P11K*P110J*(1-P110J)/n10J)*u110J*u111K

                               A111.000=(-P1*(1-P1)*P11J*P00K*P111J*P000K/N)*u111J*u000K
                               A111.001=(-P1*(1-P1)*P11J*P00K*P111J*P000K/N)*u111J*u001K
                               A111.010=(-P1*(1-P1)*P11J*P00K*P111J*P000K/N)*u111J*u010K
                               A111.011=(-P1*(1-P1)*P11J*P00K*P111J*P000K/N)*u111J*u011K
                               A111.100=( P1*(1-P1)*P11J*P10K*P111J*P100K/N-I[1]*P1*P1*P11J*(1-P11J)*P111J*P100K/n1J)*u111J*u100K
                               A111.101=( P1*(1-P1)*P11J*P10K*P111J*P101K/N-I[1]*P1*P1*P11J*(1-P11J)*P111J*P101K/n1J)*u111J*u101K
                               A111.110=( P1*(1-P1)*P11J*P11K*P111J*P110K/N+I[1]*P1*P1*P11J*(1-P11J)*P111J*P110K/n1J-I[1]*I[3]*P1*P1*P11J*P11K*P111J*(1-P111J)/n10J)*u111J*u110K
                               A111.111=( P1*(1-P1)*P11J*P11K*P111J*P111K/N+I[1]*P1*P1*P11J*(1-P11J)*P111J*P111K/n1J+I[1]*I[3]*P1*P1*P11J*P11K*P111J*(1-P111J)/n10J)*u111J*u111K

                               A=sum(A000.000,A000.001,A000.010,A000.011,A000.100,A000.101,A000.110,A000.111,
                                     A001.000,A001.001,A001.010,A001.011,A001.100,A001.101,A001.110,A001.111,
                                     A010.000,A010.001,A010.010,A010.011,A010.100,A010.101,A010.110,A010.111,
                                     A011.000,A011.001,A011.010,A011.011,A011.100,A011.101,A011.110,A011.111,
                                     A100.000,A100.001,A100.010,A100.011,A100.100,A100.101,A100.110,A100.111,
                                     A101.000,A101.001,A101.010,A101.011,A101.100,A101.101,A101.110,A101.111,
                                     A110.000,A110.001,A110.010,A110.011,A110.100,A110.101,A110.110,A110.111,
                                     A111.000,A111.001,A111.010,A111.011,A111.100,A111.101,A111.110,A111.111)

                               if (FA=="binomial"){if (I[1]==1 & I[3]==1 & I[ 7]==1) {B000.000=P0*P0*P00J*P00K*P000J*P000K*bv000J/n000J}  else {B000.000=0}
                                                   if (I[1]==1 & I[3]==1 & I[ 8]==1) {B001.001=P0*P0*P00J*P00K*P001J*P001K*bv001J/n001J}  else {B001.001=0}
                                                   if (I[1]==1 & I[4]==1 & I[ 9]==1) {B010.010=P0*P0*P01J*P01K*P010J*P010K*bv010J/n010J}  else {B010.010=0}
                                                   if (I[1]==1 & I[4]==1 & I[10]==1) {B011.011=P0*P0*P01J*P01K*P011J*P011K*bv011J/n011J}  else {B011.011=0}
                                                   if (I[1]==1 & I[5]==1 & I[11]==1) {B100.100=P1*P1*P10J*P10K*P100J*P100K*bv100J/n100J}  else {B100.100=0}
                                                   if (I[1]==1 & I[5]==1 & I[12]==1) {B101.101=P1*P1*P10J*P10K*P101J*P101K*bv101J/n101J}  else {B101.101=0}
                                                   if (I[1]==1 & I[6]==1 & I[13]==1) {B110.110=P1*P1*P11J*P11K*P110J*P110K*bv110J/n110J}  else {B110.110=0}
                                                   if (I[1]==1 & I[6]==1 & I[14]==1) {B111.111=P1*P1*P11J*P11K*P111J*P111K*bv111J/n111J}  else {B111.111=0}
                                                   }
                                      else if (C) {if (I[1]==1 & I[3]==1 & I[ 7]==1) {B000.000=P0*P0*P00J*P00K*P000J*P000K*comvar/n000J} else {B000.000=0}
                                                   if (I[1]==1 & I[3]==1 & I[ 8]==1) {B001.001=P0*P0*P01J*P01K*P001J*P001K*comvar/n001J} else {B001.001=0}
                                                   if (I[1]==1 & I[4]==1 & I[ 9]==1) {B010.010=P1*P1*P10J*P10K*P010J*P010K*comvar/n010J} else {B010.010=0}
                                                   if (I[1]==1 & I[4]==1 & I[10]==1) {B011.011=P1*P1*P11J*P11K*P011J*P011K*comvar/n011J} else {B011.011=0}
                                                   if (I[1]==1 & I[5]==1 & I[11]==1) {B100.100=P0*P0*P00J*P00K*P100J*P100K*comvar/n100J} else {B100.100=0}
                                                   if (I[1]==1 & I[5]==1 & I[12]==1) {B101.101=P0*P0*P01J*P01K*P101J*P101K*comvar/n101J} else {B101.101=0}
                                                   if (I[1]==1 & I[6]==1 & I[13]==1) {B110.110=P1*P1*P10J*P10K*P110J*P110K*comvar/n110J} else {B110.110=0}
                                                   if (I[1]==1 & I[6]==1 & I[14]==1) {B111.111=P1*P1*P11J*P11K*P111J*P111K*comvar/n111J} else {B111.111=0}
                                                   }
                                             else {if (I[1]==1 & I[3]==1 & I[ 7]==1) {B000.000=P0*P0*P00J*P00K*P000J*P000K*v000J/n000J}  else {B000.000=0}
                                                   if (I[1]==1 & I[3]==1 & I[ 8]==1) {B001.001=P0*P0*P01J*P01K*P001J*P001K*v001J/n001J}  else {B001.001=0}
                                                   if (I[1]==1 & I[4]==1 & I[ 9]==1) {B010.010=P1*P1*P10J*P10K*P010J*P010K*v010J/n010J}  else {B010.010=0}
                                                   if (I[1]==1 & I[4]==1 & I[10]==1) {B011.011=P1*P1*P11J*P11K*P011J*P011K*v011J/n011J}  else {B011.011=0}
                                                   if (I[1]==1 & I[5]==1 & I[11]==1) {B100.100=P0*P0*P00J*P00K*P100J*P100K*v100J/n100J}  else {B100.100=0}
                                                   if (I[1]==1 & I[5]==1 & I[12]==1) {B101.101=P0*P0*P01J*P01K*P101J*P101K*v101J/n101J}  else {B101.101=0}
                                                   if (I[1]==1 & I[6]==1 & I[13]==1) {B110.110=P1*P1*P10J*P10K*P110J*P110K*v110J/n110J}  else {B110.110=0}
                                                   if (I[1]==1 & I[6]==1 & I[14]==1) {B111.111=P1*P1*P11J*P11K*P111J*P111K*v111J/n111J}  else {B111.111=0}
                                                   }
                               B=sum(B000.000,B001.001,B010.010,B011.011,
                                     B100.100,B101.101,B110.110,B111.111)
                               vc=A+B
                               }
     return(vc)
}



