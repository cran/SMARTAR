#tcov() to get var-cov component for effect size
#input two strategies and sequence information matrix, also the method to estimate

tcov=function(ats1,ats2,sim,family="normal",method="Gest"){
     ATS1=ats1; ATS2=ats2; Smat=sim; FA=family; MA=method
     if (is.null(Smat$O1)) {Base=0} else {Base=1}
     Nstage=nstage(data=Smat)

     I=overlap(ats1=ATS1,ats2=ATS2,nstage=Nstage,baseline=Base)
     if (Nstage==1 & Base==0) {pi1=Smat$PI1[which(Smat$A1==ATS1[1])]   #P(A1==1)
                               uJ=Smat$MEAN[which(Smat$A1==ATS1[1])]
                               if (FA=="normal") {vJ=Smat$SD[which(Smat$A1==ATS1[1])]^2} else
                               if (FA=="binomial") {vJ=uJ*(1-uJ)}
                               tvc=I[1]*vJ/pi1
                               } else
     if (Nstage==1 & Base==1) {pr=mean(Smat$P1[which(Smat$O1==1)])                   #P(O1==1)
                               pi0=Smat$PI1[which(Smat$O1==1 & Smat$A1==ATS1[1])]    #P(A1==1|O1==0)
                               pi1=Smat$PI1[which(Smat$O1==1 & Smat$A1==ATS1[1])]    #P(A1==1|O1==1)

                               u0J=Smat$MEAN[which(Smat$O1==0 & Smat$A1==ATS1[1])]; if (is.na(u0J)) {u0J=0}                        # within-sequence true mean (EY)
                               u1J=Smat$MEAN[which(Smat$O1==1 & Smat$A1==ATS1[2])]; if (is.na(u1J)) {u1J=0}
                               u0K=Smat$MEAN[which(Smat$O1==0 & Smat$A1==ATS2[1])]; if (is.na(u0K)) {u0K=0}
                               u1K=Smat$MEAN[which(Smat$O1==1 & Smat$A1==ATS2[2])]; if (is.na(u1K)) {u1K=0}

                               if (FA=="normal"){v0J=Smat$SD[which(Smat$O1==0 & Smat$A1==ATS1[1])]^2                         #within-sequence variance for cont outcome
                                                 v1J=Smat$SD[which(Smat$O1==1 & Smat$A1==ATS1[2])]^2
                                                 v0K=Smat$SD[which(Smat$O1==0 & Smat$A1==ATS2[1])]^2
                                                 v1K=Smat$SD[which(Smat$O1==1 & Smat$A1==ATS2[2])]^2} else
                               if (FA=="binomial"){v0J=u0J*(1-u0J)
                                                   v1J=u1J*(1-u1J)
                                                   v0K=u0K*(1-u0K)
                                                   v1K=u1K*(1-u1K)}
                               if (is.na(v0J)) {v0J=0}
                               if (is.na(v1J)) {v1J=0}
                               if (is.na(v0K)) {v0K=0}
                               if (is.na(v1K)) {v1K=0}

                               if (MA=="Gest") {B1=1*pr*(1-pr)*(u1J-u0J)*(u1K-u0K)          #MLE for cont outcome without com assumption
                                                    B2=I[1]*(1-pr)*v0J/pi0
                                                    B3=I[2]*pr*v1J/pi1
                                                    tvc=B1+B2+B3}
                               } else
     if (Nstage==2 & Base==0) {pi1=mean(Smat$PI1[which(Smat$A1==ATS1[1])])
                               pr =mean(Smat$P2[which(Smat$A1==ATS1[1] & Smat$O2==1)])
                               pi20J=mean(Smat$PI2[which(Smat$A1==ATS1[1] & Smat$O2==0 & Smat$A2==ATS1[2])])
                               pi21J=mean(Smat$PI2[which(Smat$A1==ATS1[1] & Smat$O2==1 & Smat$A2==ATS1[3])])

                               u0J=Smat$MEAN[which(Smat$A1==ATS1[1] & Smat$O2==0 & Smat$A2==ATS1[2])]
                               u1J=Smat$MEAN[which(Smat$A1==ATS1[1] & Smat$O2==1 & Smat$A2==ATS1[3])]
                               u0K=Smat$MEAN[which(Smat$A1==ATS2[1] & Smat$O2==0 & Smat$A2==ATS2[2])]
                               u1K=Smat$MEAN[which(Smat$A1==ATS2[1] & Smat$O2==1 & Smat$A2==ATS2[3])]

                               if (FA=="normal"){v0J=Smat$SD[which(Smat$A1==ATS1[1] & Smat$O2==0 & Smat$A2==ATS1[2])]^2
                                                 v1J=Smat$SD[which(Smat$A1==ATS1[1] & Smat$O2==1 & Smat$A2==ATS1[3])]^2
                                                 v0K=Smat$SD[which(Smat$A1==ATS2[1] & Smat$O2==0 & Smat$A2==ATS2[2])]^2
                                                 v1K=Smat$SD[which(Smat$A1==ATS2[1] & Smat$O2==1 & Smat$A2==ATS2[3])]^2} else
                               if (FA=="binomial"){v0J=u0J*(1-u0J)
                                                   v1J=u1J*(1-u1J)
                                                   v0K=u0K*(1-u0K)
                                                   v1K=u1K*(1-u1K)}

                               if (is.na(u0J)) {u0J=0}
                               if (is.na(u1J)) {u1J=0}
                               if (is.na(u0K)) {u0K=0}
                               if (is.na(u1K)) {u1K=0}

                               if (length(v0J)==0L) {v0J=0}
                               if (length(v1J)==0L) {v1J=0}
                               if (length(v0K)==0L) {v0K=0}
                               if (length(v1K)==0L) {v1K=0}

                               if (is.nan(pr)) {pr=0}
                               if (is.nan(pi1)) {pi1=1}
                               if (is.nan(pi20J)) {pi20J=1}
                               if (is.nan(pi21J)) {pi21J=1}

                               if (MA=="Gest"){B1=I[1]*(u1J-u0J)*(u1K-u0K)*pr*(1-pr)/pi1
                                               B2=I[2]*v0J*(1-pr)/(pi1*pi20J)
                                               B3=I[3]*v1J*pr    /(pi1*pi21J)
                                               tvc=B1+B2+B3}
                               } else
     if (Nstage==2 & Base==1) {P0=mean(Smat$P1[which(Smat$O1==0)]); if (is.nan(P0)) {P0=0}
                               P1=mean(Smat$P1[which(Smat$O1==1)]); if (is.nan(P1)) {P1=0}

                               pi0J=mean(Smat$PI1[which(Smat$O1==0 & Smat$A1==ATS1[1])]); if (is.nan(pi0J)) {pi0J=1}
                               pi1J=mean(Smat$PI1[which(Smat$O1==1 & Smat$A1==ATS1[2])]); if (is.nan(pi1J)) {pi1J=1}

                               P00J=mean(Smat$P2[which(Smat$O1==0 & Smat$A1==ATS1[1] & Smat$O2==0)]); if (is.nan(P00J)) {P00J=0}
                               P01J=mean(Smat$P2[which(Smat$O1==0 & Smat$A1==ATS1[1] & Smat$O2==1)]); if (is.nan(P01J)) {P01J=0}
                               P10J=mean(Smat$P2[which(Smat$O1==1 & Smat$A1==ATS1[2] & Smat$O2==0)]); if (is.nan(P10J)) {P10J=0}
                               P11J=mean(Smat$P2[which(Smat$O1==1 & Smat$A1==ATS1[2] & Smat$O2==1)]); if (is.nan(P11J)) {P11J=0}

                               P00K=mean(Smat$P2[which(Smat$O1==0 & Smat$A1==ATS2[1] & Smat$O2==0)]); if (is.nan(P00K)) {P00K=0}
                               P01K=mean(Smat$P2[which(Smat$O1==0 & Smat$A1==ATS2[1] & Smat$O2==1)]); if (is.nan(P01K)) {P01K=0}
                               P10K=mean(Smat$P2[which(Smat$O1==1 & Smat$A1==ATS2[2] & Smat$O2==0)]); if (is.nan(P10K)) {P10K=0}
                               P11K=mean(Smat$P2[which(Smat$O1==1 & Smat$A1==ATS2[2] & Smat$O2==1)]); if (is.nan(P11K)) {P11K=0}

                               pi00J=mean(Smat$PI2[which(Smat$O1==0 & Smat$A1==ATS1[1] & Smat$O2==0 & Smat$A2==ATS1[3])]); if (is.nan(pi00J)) {pi00J=1}
                               pi01J=mean(Smat$PI2[which(Smat$O1==0 & Smat$A1==ATS1[1] & Smat$O2==1 & Smat$A2==ATS1[4])]); if (is.nan(pi01J)) {pi01J=1}
                               pi10J=mean(Smat$PI2[which(Smat$O1==1 & Smat$A1==ATS1[2] & Smat$O2==0 & Smat$A2==ATS1[5])]); if (is.nan(pi10J)) {pi10J=1}
                               pi11J=mean(Smat$PI2[which(Smat$O1==1 & Smat$A1==ATS1[2] & Smat$O2==1 & Smat$A2==ATS1[6])]); if (is.nan(pi11J)) {pi11J=1}

                               u00J=Smat$MEAN[which(Smat$O1==0 & Smat$A1==ATS1[1] & Smat$O2==0 & Smat$A2==ATS1[3])]; if (is.na(u00J)) {u00J=0}
                               u01J=Smat$MEAN[which(Smat$O1==0 & Smat$A1==ATS1[1] & Smat$O2==1 & Smat$A2==ATS1[4])]; if (is.na(u01J)) {u01J=0}
                               u10J=Smat$MEAN[which(Smat$O1==1 & Smat$A1==ATS1[2] & Smat$O2==0 & Smat$A2==ATS1[5])]; if (is.na(u10J)) {u10J=0}
                               u11J=Smat$MEAN[which(Smat$O1==1 & Smat$A1==ATS1[2] & Smat$O2==1 & Smat$A2==ATS1[6])]; if (is.na(u11J)) {u11J=0}

                               u00K=Smat$MEAN[which(Smat$O1==0 & Smat$A1==ATS2[1] & Smat$O2==0 & Smat$A2==ATS2[3])]; if (is.na(u00K)) {u00K=0}
                               u01K=Smat$MEAN[which(Smat$O1==0 & Smat$A1==ATS2[1] & Smat$O2==1 & Smat$A2==ATS2[4])]; if (is.na(u01K)) {u01K=0}
                               u10K=Smat$MEAN[which(Smat$O1==1 & Smat$A1==ATS2[2] & Smat$O2==0 & Smat$A2==ATS2[5])]; if (is.na(u10K)) {u10K=0}
                               u11K=Smat$MEAN[which(Smat$O1==1 & Smat$A1==ATS2[2] & Smat$O2==1 & Smat$A2==ATS2[6])]; if (is.na(u11K)) {u11K=0}

                               if (FA=="normal"){v00J=Smat$SD[which(Smat$O1==0 & Smat$A1==ATS1[1] & Smat$O2==0 & Smat$A2==ATS1[3])]^2
                                                 v01J=Smat$SD[which(Smat$O1==0 & Smat$A1==ATS1[1] & Smat$O2==1 & Smat$A2==ATS1[4])]^2
                                                 v10J=Smat$SD[which(Smat$O1==1 & Smat$A1==ATS1[2] & Smat$O2==0 & Smat$A2==ATS1[5])]^2
                                                 v11J=Smat$SD[which(Smat$O1==1 & Smat$A1==ATS1[2] & Smat$O2==1 & Smat$A2==ATS1[6])]^2} else
                               if (FA=="binomial"){v00J=u00J*(1-u00J)
                                                   v01J=u01J*(1-u01J)
                                                   v10J=u10J*(1-u10J)
                                                   v11J=u11J*(1-u11J)}
                               if (length(v00J)==0L) {v00J=0}
                               if (length(v01J)==0L) {v01J=0}
                               if (length(v10J)==0L) {v10J=0}
                               if (length(v11J)==0L) {v11J=0}


                               A0000= (P0*(1-P0)*P00J*P00K+I[1]*P0*P0*P00J*(1-P00J)/(P0*pi0J))*u00J*u00K
                               A0001= (P0*(1-P0)*P00J*P01K-I[1]*P0*P0*P00J*(1-P00J)/(P0*pi0J))*u00J*u01K
                               A0010=-(P0*(1-P0)*P00J*P10K)*u00J*u10K
                               A0011=-(P0*(1-P0)*P00J*P11K)*u00J*u11K

                               A0100= (P0*(1-P0)*P01J*P00K-I[1]*P0*P0*P01J*(1-P01J)/(P0*pi0J))*u01J*u00K
                               A0101= (P0*(1-P0)*P01J*P01K+I[1]*P0*P0*P01J*(1-P01J)/(P0*pi0J))*u01J*u01K
                               A0110=-(P0*(1-P0)*P01J*P10K)*u01J*u10K
                               A0111=-(P0*(1-P0)*P01J*P11K)*u01J*u11K

                               A1000=-(P1*(1-P1)*P10J*P00K)*u10J*u10K
                               A1001=-(P1*(1-P1)*P10J*P01K)*u10J*u11K
                               A1010= (P1*(1-P1)*P10J*P10K+I[1]*P1*P1*P10J*(1-P10J)/(P1*pi1J))*u10J*u10K
                               A1011= (P1*(1-P1)*P10J*P11K-I[1]*P1*P1*P10J*(1-P10J)/(P1*pi1J))*u10J*u11K

                               A1100=-(P1*(1-P1)*P10J*P00K)*u10J*u10K
                               A1101=-(P1*(1-P1)*P10J*P01K)*u10J*u11K
                               A1110= (P1*(1-P1)*P11J*P10K+I[1]*P1*P1*P11J*(1-P11J)/(P1*pi1J))*u11J*u10K
                               A1111= (P1*(1-P1)*P10J*P11K-I[1]*P1*P1*P11J*(1-P11J)/(P1*pi1J))*u11J*u11K
                               A=sum(A0000,A0001,A0010,A0011,
                                     A0100,A0101,A0110,A0111,
                                     A1000,A1001,A1010,A1011,
                                     A1100,A1101,A1110,A1111)

                               if (I[1]==1 & I[3]==1) {B0000=P0*P0*P00J*P00J*v00J/(P0*pi0J*P00J*pi00J)} else {B0000=0}
                               if (I[1]==1 & I[4]==1) {B0101=P0*P0*P01J*P01J*v01J/(P0*pi0J*P01J*pi01J)} else {B0101=0}
                               if (I[2]==1 & I[5]==1) {B1010=P1*P1*P10J*P10J*v10J/(P1*pi1J*P10J*pi10J)} else {B1010=0}
                               if (I[2]==1 & I[6]==1) {B1111=P1*P1*P11J*P11J*v11J/(P1*pi1J*P11J*pi11J)} else {B1111=0}
                               B=sum(B0000,B0101,B1010,B1111)
                               tvc=A+B
                               } else
     if (Nstage==3 & Base==0) {piJ=mean(Smat$PI1[which(Smat$A1==ATS1[1])]); if (is.nan(piJ)) {piJ=0} # A1=dJ[1]

                               P0J=mean(Smat$P2[which(Smat$A1==ATS1[1] & Smat$O2==0)]); if (is.nan(P0J)) {P0J=0}
                               P1J=mean(Smat$P2[which(Smat$A1==ATS1[1] & Smat$O2==1)]); if (is.nan(P1J)) {P1J=0}
                               P0K=mean(Smat$P2[which(Smat$A1==ATS2[1] & Smat$O2==0)]); if (is.nan(P0K)) {P0K=0}
                               P1K=mean(Smat$P2[which(Smat$A1==ATS2[1] & Smat$O2==1)]); if (is.nan(P1K)) {P1K=0}

                               pi0J=mean(Smat$PI2[which(Smat$A1==ATS1[1] & Smat$O2==0 & Smat$A2==ATS1[2])]); if (is.nan(pi0J)) {pi0J=0}
                               pi1J=mean(Smat$PI2[which(Smat$A1==ATS1[1] & Smat$O2==1 & Smat$A2==ATS1[3])]); if (is.nan(pi1J)) {pi1J=0}

                               P00J=mean(Smat$P3[which(Smat$A1==ATS1[1] & Smat$O2==0 & Smat$A2==ATS1[2] & Smat$O3==0)]); if (is.nan(P00J)) {P00J=0}
                               P01J=mean(Smat$P3[which(Smat$A1==ATS1[1] & Smat$O2==0 & Smat$A2==ATS1[2] & Smat$O3==1)]); if (is.nan(P01J)) {P01J=0}
                               P10J=mean(Smat$P3[which(Smat$A1==ATS1[1] & Smat$O2==1 & Smat$A2==ATS1[3] & Smat$O3==0)]); if (is.nan(P10J)) {P10J=0}
                               P11J=mean(Smat$P3[which(Smat$A1==ATS1[1] & Smat$O2==1 & Smat$A2==ATS1[3] & Smat$O3==1)]); if (is.nan(P11J)) {P11J=0}

                               P00K=mean(Smat$P3[which(Smat$A1==ATS2[1] & Smat$O2==0 & Smat$A2==ATS2[2] & Smat$O3==0)]); if (is.nan(P00K)) {P00K=0}
                               P01K=mean(Smat$P3[which(Smat$A1==ATS2[1] & Smat$O2==0 & Smat$A2==ATS2[2] & Smat$O3==1)]); if (is.nan(P01K)) {P01K=0}
                               P10K=mean(Smat$P3[which(Smat$A1==ATS2[1] & Smat$O2==1 & Smat$A2==ATS2[3] & Smat$O3==0)]); if (is.nan(P10K)) {P10K=0}
                               P11K=mean(Smat$P3[which(Smat$A1==ATS2[1] & Smat$O2==1 & Smat$A2==ATS2[3] & Smat$O3==1)]); if (is.nan(P11K)) {P11K=0}

                               pi00J=mean(Smat$PI3[which(Smat$A1==ATS1[1] & Smat$O2==0 & Smat$A2==ATS1[2] & Smat$O3==0 & Smat$A3==ATS1[4])]); if (is.nan(pi00J)) {pi00J=0}
                               pi01J=mean(Smat$PI3[which(Smat$A1==ATS1[1] & Smat$O2==0 & Smat$A2==ATS1[2] & Smat$O3==1 & Smat$A3==ATS1[5])]); if (is.nan(pi01J)) {pi01J=0}
                               pi10J=mean(Smat$PI3[which(Smat$A1==ATS1[1] & Smat$O2==1 & Smat$A2==ATS1[3] & Smat$O3==0 & Smat$A3==ATS1[6])]); if (is.nan(pi10J)) {pi10J=0}
                               pi11J=mean(Smat$PI3[which(Smat$A1==ATS1[1] & Smat$O2==1 & Smat$A2==ATS1[3] & Smat$O3==1 & Smat$A3==ATS1[7])]); if (is.nan(pi11J)) {pi11J=0}

                               u00J=Smat$MEAN[which(Smat$A1==ATS1[1] & Smat$O2==0 & Smat$A2==ATS1[2] & Smat$O3==0 & Smat$A3==ATS1[4])]; if (is.nan(u00J)) {u00J=0}        # within-sequence true mean (EY)
                               u01J=Smat$MEAN[which(Smat$A1==ATS1[1] & Smat$O2==0 & Smat$A2==ATS1[2] & Smat$O3==0 & Smat$A3==ATS1[5])]; if (is.nan(u01J)) {u01J=0}
                               u10J=Smat$MEAN[which(Smat$A1==ATS1[1] & Smat$O2==0 & Smat$A2==ATS1[3] & Smat$O3==0 & Smat$A3==ATS1[6])]; if (is.nan(u10J)) {u10J=0}
                               u11J=Smat$MEAN[which(Smat$A1==ATS1[1] & Smat$O2==0 & Smat$A2==ATS1[3] & Smat$O3==0 & Smat$A3==ATS1[7])]; if (is.nan(u11J)) {u11J=0}

                               u00K=Smat$MEAN[which(Smat$A1==ATS2[1] & Smat$O2==0 & Smat$A2==ATS2[2] & Smat$O3==0 & Smat$A3==ATS2[4])]; if (is.nan(u00K)) {u00K=0}        # within-sequence true mean (EY)
                               u01K=Smat$MEAN[which(Smat$A1==ATS2[1] & Smat$O2==0 & Smat$A2==ATS2[2] & Smat$O3==0 & Smat$A3==ATS2[5])]; if (is.nan(u01K)) {u01K=0}
                               u10K=Smat$MEAN[which(Smat$A1==ATS2[1] & Smat$O2==0 & Smat$A2==ATS2[3] & Smat$O3==0 & Smat$A3==ATS2[6])]; if (is.nan(u10K)) {u10K=0}
                               u11K=Smat$MEAN[which(Smat$A1==ATS2[1] & Smat$O2==0 & Smat$A2==ATS2[3] & Smat$O3==0 & Smat$A3==ATS2[7])]; if (is.nan(u11K)) {u11K=0}

                               if (FA=="normal"){v00J=Smat$SD[which(Smat$A1==ATS1[1] & Smat$O2==0 & Smat$A2==ATS1[2] & Smat$O3==0 & Smat$A3==ATS1[4])]^2        # within-sequence true mean (EY)
                                                 v01J=Smat$SD[which(Smat$A1==ATS1[1] & Smat$O2==0 & Smat$A2==ATS1[2] & Smat$O3==0 & Smat$A3==ATS1[5])]^2
                                                 v10J=Smat$SD[which(Smat$A1==ATS1[1] & Smat$O2==0 & Smat$A2==ATS1[3] & Smat$O3==0 & Smat$A3==ATS1[6])]^2
                                                 v11J=Smat$SD[which(Smat$A1==ATS1[1] & Smat$O2==0 & Smat$A2==ATS1[3] & Smat$O3==0 & Smat$A3==ATS1[7])]^2}
                               if (FA=="binomial") {bv00J=u00J*(1-u00J);  if (is.na(bv00J)) {bv00J=0}                              #within-sequence var est for bin outcome
                                                    bv01J=u01J*(1-u01J);  if (is.na(bv01J)) {bv01J=0}
                                                    bv10J=u10J*(1-u10J);  if (is.na(bv10J)) {bv10J=0}
                                                    bv11J=u11J*(1-u11J);  if (is.na(bv11J)) {bv11J=0}}
                               if (length(v00J)==0L) {v00J=0}
                               if (length(v01J)==0L) {v01J=0}
                               if (length(v10J)==0L) {v10J=0}
                               if (length(v11J)==0L) {v11J=0}

                               A0000= (I[1]*P0J*(1-P0J)*P00J*P00K/piJ+I[2]*P0J*P0J*P00J*(1-P00J)/(piJ*P0J*pi00J))*u00J*u00K
                               A0001= (I[1]*P0J*(1-P0J)*P00J*P01K/piJ-I[2]*P0J*P0J*P00J*(1-P00J)/(piJ*P0J*pi00J))*u00J*u01K
                               A0010=-(I[1]*P0J*(1-P0J)*P00J*P10K/piJ)*u00J*u10K
                               A0011=-(I[1]*P0J*(1-P0J)*P00J*P11K/piJ)*u00J*u11K

                               A0100= (I[1]*P0J*(1-P0J)*P01J*P00K/piJ-I[2]*P0J*P0J*P01J*(1-P01J)/(piJ*P0J*pi01J))*u01J*u00K
                               A0101= (I[1]*P0J*(1-P0J)*P01J*P01K/piJ+I[2]*P0J*P0J*P01J*(1-P01J)/(piJ*P0J*pi01J))*u01J*u01K
                               A0110=-(I[1]*P0J*(1-P0J)*P01J*P10K/piJ)*u01J*u10K
                               A0111=-(I[1]*P0J*(1-P0J)*P01J*P11K/piJ)*u01J*u11K

                               A1000=-(I[1]*P1J*(1-P1J)*P10J*P00K/piJ)*u10J*u10K
                               A1001=-(I[1]*P1J*(1-P1J)*P10J*P01K/piJ)*u10J*u11K
                               A1010= (I[1]*P1J*(1-P1J)*P10J*P10K/piJ+I[3]*P1J*P1J*P10J*(1-P10J)/(piJ*P1J*pi10J))*u10J*u10K
                               A1011= (I[1]*P1J*(1-P1J)*P10J*P11K/piJ-I[3]*P1J*P1J*P10J*(1-P10J)/(piJ*P1J*pi10J))*u10J*u11K

                               A1100=-(I[1]*P1J*(1-P1J)*P10J*P00K/piJ)*u10J*u10K
                               A1101=-(I[1]*P1J*(1-P1J)*P10J*P01K/piJ)*u10J*u11K
                               A1110= (I[1]*P1J*(1-P1J)*P11J*P10K/piJ+I[3]*P1J*P1J*P11J*(1-P11J)/(piJ*P1J*pi11J))*u11J*u10K
                               A1111= (I[1]*P1J*(1-P1J)*P10J*P11K/piJ-I[3]*P1J*P1J*P11J*(1-P11J)/(piJ*P1J*pi11J))*u11J*u11K
                               A=sum(A0000,A0001,A0010,A0011,
                                     A0100,A0101,A0110,A0111,
                                     A1000,A1001,A1010,A1011,
                                     A1100,A1101,A1110,A1111)

                               if (FA=="binomial"){if (I[1]==1 & I[2]==1 & I[4]==1) {B0000=P0J*P0J*P00J*P00K*bv00J/(piJ*P0J*pi00J*P00J*pi00J)}  else {B0000=0}
                                                   if (I[1]==1 & I[2]==1 & I[5]==1) {B0101=P0J*P0J*P01J*P01K*bv01J/(piJ*P0J*pi00J*P01J*pi01J)}  else {B0101=0}
                                                   if (I[1]==1 & I[3]==1 & I[6]==1) {B1010=P1J*P1J*P10J*P10K*bv10J/(piJ*P1J*pi10J*P10J*pi10J)}  else {B1010=0}
                                                   if (I[1]==1 & I[3]==1 & I[7]==1) {B1111=P1J*P1J*P11J*P11K*bv11J/(piJ*P1J*pi11J*P11J*pi11J)}  else {B1111=0}
                                                   } else
                               if (FA=="normal")  {if (I[1]==1 & I[2]==1 & I[4]==1) {B0000=P0J*P0J*P00J*P00K*v00J/(piJ*P0J*pi00J*P00J*pi00J)}   else {B0000=0}
                                                   if (I[1]==1 & I[2]==1 & I[5]==1) {B0101=P0J*P0J*P01J*P01K*v01J/(piJ*P0J*pi00J*P01J*pi01J)}   else {B0101=0}
                                                   if (I[1]==1 & I[3]==1 & I[6]==1) {B1010=P1J*P1J*P10J*P10K*v10J/(piJ*P1J*pi10J*P10J*pi10J)}   else {B1010=0}
                                                   if (I[1]==1 & I[3]==1 & I[7]==1) {B1111=P1J*P1J*P11J*P11K*v11J/(piJ*P1J*pi11J*P11J*pi11J)}   else {B1111=0}
                                                   }
                               B=sum(B0000,B0101,B1010,B1111)
                               tvc=A+B
                               } else
     if (Nstage==3 & Base==1) {P0=mean(Smat$P1[which(Smat$O1==0)]); if (is.nan(P0)) {P0=0}
                               P1=mean(Smat$P1[which(Smat$O1==1)]); if (is.nan(P1)) {P1=0}

                               pi0J=mean(Smat$PI1[which(Smat$O1==0 & Smat$A1==ATS1[1])]); if (is.nan(pi0J)) {pi0J=1}
                               pi1J=mean(Smat$PI1[which(Smat$O1==1 & Smat$A1==ATS1[2])]); if (is.nan(pi1J)) {pi1J=1}

                               P00J=mean(Smat$P2[which(Smat$O1==0 & Smat$A1==ATS1[1] & Smat$O2==0)]); if (is.nan(P00J)) {P00J=0}
                               P01J=mean(Smat$P2[which(Smat$O1==0 & Smat$A1==ATS1[1] & Smat$O2==1)]); if (is.nan(P01J)) {P01J=0}
                               P10J=mean(Smat$P2[which(Smat$O1==1 & Smat$A1==ATS1[2] & Smat$O2==0)]); if (is.nan(P10J)) {P10J=0}
                               P11J=mean(Smat$P2[which(Smat$O1==1 & Smat$A1==ATS1[2] & Smat$O2==1)]); if (is.nan(P11J)) {P11J=0}

                               P00K=mean(Smat$P2[which(Smat$O1==0 & Smat$A1==ATS2[1] & Smat$O2==0)]); if (is.nan(P00K)) {P00K=0}
                               P01K=mean(Smat$P2[which(Smat$O1==0 & Smat$A1==ATS2[1] & Smat$O2==1)]); if (is.nan(P01K)) {P01K=0}
                               P10K=mean(Smat$P2[which(Smat$O1==1 & Smat$A1==ATS2[2] & Smat$O2==0)]); if (is.nan(P10K)) {P10K=0}
                               P11K=mean(Smat$P2[which(Smat$O1==1 & Smat$A1==ATS2[2] & Smat$O2==1)]); if (is.nan(P11K)) {P11K=0}

                               pi00J=mean(Smat$PI2[which(Smat$O1==0 & Smat$A1==ATS1[1] & Smat$O2==0 & Smat$A2==ATS1[3])]); if (is.nan(pi00J)) {pi00J=1}
                               pi01J=mean(Smat$PI2[which(Smat$O1==0 & Smat$A1==ATS1[1] & Smat$O2==1 & Smat$A2==ATS1[4])]); if (is.nan(pi01J)) {pi01J=1}
                               pi10J=mean(Smat$PI2[which(Smat$O1==1 & Smat$A1==ATS1[2] & Smat$O2==0 & Smat$A2==ATS1[5])]); if (is.nan(pi10J)) {pi10J=1}
                               pi11J=mean(Smat$PI2[which(Smat$O1==1 & Smat$A1==ATS1[2] & Smat$O2==1 & Smat$A2==ATS1[6])]); if (is.nan(pi11J)) {pi11J=1}

                               P000J=mean(Smat$P3[which(Smat$O1==0 & Smat$A1==ATS1[1] & Smat$O2==0 & Smat$A2==ATS1[3] & Smat$O3==0)]); if (is.nan(P000J)) {P000J=0}
                               P001J=mean(Smat$P3[which(Smat$O1==0 & Smat$A1==ATS1[1] & Smat$O2==0 & Smat$A2==ATS1[3] & Smat$O3==1)]); if (is.nan(P001J)) {P001J=0}
                               P010J=mean(Smat$P3[which(Smat$O1==0 & Smat$A1==ATS1[1] & Smat$O2==1 & Smat$A2==ATS1[4] & Smat$O3==0)]); if (is.nan(P010J)) {P010J=0}
                               P011J=mean(Smat$P3[which(Smat$O1==0 & Smat$A1==ATS1[1] & Smat$O2==1 & Smat$A2==ATS1[4] & Smat$O3==1)]); if (is.nan(P011J)) {P011J=0}
                               P100J=mean(Smat$P3[which(Smat$O1==1 & Smat$A1==ATS1[2] & Smat$O2==0 & Smat$A2==ATS1[5] & Smat$O3==0)]); if (is.nan(P100J)) {P100J=0}
                               P101J=mean(Smat$P3[which(Smat$O1==1 & Smat$A1==ATS1[2] & Smat$O2==0 & Smat$A2==ATS1[5] & Smat$O3==1)]); if (is.nan(P101J)) {P101J=0}
                               P110J=mean(Smat$P3[which(Smat$O1==1 & Smat$A1==ATS1[2] & Smat$O2==1 & Smat$A2==ATS1[6] & Smat$O3==0)]); if (is.nan(P110J)) {P110J=0}
                               P111J=mean(Smat$P3[which(Smat$O1==1 & Smat$A1==ATS1[2] & Smat$O2==1 & Smat$A2==ATS1[6] & Smat$O3==1)]); if (is.nan(P111J)) {P111J=0}

                               P000K=mean(Smat$P3[which(Smat$O1==0 & Smat$A1==ATS2[1] & Smat$O2==0 & Smat$A2==ATS2[3] & Smat$O3==0)]); if (is.nan(P000K)) {P000K=0}
                               P001K=mean(Smat$P3[which(Smat$O1==0 & Smat$A1==ATS2[1] & Smat$O2==0 & Smat$A2==ATS2[3] & Smat$O3==1)]); if (is.nan(P001K)) {P001K=0}
                               P010K=mean(Smat$P3[which(Smat$O1==0 & Smat$A1==ATS2[1] & Smat$O2==1 & Smat$A2==ATS2[4] & Smat$O3==0)]); if (is.nan(P010K)) {P010K=0}
                               P011K=mean(Smat$P3[which(Smat$O1==0 & Smat$A1==ATS2[1] & Smat$O2==1 & Smat$A2==ATS2[4] & Smat$O3==1)]); if (is.nan(P011K)) {P011K=0}
                               P100K=mean(Smat$P3[which(Smat$O1==1 & Smat$A1==ATS2[2] & Smat$O2==0 & Smat$A2==ATS2[5] & Smat$O3==0)]); if (is.nan(P100K)) {P100K=0}
                               P101K=mean(Smat$P3[which(Smat$O1==1 & Smat$A1==ATS2[2] & Smat$O2==0 & Smat$A2==ATS2[5] & Smat$O3==1)]); if (is.nan(P101K)) {P101K=0}
                               P110K=mean(Smat$P3[which(Smat$O1==1 & Smat$A1==ATS2[2] & Smat$O2==1 & Smat$A2==ATS2[6] & Smat$O3==0)]); if (is.nan(P110K)) {P110K=0}
                               P111K=mean(Smat$P3[which(Smat$O1==1 & Smat$A1==ATS2[2] & Smat$O2==1 & Smat$A2==ATS2[6] & Smat$O3==1)]); if (is.nan(P111K)) {P111K=0}

                               pi000J=mean(Smat$PI3[which(Smat$O1==0 & Smat$A1==ATS1[1] & Smat$O2==0 & Smat$A2==ATS1[3] & Smat$O3==0 & Smat$A3==ATS1[7])]); if (is.nan(P000J)) {P000J=0}
                               pi001J=mean(Smat$PI3[which(Smat$O1==0 & Smat$A1==ATS1[1] & Smat$O2==0 & Smat$A2==ATS1[3] & Smat$O3==1 & Smat$A3==ATS1[8])]); if (is.nan(P001J)) {P001J=0}
                               pi010J=mean(Smat$PI3[which(Smat$O1==0 & Smat$A1==ATS1[1] & Smat$O2==1 & Smat$A2==ATS1[4] & Smat$O3==0 & Smat$A3==ATS1[9])]); if (is.nan(P010J)) {P010J=0}
                               pi011J=mean(Smat$PI3[which(Smat$O1==0 & Smat$A1==ATS1[1] & Smat$O2==1 & Smat$A2==ATS1[4] & Smat$O3==1 & Smat$A3==ATS1[10])]);if (is.nan(P011J)) {P011J=0}
                               pi100J=mean(Smat$PI3[which(Smat$O1==1 & Smat$A1==ATS1[2] & Smat$O2==0 & Smat$A2==ATS1[5] & Smat$O3==0 & Smat$A3==ATS1[11])]);if (is.nan(P100J)) {P100J=0}
                               pi101J=mean(Smat$PI3[which(Smat$O1==1 & Smat$A1==ATS1[2] & Smat$O2==0 & Smat$A2==ATS1[5] & Smat$O3==1 & Smat$A3==ATS1[12])]);if (is.nan(P101J)) {P101J=0}
                               pi110J=mean(Smat$PI3[which(Smat$O1==1 & Smat$A1==ATS1[2] & Smat$O2==1 & Smat$A2==ATS1[6] & Smat$O3==0 & Smat$A3==ATS1[13])]);if (is.nan(P110J)) {P110J=0}
                               pi111J=mean(Smat$PI3[which(Smat$O1==1 & Smat$A1==ATS1[2] & Smat$O2==1 & Smat$A2==ATS1[6] & Smat$O3==1 & Smat$A3==ATS1[14])]);if (is.nan(P111J)) {P111J=0}

                               u000J=Smat$MEAN[which(Smat$O1==0 & Smat$A1==ATS1[1] & Smat$O2==0 & Smat$A2==ATS1[3] & Smat$O3==0 & Smat$A3==ATS1[7])];  if (is.na(u000J)) {u000J=0}
                               u001J=Smat$MEAN[which(Smat$O1==0 & Smat$A1==ATS1[1] & Smat$O2==0 & Smat$A2==ATS1[3] & Smat$O3==1 & Smat$A3==ATS1[8])];  if (is.na(u001J)) {u001J=0}
                               u010J=Smat$MEAN[which(Smat$O1==0 & Smat$A1==ATS1[1] & Smat$O2==1 & Smat$A2==ATS1[4] & Smat$O3==0 & Smat$A3==ATS1[9])];  if (is.na(u010J)) {u010J=0}
                               u011J=Smat$MEAN[which(Smat$O1==0 & Smat$A1==ATS1[1] & Smat$O2==1 & Smat$A2==ATS1[4] & Smat$O3==1 & Smat$A3==ATS1[10])]; if (is.na(u011J)) {u011J=0}
                               u100J=Smat$MEAN[which(Smat$O1==1 & Smat$A1==ATS1[2] & Smat$O2==0 & Smat$A2==ATS1[5] & Smat$O3==0 & Smat$A3==ATS1[11])]; if (is.na(u100J)) {u100J=0}
                               u101J=Smat$MEAN[which(Smat$O1==1 & Smat$A1==ATS1[2] & Smat$O2==0 & Smat$A2==ATS1[5] & Smat$O3==1 & Smat$A3==ATS1[12])]; if (is.na(u101J)) {u101J=0}
                               u110J=Smat$MEAN[which(Smat$O1==1 & Smat$A1==ATS1[2] & Smat$O2==1 & Smat$A2==ATS1[6] & Smat$O3==0 & Smat$A3==ATS1[13])]; if (is.na(u110J)) {u110J=0}
                               u111J=Smat$MEAN[which(Smat$O1==1 & Smat$A1==ATS1[2] & Smat$O2==1 & Smat$A2==ATS1[6] & Smat$O3==1 & Smat$A3==ATS1[14])]; if (is.na(u111J)) {u111J=0}

                               u000K=Smat$MEAN[which(Smat$O1==0 & Smat$A1==ATS2[1] & Smat$O2==0 & Smat$A2==ATS2[3] & Smat$O3==0 & Smat$A3==ATS2[7])];  if (is.na(u000K)) {u000K=0}
                               u001K=Smat$MEAN[which(Smat$O1==0 & Smat$A1==ATS2[1] & Smat$O2==0 & Smat$A2==ATS2[3] & Smat$O3==1 & Smat$A3==ATS2[8])];  if (is.na(u001K)) {u001K=0}
                               u010K=Smat$MEAN[which(Smat$O1==0 & Smat$A1==ATS2[1] & Smat$O2==1 & Smat$A2==ATS2[4] & Smat$O3==0 & Smat$A3==ATS2[9])];  if (is.na(u010K)) {u010K=0}
                               u011K=Smat$MEAN[which(Smat$O1==0 & Smat$A1==ATS2[1] & Smat$O2==1 & Smat$A2==ATS2[4] & Smat$O3==1 & Smat$A3==ATS2[10])]; if (is.na(u011K)) {u011K=0}
                               u100K=Smat$MEAN[which(Smat$O1==1 & Smat$A1==ATS2[2] & Smat$O2==0 & Smat$A2==ATS2[5] & Smat$O3==0 & Smat$A3==ATS2[11])]; if (is.na(u100K)) {u100K=0}
                               u101K=Smat$MEAN[which(Smat$O1==1 & Smat$A1==ATS2[2] & Smat$O2==0 & Smat$A2==ATS2[5] & Smat$O3==1 & Smat$A3==ATS2[12])]; if (is.na(u101K)) {u101K=0}
                               u110K=Smat$MEAN[which(Smat$O1==1 & Smat$A1==ATS2[2] & Smat$O2==1 & Smat$A2==ATS2[6] & Smat$O3==0 & Smat$A3==ATS2[13])]; if (is.na(u110K)) {u110K=0}
                               u111K=Smat$MEAN[which(Smat$O1==1 & Smat$A1==ATS2[2] & Smat$O2==1 & Smat$A2==ATS2[6] & Smat$O3==1 & Smat$A3==ATS2[14])]; if (is.na(u111K)) {u111K=0}

                               if (FA=="normal"){v000J=Smat$SD[which(Smat$O1==0 & Smat$A1==ATS1[1] & Smat$O2==0 & Smat$A2==ATS1[3] & Smat$O3==0 & Smat$A3==ATS1[7])]^2
                                                 v001J=Smat$SD[which(Smat$O1==0 & Smat$A1==ATS1[1] & Smat$O2==0 & Smat$A2==ATS1[3] & Smat$O3==1 & Smat$A3==ATS1[8])]^2
                                                 v010J=Smat$SD[which(Smat$O1==0 & Smat$A1==ATS1[1] & Smat$O2==1 & Smat$A2==ATS1[4] & Smat$O3==0 & Smat$A3==ATS1[9])]^2
                                                 v011J=Smat$SD[which(Smat$O1==0 & Smat$A1==ATS1[1] & Smat$O2==1 & Smat$A2==ATS1[4] & Smat$O3==1 & Smat$A3==ATS1[10])]^2
                                                 v100J=Smat$SD[which(Smat$O1==1 & Smat$A1==ATS1[2] & Smat$O2==0 & Smat$A2==ATS1[5] & Smat$O3==0 & Smat$A3==ATS1[11])]^2
                                                 v101J=Smat$SD[which(Smat$O1==1 & Smat$A1==ATS1[2] & Smat$O2==0 & Smat$A2==ATS1[5] & Smat$O3==1 & Smat$A3==ATS1[12])]^2
                                                 v110J=Smat$SD[which(Smat$O1==1 & Smat$A1==ATS1[2] & Smat$O2==1 & Smat$A2==ATS1[6] & Smat$O3==0 & Smat$A3==ATS1[13])]^2
                                                 v111J=Smat$SD[which(Smat$O1==1 & Smat$A1==ATS1[2] & Smat$O2==1 & Smat$A2==ATS1[6] & Smat$O3==1 & Smat$A3==ATS1[14])]^2} else
                               if (FA=="binomial") {bv000J=u000J*(1-u000J)
                                                    bv001J=u001J*(1-u001J)
                                                    bv010J=u010J*(1-u010J)
                                                    bv011J=u011J*(1-u011J)
                                                    bv100J=u100J*(1-u100J)
                                                    bv101J=u101J*(1-u101J)
                                                    bv110J=u110J*(1-u110J)
                                                    bv111J=u111J*(1-u111J)}
                               if (length(v000J)==0L) {v000J=0}
                               if (length(v001J)==0L) {v001J=0}
                               if (length(v010J)==0L) {v010J=0}
                               if (length(v011J)==0L) {v011J=0}
                               if (length(v100J)==0L) {v100J=0}
                               if (length(v101J)==0L) {v101J=0}
                               if (length(v110J)==0L) {v110J=0}
                               if (length(v111J)==0L) {v111J=0}

                               A000.000=( P0*(1-P0)*P00J*P00K*P000J*P000K+I[1]*P0*P0*P00J*(1-P00J)*P000J*P000K/(P0*pi0J)+I[1]*I[3]*P0*P0*P00J*P00K*P000J*(1-P000J)/(P0*pi0J*P00J*pi00J))*u000J*u000K
                               A000.001=( P0*(1-P0)*P00J*P00K*P000J*P001K+I[1]*P0*P0*P00J*(1-P00J)*P000J*P001K/(P0*pi0J)-I[1]*I[3]*P0*P0*P00J*P00K*P000J*(1-P000J)/(P0*pi0J*P00J*pi00J))*u000J*u001K
                               A000.010=( P0*(1-P0)*P00J*P10K*P000J*P010K-I[1]*P0*P0*P00J*(1-P00J)*P000J*P010K/(P0*pi0J))*u000J*u010K
                               A000.011=( P0*(1-P0)*P00J*P01K*P000J*P011K-I[1]*P0*P0*P00J*(1-P00J)*P000J*P011K/(P0*pi0J))*u000J*u011K
                               A000.100=(-P0*(1-P0)*P00J*P10K*P000J*P100K)*u000J*u100K
                               A000.101=(-P0*(1-P0)*P00J*P10K*P000J*P101K)*u000J*u101K
                               A000.110=(-P0*(1-P0)*P00J*P11K*P000J*P110K)*u000J*u110K
                               A000.111=(-P0*(1-P0)*P00J*P11K*P000J*P111K)*u000J*u111K

                               A001.000=( P0*(1-P0)*P00J*P00K*P001J*P000K+I[1]*P0*P0*P00J*(1-P00J)*P001J*P000K/(P0*pi0J)-I[1]*I[3]*P0*P0*P00J*P00K*P001J*(1-P001J)/(P0*pi0J*P00J*pi00J))*u001J*u000K
                               A001.001=( P0*(1-P0)*P00J*P00K*P001J*P001K+I[1]*P0*P0*P00J*(1-P00J)*P001J*P001K/(P0*pi0J)+I[1]*I[3]*P0*P0*P00J*P00K*P001J*(1-P001J)/(P0*pi0J*P00J*pi00J))*u001J*u001K
                               A001.010=( P0*(1-P0)*P00J*P10K*P001J*P010K-I[1]*P0*P0*P00J*(1-P00J)*P001J*P010K/(P0*pi0J))*u001J*u010K
                               A001.011=( P0*(1-P0)*P00J*P10K*P001J*P011K-I[1]*P0*P0*P00J*(1-P00J)*P001J*P011K/(P0*pi0J))*u001J*u011K
                               A001.100=(-P0*(1-P0)*P00J*P10K*P001J*P100K)*u001J*u100K
                               A001.101=(-P0*(1-P0)*P00J*P10K*P001J*P101K)*u001J*u101K
                               A001.110=(-P0*(1-P0)*P00J*P11K*P001J*P110K)*u001J*u110K
                               A001.111=(-P0*(1-P0)*P00J*P11K*P001J*P111K)*u001J*u111K

                               A010.000=( P0*(1-P0)*P01J*P00K*P010J*P000K-I[1]*P0*P0*P01J*(1-P01J)*P010J*P000K/(P0*pi0J))*u010J*u000K
                               A010.001=( P0*(1-P0)*P01J*P00K*P010J*P001K-I[1]*P0*P0*P01J*(1-P01J)*P010J*P001K/(P0*pi0J))*u010J*u001K
                               A010.010=( P0*(1-P0)*P01J*P01K*P010J*P010K+I[1]*P0*P0*P01J*(1-P01J)*P010J*P010K/(P0*pi0J)+I[1]*I[3]*P0*P0*P01J*P01K*P010J*(1-P010J)/(P0*pi0J*P01J*pi01J))*u010J*u010K
                               A010.011=( P0*(1-P0)*P01J*P01K*P010J*P011K+I[1]*P0*P0*P01J*(1-P01J)*P010J*P011K/(P0*pi0J)-I[1]*I[3]*P0*P0*P01J*P01K*P010J*(1-P010J)/(P0*pi0J*P01J*pi01J))*u010J*u011K
                               A010.100=(-P0*(1-P0)*P01J*P10K*P010J*P100K)*u010J*u100K
                               A010.101=(-P0*(1-P0)*P01J*P10K*P010J*P101K)*u010J*u101K
                               A010.110=(-P0*(1-P0)*P01J*P11K*P010J*P110K)*u010J*u110K
                               A010.111=(-P0*(1-P0)*P01J*P11K*P010J*P111K)*u010J*u111K

                               A011.000=( P0*(1-P0)*P01J*P00K*P011J*P000K-I[1]*P0*P0*P01J*(1-P01J)*P011J*P000K/(P0*pi0J))*u011J*u000K
                               A011.001=( P0*(1-P0)*P01J*P00K*P011J*P001K-I[1]*P0*P0*P01J*(1-P01J)*P011J*P001K/(P0*pi0J))*u011J*u001K
                               A011.010=( P0*(1-P0)*P01J*P01K*P011J*P010K+I[1]*P0*P0*P01J*(1-P01J)*P011J*P010K/(P0*pi0J)-I[1]*I[3]*P0*P0*P01J*P01K*P011J*(1-P011J)/(P0*pi0J*P01J*pi01J))*u011J*u010K
                               A011.011=( P0*(1-P0)*P01J*P01K*P011J*P011K+I[1]*P0*P0*P01J*(1-P01J)*P011J*P011K/(P0*pi0J)+I[1]*I[3]*P0*P0*P01J*P01K*P011J*(1-P011J)/(P0*pi0J*P01J*pi01J))*u011J*u011K
                               A011.100=(-P0*(1-P0)*P01J*P10K*P011J*P100K)*u011J*u100K
                               A011.101=(-P0*(1-P0)*P01J*P10K*P011J*P101K)*u011J*u101K
                               A011.110=(-P0*(1-P0)*P01J*P11K*P011J*P110K)*u011J*u110K
                               A011.111=(-P0*(1-P0)*P01J*P11K*P011J*P111K)*u011J*u111K

                               A100.000=(-P1*(1-P1)*P10J*P00K*P100J*P000K)*u100J*u000K
                               A100.001=(-P1*(1-P1)*P10J*P00K*P100J*P001K)*u100J*u001K
                               A100.010=(-P1*(1-P1)*P10J*P01K*P100J*P010K)*u100J*u010K
                               A100.011=(-P1*(1-P1)*P10J*P01K*P100J*P011K)*u100J*u011K
                               A100.100=( P1*(1-P1)*P10J*P10K*P100J*P100K+I[1]*P1*P1*P10J*(1-P10J)*P100J*P100K/(P1*pi1J)+I[1]*I[3]*P1*P1*P10J*P10K*P100J*(1-P100J)/(P1*pi1J*P10J*pi10J))*u100J*u100K
                               A100.101=( P1*(1-P1)*P10J*P10K*P100J*P101K+I[1]*P1*P1*P10J*(1-P10J)*P100J*P101K/(P1*pi1J)-I[1]*I[3]*P1*P1*P10J*P10K*P100J*(1-P100J)/(P1*pi1J*P10J*pi10J))*u100J*u101K
                               A100.110=( P1*(1-P1)*P10J*P11K*P100J*P110K-I[1]*P1*P1*P10J*(1-P10J)*P100J*P110K/(P1*pi1J))*u100J*u110K
                               A100.111=( P1*(1-P1)*P10J*P11K*P100J*P111K-I[1]*P1*P1*P10J*(1-P10J)*P100J*P111K/(P1*pi1J))*u100J*u111K

                               A101.000=(-P1*(1-P1)*P10J*P00K*P101J*P000K)*u101J*u000K
                               A101.001=(-P1*(1-P1)*P10J*P00K*P101J*P001K)*u101J*u001K
                               A101.010=(-P1*(1-P1)*P10J*P01K*P101J*P010K)*u101J*u010K
                               A101.011=(-P1*(1-P1)*P10J*P01K*P101J*P011K)*u101J*u011K
                               A101.100=( P1*(1-P1)*P10J*P10K*P101J*P100K+I[1]*P1*P1*P10J*(1-P10J)*P101J*P100K/(P1*pi1J)-I[1]*I[3]*P1*P1*P10J*P10K*P101J*(1-P101J)/(P1*pi1J*P10J*pi10J))*u101J*u100K
                               A101.101=( P1*(1-P1)*P10J*P10K*P101J*P101K+I[1]*P1*P1*P10J*(1-P10J)*P101J*P101K/(P1*pi1J)+I[1]*I[3]*P1*P1*P10J*P10K*P101J*(1-P101J)/(P1*pi1J*P10J*pi10J))*u101J*u101K
                               A101.110=( P1*(1-P1)*P10J*P11K*P101J*P110K+I[1]*P1*P1*P10J*(1-P10J)*P101J*P110K/(P1*pi1J))*u101J*u110K
                               A101.111=( P1*(1-P1)*P10J*P11K*P101J*P111K+I[1]*P1*P1*P10J*(1-P10J)*P101J*P111K/(P1*pi1J))*u101J*u111K

                               A110.000=(-P1*(1-P1)*P11J*P00K*P110J*P000K)*u110J*u000K
                               A110.001=(-P1*(1-P1)*P11J*P00K*P110J*P001K)*u110J*u001K
                               A110.010=(-P1*(1-P1)*P11J*P01K*P110J*P010K)*u110J*u010K
                               A110.011=(-P1*(1-P1)*P11J*P01K*P110J*P011K)*u110J*u011K
                               A110.100=( P1*(1-P1)*P11J*P10K*P110J*P100K-I[1]*P1*P1*P11J*(1-P11J)*P110J*P100K/(P1*pi1J))*u110J*u100K
                               A110.101=( P1*(1-P1)*P11J*P10K*P110J*P101K-I[1]*P1*P1*P11J*(1-P11J)*P110J*P101K/(P1*pi1J))*u110J*u101K
                               A110.110=( P1*(1-P1)*P11J*P11K*P110J*P110K+I[1]*P1*P1*P11J*(1-P11J)*P110J*P110K/(P1*pi1J)-I[1]*I[3]*P1*P1*P11J*P11K*P110J*(1-P110J)/(P1*pi1J*P11J*pi11J))*u110J*u110K
                               A110.111=( P1*(1-P1)*P11J*P11K*P110J*P111K+I[1]*P1*P1*P11J*(1-P11J)*P111J*P111K/(P1*pi1J)-I[1]*I[3]*P1*P1*P11J*P11K*P110J*(1-P110J)/(P1*pi1J*P11J*pi11J))*u110J*u111K

                               A111.000=(-P1*(1-P1)*P11J*P00K*P111J*P000K)*u111J*u000K
                               A111.001=(-P1*(1-P1)*P11J*P00K*P111J*P000K)*u111J*u001K
                               A111.010=(-P1*(1-P1)*P11J*P00K*P111J*P000K)*u111J*u010K
                               A111.011=(-P1*(1-P1)*P11J*P00K*P111J*P000K)*u111J*u011K
                               A111.100=( P1*(1-P1)*P11J*P10K*P111J*P100K-I[1]*P1*P1*P11J*(1-P11J)*P111J*P100K/(P1*pi1J))*u111J*u100K
                               A111.101=( P1*(1-P1)*P11J*P10K*P111J*P101K-I[1]*P1*P1*P11J*(1-P11J)*P111J*P101K/(P1*pi1J))*u111J*u101K
                               A111.110=( P1*(1-P1)*P11J*P11K*P111J*P110K+I[1]*P1*P1*P11J*(1-P11J)*P111J*P110K/(P1*pi1J)-I[1]*I[3]*P1*P1*P11J*P11K*P111J*(1-P111J)/(P1*pi1J*P11J*pi11J))*u111J*u110K
                               A111.111=( P1*(1-P1)*P11J*P11K*P111J*P111K+I[1]*P1*P1*P11J*(1-P11J)*P111J*P111K/(P1*pi1J)+I[1]*I[3]*P1*P1*P11J*P11K*P111J*(1-P111J)/(P1*pi1J*P11J*pi11J))*u111J*u111K

                               A=sum(A000.000,A000.001,A000.010,A000.011,A000.100,A000.101,A000.110,A000.111,
                                     A001.000,A001.001,A001.010,A001.011,A001.100,A001.101,A001.110,A001.111,
                                     A010.000,A010.001,A010.010,A010.011,A010.100,A010.101,A010.110,A010.111,
                                     A011.000,A011.001,A011.010,A011.011,A011.100,A011.101,A011.110,A011.111,
                                     A100.000,A100.001,A100.010,A100.011,A100.100,A100.101,A100.110,A100.111,
                                     A101.000,A101.001,A101.010,A101.011,A101.100,A101.101,A101.110,A101.111,
                                     A110.000,A110.001,A110.010,A110.011,A110.100,A110.101,A110.110,A110.111,
                                     A111.000,A111.001,A111.010,A111.011,A111.100,A111.101,A111.110,A111.111)

                               if (FA=="binomial"){if (I[1]==1 & I[3]==1 & I[ 7]==1) {B000.000=P0*P0*P00J*P00K*P000J*P000K*bv000J/(P0*pi0J*P00J*pi00J*P000J*pi000J)}  else {B000.000=0}
                                                   if (I[1]==1 & I[3]==1 & I[ 8]==1) {B001.001=P0*P0*P01J*P01K*P001J*P001K*bv001J/(P0*pi0J*P00J*pi00J*P001J*pi001J)}  else {B001.001=0}
                                                   if (I[1]==1 & I[4]==1 & I[ 9]==1) {B010.010=P0*P0*P01J*P01K*P010J*P010K*bv010J/(P0*pi0J*P01J*pi01J*P010J*pi010J)}  else {B010.010=0}
                                                   if (I[1]==1 & I[4]==1 & I[10]==1) {B011.011=P0*P0*P01J*P01K*P011J*P011K*bv011J/(P0*pi0J*P01J*pi01J*P011J*pi011J)}  else {B011.011=0}
                                                   if (I[1]==1 & I[5]==1 & I[11]==1) {B100.100=P1*P1*P10J*P10K*P100J*P100K*bv100J/(P1*pi1J*P10J*pi10J*P100J*pi100J)}  else {B100.100=0}
                                                   if (I[1]==1 & I[5]==1 & I[12]==1) {B101.101=P1*P1*P10J*P10K*P101J*P101K*bv101J/(P1*pi1J*P10J*pi10J*P101J*pi101J)}  else {B101.101=0}
                                                   if (I[1]==1 & I[6]==1 & I[13]==1) {B110.110=P1*P1*P11J*P11K*P110J*P110K*bv110J/(P1*pi1J*P11J*pi11J*P110J*pi110J)}  else {B110.110=0}
                                                   if (I[1]==1 & I[6]==1 & I[14]==1) {B111.111=P1*P1*P11J*P11K*P111J*P111K*bv111J/(P1*pi1J*P11J*pi11J*P111J*pi111J)}  else {B111.111=0}
                                                   }
                                             else {if (I[1]==1 & I[3]==1 & I[ 7]==1) {B000.000=P0*P0*P00J*P00K*P000J*P000K*v000J/(P0*pi0J*P00J*pi00J*P000J*pi000J)}  else {B000.000=0}
                                                   if (I[1]==1 & I[3]==1 & I[ 8]==1) {B001.001=P0*P0*P01J*P01K*P001J*P001K*v001J/(P0*pi0J*P00J*pi00J*P001J*pi001J)}  else {B001.001=0}
                                                   if (I[1]==1 & I[4]==1 & I[ 9]==1) {B010.010=P0*P0*P01J*P01K*P010J*P010K*v010J/(P0*pi0J*P01J*pi01J*P010J*pi010J)}  else {B010.010=0}
                                                   if (I[1]==1 & I[4]==1 & I[10]==1) {B011.011=P0*P0*P01J*P01K*P011J*P011K*v011J/(P0*pi0J*P01J*pi01J*P011J*pi011J)}  else {B011.011=0}
                                                   if (I[1]==1 & I[5]==1 & I[11]==1) {B100.100=P1*P1*P10J*P10K*P100J*P100K*v100J/(P1*pi1J*P10J*pi10J*P100J*pi100J)}  else {B100.100=0}
                                                   if (I[1]==1 & I[5]==1 & I[12]==1) {B101.101=P1*P1*P10J*P10K*P101J*P101K*v101J/(P1*pi1J*P10J*pi10J*P101J*pi101J)}  else {B101.101=0}
                                                   if (I[1]==1 & I[6]==1 & I[13]==1) {B110.110=P1*P1*P11J*P11K*P110J*P110K*v110J/(P1*pi1J*P11J*pi11J*P110J*pi110J)}  else {B110.110=0}
                                                   if (I[1]==1 & I[6]==1 & I[14]==1) {B111.111=P1*P1*P11J*P11K*P111J*P111K*v111J/(P1*pi1J*P11J*pi11J*P111J*pi111J)}  else {B111.111=0}
                                                   }
                               B=sum(B000.000,B001.001,B010.010,B011.011,
                                     B100.100,B101.101,B110.110,B111.111)
                               tvc=A+B
                               }
     return(tvc)
}
