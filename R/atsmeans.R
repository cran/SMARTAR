#' Identify adaptive treatment strategy and estimate strategy values
#'
#' Return a message that lists all the adaptive treatment strategy embedded in SMART design. It also gives the estiamted strategy values and the variance-covariance matrix of estimated values.
#' @param data Input data frame of the sequential randomized trial (SMART) data used for analysis. The data should include the variables of stage-specific treatments (At; t=1,2,3), intermediate evaluation (Ot; t=1,2,3) and final primary outcome (Y), where t represent the number of stages embedded in design. If stage-1 treatment (A1) takes into account the information of baseline evaluation, O1 needed to be include in data, otherwise not.
#' @param family A character string to specify the type of final primary outcome. The default is family=“normal”, which refers to the continuous primary outcome. If family=”binomial” then the primary outcome will be treated as binary variable.
#' @param method A character string to specify the method of estimation. If method=``Gest'' then G-computation method is used. If method=``IPW'' then Inversed Probability Weighting method is used.
#' @param common If common=TRUE, the pooled variance across all the treatment sequences are used in estimation. Otherwise use the sequence-specific variance. The default is common=FALSE.
#' @param conf If conf=TRUE, output confidence intervals for estimate strategy values. The default is conf=FALSE.
#' @param alpha Type I error rate control for confidence interval. The default is alpha=0.05.
#' @param plot If plot=TRUE, output the graphs of treatment effects with CIs.
#' @return
#'
#' An object of ``value” is return, which contain the index of all the adaptive treatment strategies, strategy-specific sample sizes and estimated values with standardized errors.
##' \itemize{
##'    \item ATS: Index of adaptive treatment strategy from 1 to G, where G is total number of strategies defined in SMART
##'    \item ds: Stage-specific decision makings given certain histories corresponding to each strategy. The number of columns of ``ds'' is defined by strategy and details are shown in the output.
##'    \item N: Number of subjects following a strategy.
##'    \item value: Estimated strategy values.
##'    \item se: standard errors of estimation
##'    \item lower.CI: Lower bound of (1-alpha) level confidence interval for strategy values
##'    \item upper.CI: Upper bound of (1-alpha) level confidence interval for strategy values
##'    }
#' An object of ``vmat'' is return, which is variance-covariance matrix of estimated strategy values
#'
#' @references Lavori P.W. and Dawson R. (2007). Improving the efficiency of estimation in randomization trials of adaptive treatment strategies. \emph{Clinical Trials}, 4: 297-308.
#' @references Ko and Wahed A.S. (2015). Design of sequentially randomization trials for testing adaptive treatment strategies. \emph{Statistics in Medicine}, 31, 812-830.
#'
#' @examples
#' #generate a pesudo SMART data
#' N=8000                       # total sample size
#' A1=O2=A2=Y=rep(NA,N)
#' Dat=data.frame(A1,O2,A2,Y)
#'
#' # stage-1 treatment
#' Dat$A1=sample(c(0,1),size=N,prob=c(1,1),replace=TRUE)
#'
#' # intermediate outcome
#' n0=length(Dat$A1[which(Dat$A1==0)]); n1=N-n0
#' Dat$O2[which(Dat$A1==0)]=rbinom(n=n0,size=1,p=0.6)
#' Dat$O2[which(Dat$A1==1)]=rbinom(n=n1,size=1,p=0.4)
#'
#' # stage-2 treatment
#' n00=nrow(Dat[which(Dat$A1==0 & Dat$O2==0),]); n01=n0-n00
#' n10=nrow(Dat[which(Dat$A1==1 & Dat$O2==0),]); n11=n1-n10
#' Dat$A2[which(Dat$A1==0 & Dat$O2==0)]=sample(c(0,1),size=n00,prob=c(1,1),replace=TRUE)
#' Dat$A2[which(Dat$A1==0 & Dat$O2==1)]=sample(c(0,1),size=n01,prob=c(1,1),replace=TRUE)
#' Dat$A2[which(Dat$A1==1 & Dat$O2==0)]=sample(c(0,1),size=n10,prob=c(1,1),replace=TRUE)
#' Dat$A2[which(Dat$A1==1 & Dat$O2==1)]=sample(c(0,1),size=n11,prob=c(1,1),replace=TRUE)
#'
#' n000=nrow(Dat[which(Dat$A1==0 & Dat$O2==0 & Dat$A2==0),]); n001=n00-n000
#' n010=nrow(Dat[which(Dat$A1==0 & Dat$O2==1 & Dat$A2==0),]); n011=n01-n010
#' n100=nrow(Dat[which(Dat$A1==1 & Dat$O2==0 & Dat$A2==0),]); n101=n10-n100
#' n110=nrow(Dat[which(Dat$A1==1 & Dat$O2==1 & Dat$A2==0),]); n111=n11-n110
#'
#' # final primary outcome
#' Dat$Y[which(Dat$A1==0 & Dat$O2==0 & Dat$A2==0)]=rnorm(n=n000,mean=7,sd=6)
#' Dat$Y[which(Dat$A1==0 & Dat$O2==0 & Dat$A2==1)]=rnorm(n=n001,mean=7,sd=6)
#' Dat$Y[which(Dat$A1==0 & Dat$O2==1 & Dat$A2==0)]=rnorm(n=n010,mean=7,sd=8)
#' Dat$Y[which(Dat$A1==0 & Dat$O2==1 & Dat$A2==1)]=rnorm(n=n011,mean=7,sd=8)
#' Dat$Y[which(Dat$A1==1 & Dat$O2==0 & Dat$A2==0)]=rnorm(n=n100,mean=7,sd=6)
#' Dat$Y[which(Dat$A1==1 & Dat$O2==0 & Dat$A2==1)]=rnorm(n=n101,mean=7,sd=6)
#' Dat$Y[which(Dat$A1==1 & Dat$O2==1 & Dat$A2==0)]=rnorm(n=n110,mean=7,sd=8)
#' Dat$Y[which(Dat$A1==1 & Dat$O2==1 & Dat$A2==1)]=rnorm(n=n111,mean=7,sd=8)
#'
#' atsmeans(data=Dat,family="normal",method="Gest",conf=TRUE,common=TRUE,alpha=0.05,plot=TRUE)
#'
#' @export

atsmeans=function(data,family="normal",method="Gest",common=FALSE,conf=TRUE,alpha=0.05,plot=FALSE){
         D=as.data.frame(data); FA=family; Ma=method; Com=common; Al=alpha
         if (is.null(D$O1)) {Base=0} else {Base=1}
         Nstage=nstage(data=D)

         Umat=em(data=D,method=Ma)
         Val=Umat$value
         Vmat=round(evcmat(data=D,family=FA,method=Ma,common=Com),4)
         se=round(sqrt(diag(Vmat)),2)
         CIs=round(atsci(eumat=Umat,evmat=Vmat,alpha=Al),2)
         if (conf==FALSE) {Umat=data.frame(Umat,se)} else
         if (conf==TRUE) {Umat=data.frame(Umat,se,CIs)}
         message(paste("$value: estimated strategy values (with confidence intervals)",
                       "$vmat: variance-covariance matrix of estimated strategy values \n",
                       sep="\n"))
         if (Nstage==1 & Base==0) {message("A strategy is defined as a single-stage decision making (d0) for A1 at baseline")
                                   opar=par(mar=c(4,4,4,3))
                                   on.exit(par(opar))} else
         if (Nstage==1 & Base==1) {message(paste("A strategy is defined as a vector of single-stage decision makings (d0,d1),",
                                           "each of which corresponds to a possible outcome of baseline evulation (O1). \n",
                                           "d0 is the stage-1 decision making for A1, conditioning on O1=0",
                                           "d1 is the stage-1 decision making for A1, conditioning on O1=1",
                                           sep="\n"))
                                   opar=par(mar=c(5,4,4,3))
                                   on.exit(par(opar))
                                   } else
          if (Nstage==2 & Base==0) {message(paste("A strategy is defined as a vector of decision makings (d0;d00,d01) for 2 stages \n",
                                                  "d0 is the stage-1 decision making for A1",
                                                  "d00 is the stage-2 decision making for A2, conditioning on A1=d0 and O2=0",
                                                  "d01 is the stage-2 decision making for A2, conditioning on A1=d0 and O2=0",
                                                  sep="\n"))
                                    opar=par(mar=c(6,4,4,3))
                                    on.exit(par(opar))
                                    } else                                     
                                     
         if (Nstage==2 & Base==1) {message(paste("A strategy is defined as a vector of decision makings (d0,d1;d00,d01,d10,d11) \n",
                                                 "d0 is the stage-1 decision making conditioning on O1=0",
                                                 "d1 is the stage-1 decision making conditioning on O1=1",
                                                 "d00 is the stage-2 decision making conditioning on A1=d0 and O2=0",
                                                 "d01 is the stage-2 decision making conditioning on A1=d0 and O2=0",
                                                 "d00 is the stage-2 decision making conditioning on A1=d1 and O2=0",
                                                 "d01 is the stage-2 decision making conditioning on A1=d1 and O2=0",
                                                 sep="\n"))
                                   opar=par(mar=c(7,4,4,3))
                                   on.exit(par(opar))
                                   } else
         if (Nstage==3 & Base==0) {message(paste("A strategy is defined as a vector of decision makings (d0;d00,d01;d000,d001,d010,d111) \n",
                                                 "d0 is the stage-1 decision making",
                                                 "d00 is the stage-2 decision making conditioning on A1=d0 and O2=0",
                                                 "d01 is the stage-2 decision making conditioning on A1=d0 and O2=0",
                                                 "d000 is the stage-3 decision making conditioning on A1=d0, O2=0, A3=d00 and O3=0",
                                                 "d001 is the stage-3 decision making conditioning on A1=d0, O2=0, A3=d00 and O3=1",
                                                 "d010 is the stage-3 decision making conditioning on A1=d0, O2=1, A3=d01 and O3=0",
                                                 "d011 is the stage-3 decision making conditioning on A1=d0, O2=1, A3=d01 and O3=1",
                                                 sep="\n"))
                                   opar=par(mar=c(8,4,4,3))
                                   on.exit(par(opar))
                                   } else
         if (Nstage==3 & Base==1) {message(paste("A strategy is defined as a vector of decision makings (d0,d1;d00,d01,d10,d11;d000,d001,d010,d011,d100,d101,d110,d111) \n",
                                                 "d0 is the stage-1 decision making conditioning on O1=0",
                                                 "d1 is the stage-1 decision making conditioning on O1=1",
                                                 "d00 is the stage-2 decision making conditioning on A1=d0 and O2=0",
                                                 "d01 is the stage-2 decision making conditioning on A1=d0 and O2=0",
                                                 "d00 is the stage-2 decision making conditioning on A1=d1 and O2=0",
                                                 "d01 is the stage-2 decision making conditioning on A1=d1 and O2=0",
                                                 "d000 is the stage-3 decision making conditioning on A1=d0, O2=0, A3=d00 and O3=0",
                                                 "d001 is the stage-3 decision making conditioning on A1=d0, O2=0, A3=d00 and O3=1",
                                                 "d010 is the stage-3 decision making conditioning on A1=d0, O2=1, A3=d01 and O3=0",
                                                 "d011 is the stage-3 decision making conditioning on A1=d0, O2=1, A3=d01 and O3=1",
                                                 "d100 is the stage-3 decision making conditioning on A1=d1, O2=0, A3=d10 and O3=0",
                                                 "d101 is the stage-3 decision making conditioning on A1=d1, O2=0, A3=d10 and O3=1",
                                                 "d110 is the stage-3 decision making conditioning on A1=d1, O2=1, A3=d11 and O3=0",
                                                 "d111 is the stage-3 decision making conditioning on A1=d1, O2=1, A3=d11 and O3=1",
                                                 sep="\n"))
                                   opar=par(mar=c(8,4,4,3))
                                   on.exit(par(opar))
                                   }
         message(" \n")

         if (plot==TRUE) {atsciplot(uimat=Umat,nstage=Nstage,baseline=Base)}
         oopar=par(mar=c(4,4,4,4))
         on.exit(par(oopar))
         return(list(value=round(Umat,2),vmat=round(Vmat,4)))
}



