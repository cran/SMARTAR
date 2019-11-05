#' Summarize sequence-specific descriptive statistics
#'
#' Return a message that lists all the treatment sequence embedded in SMART design and summarizes all the sequence-specific descriptive statistics. It also provide design diagram of SMART and graphs of sequence-specific descriptive statistics (boxplot for continuous primary outcome and barchart for binary primary outcome).
#' @importFrom stats var
#' @param data Input data frame of the sequential randomized trial (SMART) data used for analysis. The data should include the variables of stage-specific treatments (At; t=1,2,3), intermediate evaluation (Ot; t=1,2,3) and final primary outcome (Y), where t represent the number of stages embedded in design. If stage-1 treatment (A1) takes into account the information of baseline evaluation, O1 needed to be include in data, otherwise not.
#' @param family A character string to specify the type of final primary outcome. The default is family=“normal”, which refers to the continuous primary outcome. If family=”binomial” then the primary outcome will be treated as binary variable.
#' @param plot A character string to specify the output figure. If plot=”d” then output the design diagram of SMART; If plot=”s” then output boxplot for continuous primary outcome or bar plot for binary primary outcome by sequence. The default is plot=”d”.
#' @return
#' an object of information of all the treatment sequences and sequences-specific descriptive statistics defined in a SMART data
##' \itemize{
##'    \item SEQ: Index of treatment sequences.
##'    \item O1: Baseline evaluation outcome.
##'    \item A1: Action of stage-1 treatment.
##'    \item O2: Intermeidate outcome evaluated at the end of stage 1.
##'    \item A2: Action of stage-2 treatment.
##'    \item O3: Intermeidate outcome evaluated at the end of stage 2.
##'    \item A3: Action of stage-3 treatment.
##'    \item N: Number of subjects following a certain treatment sequence.
##'    \item MEAN: Sequence-specified sample mean.
##'    \item VAR: Sequence-specified sample variance.
##'    }
#' @references Thall P., Millikan R. and Sung H.G. (2000), ``Evaluating multiple treatment courses in clinical trials,'' \emph{Statistics in Medicine}, 19, 1011-1028
#' @references Lavori P. W. and Dawson R. (2000), ``A design for testing clinical strategies: biased adaptive within-subject randomization,'' \emph{Journal of the Royal Statistical Society, Series A}, 163, 29-38.
#' @references Murphy, S. A. (2005), ``An experimental design for the development of adaptive treatment strategies,” \emph{Statistics in Medicine}, 24, 1455-1481.
#'@examples
#'
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
#' #get descriptive statistics
#' seqmeans(data=Dat,family="normal",plot="d")
#'
#' @export

seqmeans=function(data,family="normal",plot="d"){
         D=data.frame(data); FA=family
         if (is.null(D$O1)) {Base=0} else {Base=1}
         Smat=seqscan(data=D)

         N=nrow(D);
         Nstage=nstage(data=D)
         ns=nrow(Smat)

         M=V=rep(NA,ns)
         Smat=data.frame(Smat,M,V)

         SEQ=rep(NA,N)      #assign sequence index for each subject in input dataset
         D=cbind(D,SEQ)

         message("Each subject followed one of the below treatment sequences during the trial.")
         if (Nstage==1 & Base==0) {for (j in 1:ns) {D$SEQ[which(D$A1==Smat[j,2])]=j}
                                   message("A treatment sequence is defined as a scalar of (A1). \n")} else
         if (Nstage==1 & Base==1) {for (j in 1:ns) {D$SEQ[which(D$O1==Smat[j,2] & D$A1==Smat[j,3])]=j}
                                   message("A treatment sequence is defined as a vector of values (O1,A1). \n")}else
         if (Nstage==2 & Base==0) {for (j in 1:ns) {D$SEQ[which(D$A1==Smat[j,2] & D$O2==Smat[j,3] & D$A2==Smat[j,4])]=j}
                                   message("A treatment sequence is defined as a vector of values (A1,O2,A2). \n")} else
         if (Nstage==2 & Base==1) {for (j in 1:ns) {D$SEQ[which(D$O1==Smat[j,2] & D$A1==Smat[j,3] & D$O2==Smat[j,4] & D$A2==Smat[j,5])]=j}
                                   message("A treatment sequence is defined as a vector of values (O1,A1,O2,A2). \n")} else
         if (Nstage==3 & Base==0) {for (j in 1:ns) {D$SEQ[which(D$A1==Smat[j,2] & D$O2==Smat[j,3] & D$A2==Smat[j,4] & D$O3==Smat[j,5] & D$A3==Smat[j,6])]=j}
                                   message("A treatment sequence is defined as a vector of values (A1,O2,A2,O3,A3). \n")} else
         if (Nstage==3 & Base==1) {for (j in 1:ns) {D$SEQ[which(D$O1==Smat[j,2] & D$A1==Smat[j,3] & D$O2==Smat[j,4] & D$A2==Smat[j,5] & D$O3==Smat[j,6] & D$A3==Smat[j,7])]=j}
                                   message("A treatment sequence is defined as a vector of values (O1,A1,O2,A2,O3,A3). \n")}

         for (s in 1:ns){sm=mean(D$Y[which(D$SEQ==s)])
                         sv=var(D$Y[which(D$SEQ==s)])

                         if (family=="normal") {Smat$M[which(Smat$SEQ==s)]=round(sm,2)
                                                Smat$V[which(Smat$SEQ==s)]=round(sv,2)} else
                         if (family=="binomial") {Smat$M[which(Smat$SEQ==s)]=round(sm,2)
                                                  Smat$V[which(Smat$SEQ==s)]=round(sm*(1-sm),2)}
                         }
         colnames(Smat)[c(ncol(Smat)-1,ncol(Smat))]=c("MEAN","VAR")
         if (plot=="s") {seqplot(data=D,family=FA)} else {ddplot(data=D)}
         return(Smat)
}
