#' Conduct statistical tests using a SMART data
#'
#' Return a message that contains the results of statistical tests to compare the values of adaptive treatment strategies defined in a SMART data. The statistical tests include (1) a global test (2) a series of pairwise tests.
#' @param data Input data frame of the SMART data used for analysis, which include the variables of stage-1 treatments (A1), intermediate outcome (O2), stage-2 treatment (A2) and final primary outcome (Y). If stage-1 treatment takes into account baseline information, baseline information (O1) also needs to be included.
#' @param family A character string to specify the type of final primary outcome. The default is family=“normal”, which refers to the continuous primary outcome. If family=”binomial” then the primary outcome will be treated as binary variable.
#' @param method Method used to estimate the value of adaptive treatment strategy. "Gest" for G-computation method and "IPW" for Inversed Probabiliy Weight method. Default is method="Gest".
#' @param common If common=TRUE, the pooled variance across all the treatment sequences are used in estimation. Otherwise use the sequence-specific variance. The default is common=FALSE.
#' @param alpha Significant level of confidence interval. The default is alpha=0.05.
#' @param adjust A characteristic string to indicate whether the confidence intervals pairwise distance adjusted for multiple comparison. The default is adjust=”n”, which indicated no adjustment for multiple comparison. If adjust=”Bon”, the CIs are adjusted for the Bonferroni correction.
#' @param ntest Number of pairwise tests adjusted for Bonferroni correction
#' @return
#' An objects of “Strategy” is return, which lists all the adaptive treatment strategy defined in the input data with an index number.
##' \itemize{
##'    \item ATS: Index of the treatment adaptive treatment strategy defined in the input dataset
##'    \item ds: the sequence of decision makings that define an adaptive treatment strategy
##'    \item N: number of subjects following an adaptive treatment strategy in the intut dataset
##'    }
#' An objects of "Global.test" is return, which give the result of the global test.
##' \itemize{
##'    \item size: the total number of subjects in the input dataset
##'    \item nATS: the total number of adaptive treatment strategies defined in the input dataset
##'    \item df: the degrees of freedom of the global test, which is a chisquare test
##'    \item chisq: the chisquare test statistics for the global test
##'    \item Pvalue: the P-value of the global test
##'    }
#' An object of "Pairwise.test"
##' \itemize{
##'    \item label: The labels of pairwise tests. The details of strategies are shown in $Strategy
##'    \item diff: Estimated pairwise distance between treatment and control adaptive treatment strategy
##'    \item lower.CI: Lower bound of confidence interval for pairwise distance
##'    \item upper.CI: Upper bound of confidence interval for pairwise distance
##'    \item Z: Test statistics of pairwise test
##'    \item P-value: P-value of pairwise test
##'    }
#' @references Murphy, S. A. (2005),``An experimental design for the development of adaptive treatment strategies'' \emph{Statistics in Medicine}, 24, 1455-1481.
#' @references Ogbagaber S. B., Karp, J., and Wahed A.S. (2015), ``Design of sequentially randomization trials for testing adaptive treatment strategies,'' \emph{Statistics in Medicine}, DOI: 10.1002/sim.6747.
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
#'smartest(data=Dat,family="normal",method="Gest",common=FALSE,alpha=0.05,adjust="Bon")
#'
#' @export

smartest=function(data,family="normal",method="Gest",common=FALSE,alpha=0.05,adjust=NULL,ntest=NULL){
         D=data; FA=family; Ma=method; Al=alpha; Aj=adjust; Ntest=ntest; Com=common
         Dmat=atsscan(data=D)
         chi=x2test(data=D,family=FA,method=Ma,common=Com)
         pw=pwtest(data=D,family=FA,method=Ma,alpha=Al,adjust=Aj,ntest=Ntest,common=Com)
         message(paste("$Strategy provides the details of decision makings under strategy labels (ATS)",
                       "$Global.test assesses the null hypothesis of no difference across all the strategy values",
                       "$Pairwise.test compares all the pairs of strategies, of which the labels are shown in $Strategy",
                       sep="\n"))
         if (is.null(Aj)) {message("\n")} else
         if (Aj=="Bon") {message("The P values should compare to the critical value adjusted for the Bonferroni correction \n")}
         return(list(Strategy=Dmat,Global.test=chi,Pairwise.comparisons=pw))
}

