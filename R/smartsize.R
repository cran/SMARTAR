#' sample size calculation
#'
#' Return a message that contains the estimated strategy-specified means and their confidence interval, as well as the asymptotic variance-covariance matrix for these estimates.
#' @param sim A numeric matrix containing the values of treatment sequence-specific parameters to generate the SMART data, including the values of stage-specific treatments, intermediate outcome and final primary outcome.
#' @param family A character string to specify the type of final primary outcome. The default is family=“normal”, which refers to the continuous primary outcome. If family=”binomial” then the primary outcome will be treated as binary variable.
#' @param delta The standardized effect size for sample size calculation.
#' @param df The degrees of freedom for the chisquare test.
#' @param alpha Type I error rate.
#' @param beta Type II error rate.
#' @param global If TRUE then power the SMART based on a global test, otherwise power the SMART based on a pairwise test. The default is TRUE
#' @param method A characteristic string to specify the method of estimation. We provide 2 method of estimation. The default is method=”MLE”, which is the Maximized Likelihood Estimation method. If method=”IPWN” then use the Normalized Inversed Probability Weighted method. The input of this option is only needed when power the SMART based on global test but not the pairwise test.
#' @return Standardized effect size and total sample size for SMART
##' \itemize{
##'    \item delta: standardized effect size
##'    \item n: total sample size
##'    }
#' @references Murphy, S. A. (2005). An experimental design for the development of adaptive treatment strategies. Statistics in Medicine. 24(10): 1455-1481.
#' @references Ogbagaber S. B., Karp, J., and Wahed A.S. (2015). Design of sequentially randomization trials for testing adaptive treatment strategies. Statistics in Medicine. DOI: 10.1002/sim.6747.
#' @importFrom MASS ginv
#' @export

smartsize=function(sim=NULL,delta=NULL,df=NULL,alpha=0.05,beta=0.20,global=TRUE,family="normal",method="Gest"){
          Smat=sim; Al=alpha; Be=beta; MA=method; FA=family; DF=df; DELTA=delta
          if (global){
              if (is.null(DELTA)) {Dmat=atsscan(data=Smat)
                                   G=nrow(Dmat)
                                   Vmat=tvmat(sim=Smat,family=FA,method=MA)
                                   Umat=tumat(sim=Smat,family=FA)
                                   Cmat=getcmat(control=1,nats=G)
                                   MAT=ginv(Cmat%*%Vmat%*%t(Cmat))
                                   Del=t(Umat)%*%t(Cmat)%*%MAT%*%Cmat%*%Umat
                                   DF=getdf(data=Smat)
                                   Lam=getncp(alpha=Al,beta=Be,df=DF)
                                   if (Del<=0) {warning("Effect size cannot be less than or equal to 0")}
                                          else {size=Lam/Del}
                                   }
                             else {Lam=getncp(alpha=Al,beta=Be,df=DF)
                                   Del=DELTA
                                   if (Del<=0) {warning("Effect size cannot be less than or equal to 0")}
                                          else {size=Lam/Del}
                                   }
              N=ceiling(size)
              out=as.data.frame(cbind(Lam,round(Del,3),DF,N))
              colnames(out)=c("NCP","delta","df","N")
              message("The sample size is for total subjects registered in the trial. \n")
              }
          else if (global==FALSE) {Del=DELTA
                                   if (Del<=0) {warning("Effect size cannot be less than or equal to 0")}
                                          else {size=(qnorm(1-Al/2)+qnorm(1-Be))^2*2*(1/Del)^2}
                                   N=ceiling(size)
                                   out=as.data.frame(cbind(Del,N))
                                   colnames(out)=c("delta","N")
                                   message("The sample size is for subjects used in pairwise test. \n")
                                   }
          return(out)
}
