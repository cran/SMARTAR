#' @importFrom stats pnorm

pwtest=function(data,family="normal",method="Gest",alpha=0.05,common=FALSE,adjust=NULL,ntest=NULL){
       D=as.data.frame(data); FA=family; Ma=method; Com=common; Nt=ntest
       if (is.null(D$O1)) {Base=0} else {Base=1}
       Dmat=atsscan(data=D)
       G=nrow(Dmat)
       label=getag(nats=G)
       Val=em(data=D,method=Ma)$value
       Umat=matrix(Val,ncol=1)
       Vmat=evcmat(data=D,family=FA,method=Ma,common=Com)
       PEmat=LOmat=UPmat=Zmat=Pmat=matrix(rep(NA,(G-1)*G),ncol=G)
       if (is.null(adjust)) {Q=round(qnorm(1-alpha/2),2)} else
       if (adjust=="Bon" & is.null(Nt)) {Q=round(qnorm(1-alpha/((G-1)*G)),2)} else
       if (adjust=="Bon" & Nt>0) {Q=round(qnorm(1-alpha/Nt),2)}
       for (j in 1:G){
            Cmat=getcmat(control=j,nats=G)
            PEmat[,j]=Cmat%*%Umat
            DVD=Cmat%*%Vmat%*%t(Cmat)
            VAR.ig=diag(DVD)
            LOmat[,j]=PEmat[,j]-Q*sqrt(VAR.ig)
            UPmat[,j]=PEmat[,j]+Q*sqrt(VAR.ig)
            Zmat[,j]=PEmat[,j]/sqrt(VAR.ig)
            Pmat[,j]=(1-pnorm(abs(Zmat[,j])))*2
       }
       Lab=paste(label[,1],"vs.",label[,2])
       pwout=as.data.frame(cbind(round(as.vector(PEmat),2),
                                 round(as.vector(LOmat),2),
                                 round(as.vector(UPmat),2),
                                 round(as.vector(Zmat),2),
                                 round(as.vector(Pmat),4)))
       pwout=cbind.data.frame(Lab,pwout)
       colnames(pwout)=c("label","diff","lower.CI","upper.CI","Z","Pvalue")
       return(pwout)
       #return(Lab)
}
