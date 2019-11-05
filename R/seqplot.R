#' @importFrom graphics barplot boxplot title par
#seqplot() output the plot for sequence-specific descriptive statistics

seqplot=function(data,family="normal"){
        D=as.data.frame(data)
        if (is.null(D$O1)) {Base=0} else {Base=1}

        N=nrow(D)
        Nstage=nstage(data=D)
        Smat=seqscan(data=D)
        ns=nrow(Smat)

        SEQ=rep(NA,N)      #assign sequence index for each subject in input dataset
        D=cbind(D,SEQ)

        if (Nstage==1 & Base==0) {for (j in 1:ns) {D$SEQ[which(D$A1==Smat[j,2])]=j}
                                  lab=paste("(",Smat$A1,")")
                                  ti="Primary outcome by sequence (A1)"} else
        if (Nstage==1 & Base==1) {for (j in 1:ns) {D$SEQ[which(D$O1==Smat[j,2] & D$A1==Smat[j,3])]=j}
                                  lab=paste("(",Smat$O1,Smat$A1,")")
                                  ti="Primary outcome by sequence (O1,A1)"} else
        if (Nstage==2 & Base==0) {for (j in 1:ns) {D$SEQ[which(D$A1==Smat[j,2] & D$O2==Smat[j,3] & D$A2==Smat[j,4])]=j}
                                  lab=paste("(",Smat$A1,Smat$O2,Smat$A2,")")
                                  ti="Primary outcome by treatment sequence (A1,O2,A2)"} else
        if (Nstage==2 & Base==1) {for (j in 1:ns) {D$SEQ[which(D$O1==Smat[j,2] & D$A1==Smat[j,3] & D$O2==Smat[j,4] & D$A2==Smat[j,5])]=j}
                                  lab=paste("(",Smat$O1,Smat$A1,Smat$O2,Smat$A2,")")
                                  ti="Primary outcome by treatment sequence (O1,A1,O2,A2)"} else
        if (Nstage==3 & Base==0) {for (j in 1:ns) {D$SEQ[which(D$A1==Smat[j,2] & D$O2==Smat[j,3] & D$A2==Smat[j,4] & D$O3==Smat[j,5] & D$A3==Smat[j,6])]=j}
                                  lab=paste("(",Smat$A1,Smat$O2,Smat$A2,Smat$O3,Smat$A3,")")
                                  ti="Primary outcome by treatment sequence (A1,O2,A2,O2,A3)"} else
        if (Nstage==3 & Base==1) {for (j in 1:ns) {D$SEQ[which(D$O1==Smat[j,2] & D$A1==Smat[j,3] & D$O2==Smat[j,4] & D$A2==Smat[j,5] & D$O3==Smat[j,6] & D$A3==Smat[j,7])]=j}
                                  lab=paste("(",Smat$O1,Smat$A1,Smat$O2,Smat$A2,Smat$O3,Smat$A3,")")
                                  ti="Primary outcome by treatment sequence (O1,A1,O2,A2,O3,A3)"}

        oldbar=par(mar=c(4,4,4,2))
        on.exit(par(oldbar))                    
        if (family=="binomial") {
            counts=table(D$Y,D$SEQ)
            range=c(0,max(table(D$SEQ))*3/2)
            barplot(counts,col=c("yellow","forestgreen"),legend=rownames(counts),
                    names.arg=lab,las=2,ylim=range,xlab=NA,ylab="Primary outcome")} else
        if (family=="normal") {
            ran=max(D$Y)-min(D$Y)
            range=c(min(D$Y)-ran/4,max(D$Y)+ran/4)
            boxplot(D$Y~D$SEQ,horizontal=F,col="yellow",names=lab,las=2,pch=19,cex=0.8,ylim=range,
                    xlab=NA,ylab="Primary outcome")}
            title(main=ti,cex=0.8)
}

#seqplot(data=CD1)
#seqplot(data=OCD1)
#seqplot(data=BCD1,family="binomial")
#seqplot(data=BOCD1,family="binomial")

#seqplot(data=CD2)
#seqplot(data=OCD2)
#seqplot(data=BCD2,family="binomial")
#seqplot(data=BOCD2,family="binomial")

#seqplot(data=CD3)
#seqplot(data=OCD3)
#seqplot(data=BCD3,family="binomial")
#seqplot(data=BOCD3,family="binomial")

#seqplot(data=CODIACS)
#seqplot(data=smoke,family="binomial")
#seqplot(data=QOL123)
