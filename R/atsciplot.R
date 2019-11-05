#' @importFrom graphics par plot title axis lines

#atsciplot() plot the PEs and (1-a)% CIs of all strategy values
atsciplot=function(uimat,col="forestgreen",nstage,baseline=0,lim=NULL,cex=2,lwd=3){
          UImat=uimat; Nstage=nstage; Base=baseline; COL=col; L=lim; CEX=cex; LWD=lwd
          G=nrow(uimat)

          if (Nstage==1 & Base==0) {lab=paste("(",uimat$d0,")")}
          if (Nstage==1 & Base==1) {lab=paste("(",uimat$d0,uimat$d1,")")}
          if (Nstage==2 & Base==0) {lab=paste("(",uimat$d0,
                                                  uimat$d00,uimat$d01,")")}
          if (Nstage==2 & Base==1) {lab=paste("(",uimat$d0, uimat$d1,
                                                  uimat$d00,uimat$d01,uimat$d10,uimat$d11,")")}
          if (Nstage==3 & Base==0) {lab=paste("(",uimat$d0,
                                                  uimat$d00,uimat$d01,
                                                  uimat$d000,uimat$d001,uimat$d010,uimat$d011,")")}
          if (Nstage==3 & Base==1) {lab=paste("(",uimat$d0, uimat$d1,
                                                  uimat$d00,uimat$d01,uimat$d10,uimat$d11,
                                                  uimat$d000,uimat$d001,uimat$d010,uimat$d011,
                                                  uimat$d100,uimat$d101,uimat$d110,uimat$d111,")")}

          range=max(uimat$upper)-min(uimat$lower)
          mi=min(uimat$lower)-range/2
          ma=max(uimat$upper)+range/2

          if (is.null(L)) {L=c(mi,ma)} else {L=lim}
          plot(x=1:G,y=uimat$value,xlim=c(0,G+1),pch=15,col=COL,cex=CEX,xaxt="n",xlab="",ylab="",ylim=L)
          title(main="Strategy values with confidence interval (C.I.)",ylab="Strategy value (C.I.)")
          axis(side=1,1:G,labels=lab,las=2)
          CI=as.matrix(cbind(UImat$lower,UImat$upper))
          for (i in 1:G){lines(x=c(i,i),y=c(CI[i,1],CI[i,2]),lwd=LWD,col=COL)}
}
