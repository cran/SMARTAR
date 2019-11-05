#em() combine two method of estimation for DTR value
#"Gest" for method option indicates "G-estimation method"; "IPW" indicates "Inversed Probability Weighting method"

em=function(data,method="Gest"){
   D=as.data.frame(data)
   if (is.null(D$O1)) {Base=0} else {Base=1}
   Nstage=nstage(data=D)
   Dmat=atsscan(data=D)
   G=nrow(Dmat)                                               #total number of strategies defined in input dataset
   umat=NULL

   if (Nstage==1 & Base==0) {for (g in 1:G){
                                  if (method=="Gest") {eu=gemean(data=D,ats=Dmat[g,2])} else
                                  if (method=="IPW") {eu=memean(data=D,ats=Dmat[g,2])}
                                  umat=rbind(umat,c(Dmat[g,],eu))
                                  colnames(umat)=c("ATS","d0","N","value")}} else
   if (Nstage==1 & Base==1) {for (g in 1:G){
                                  if (method=="Gest") {eu=gemean(data=D,ats=Dmat[g,2:3])} else
                                  if (method=="IPW")  {eu=memean(data=D,ats=Dmat[g,2:3])}
                                  umat=rbind(umat,c(Dmat[g,],eu))
                                  colnames(umat)=c("ATS","d0","d1","N","value")}} else
   if (Nstage==2 & Base==0) {for (g in 1:G){
                                  if (method=="Gest") {eu=gemean(data=D,ats=Dmat[g,2:4])} else
                                  if (method=="IPW")  {eu=memean(data=D,ats=Dmat[g,2:4])}
                                  umat=rbind(umat,c(Dmat[g,],eu))
                                  colnames(umat)=c("ATS","d0","d00","d01","N","value")}} else
   if (Nstage==2 & Base==1) {for (g in 1:G){
                                  if (method=="Gest") {eu=gemean(data=D,ats=Dmat[g,2:7])} else
                                  if (method=="IPW")  {eu=memean(data=D,ats=Dmat[g,2:7])}
                                  umat=rbind(umat,c(Dmat[g,],eu))
                                  colnames(umat)=c("ATS","d0","d1","d00","d01","d10","d11","N","value")}} else
   if (Nstage==3 & Base==0) {for (g in 1:G){
                                  if (method=="Gest") {eu=gemean(data=D,ats=Dmat[g,2:8])} else
                                  if (method=="IPW")  {eu=memean(data=D,ats=Dmat[g,2:8])}
                                  umat=rbind(umat,c(Dmat[g,],eu))
                                  colnames(umat)=c("ATS","d0","d00","d01","d000","d001","d010","d011","N","value")}} else
   if (Nstage==3 & Base==1) {for (g in 1:G){
                                  if (method=="Gest") {eu=gemean(data=D,ats=Dmat[g,2:15])} else
                                  if (method=="IPW")  {eu=memean(data=D,ats=Dmat[g,2:15])}
                                  umat=rbind(umat,c(Dmat[g,],eu))
                                  colnames(umat)=c("ATS","d0","d1","d00","d01","d10","d11","d000","d001",
                                                   "d010","d011","d100","d101","d110","d111","N","value")}}
   umat=as.data.frame(umat)
   return(umat)
}

#em(data=CD1)
#em(data=OCD1)

#em(data=CD2)
#em(data=OCD2)

#em(data=CD3)
#em(data=OCD3)

#em(data=CODIACS)
#em(data=smoke)
#em(data=QOL123)


