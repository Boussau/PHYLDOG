options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)
print(args)

refDir = args[1]
refName = tail(unlist(strsplit(refDir,"/")),n=1)
refName = substr(refName,5,nchar(refName))

compDir = args[2]
compName = tail(unlist(strsplit(compDir,"/")),n=1)
compName = substr(compName,5,nchar(compName))

destDir = args[3]


cextant = read.table(paste(refDir, "/summary_content_extant",sep=""))
dextant = read.table(paste(refDir, "/summary_degres_extant",sep=""))

cRefStart = read.table(paste(refDir, "/summary_content_deco_starting",sep=""))
dRefStart = read.table(paste(refDir, "/summary_degres_deco_starting",sep=""))

cRefRec = read.table(paste(refDir, "/summary_content_deco_reconciled",sep=""))
dRefRec = read.table(paste(refDir, "/summary_degres_deco_reconciled",sep=""))


cCompRec = read.table(paste(compDir, "/summary_content_deco_reconciled",sep=""))
dCompRec = read.table(paste(compDir, "/summary_degres_deco_reconciled",sep=""))



pdf(paste(destDir, "/",refName,"_challengedBy_",compName,".pdf",sep=""))


## DEGREE

plot(dextant$V1,dextant$V2/sum(cextant$V2),col="black",xlab="Number of neighbors",ylab="Proportion of genes",lwd=1,type="o",cex=1,xlim=c(0,4),main=paste("Degree:",refName))
points(dRefStart$V1,dRefStart$V2/sum(cRefStart$V2),lwd=3,col="green",type="o")
points(dCompRec$V1,dCompRec$V2/sum(cCompRec$V2),lwd=3,col="blue",type="o")
points(dRefRec$V1,dRefRec$V2/sum(cRefRec$V2),lwd=3,col="red",type="o")
legend("topright", c(refName,compName,"extant","starting"),col = c("red","blue","black","green"),lwd=3,cex=1)


boxplot(cextant$V2,cRefRec$V2,cCompRec$V2,cRefStart$V2,names = c("extant",refName,compName,"starting"),ylab="Number of genes",main=paste("Genome size:",refName))

boxplot(cextant$V2,cRefRec$V2,cCompRec$V2,names = c("extant",refName,compName),ylab="Number of genes",notch=T,main=paste("Genome size:",refName))


## LK GAIN Distribution: hist & density

lk_ref = read.csv(paste(refDir, "/collected_lk.csv",sep=""), h=F, col.names=c("family","refStarting","refReconciled","refRecSeq","refRecScen"))
lk_comp = read.csv(paste(compDir, "/collected_lk.csv",sep=""), h=F, col.names=c("family","compStarting","compReconciled","compRecSeq","compRecScen"))
lk_total = merge(lk_ref,lk_comp, by="family")

#REF starting vs reconciled
hist(((lk_total$refReconciled - lk_total$refStarting)/lk_total$refStarting),nclass=20,proba=T,main=paste("Lk Gain:",refName,"from starting"),xlab="(reconciled-starting)/starting")
lines(density((lk_total$refReconciled - lk_total$refStarting)/lk_total$refStarting),col="red")
mtext("Distribution more on the left: better improvement from starting to reconciled")

#REF vs COMP
hist(((lk_total$refReconciled - lk_total$compReconciled)/lk_total$compReconciled),nclass=20,proba=T,main=paste("Lk Gain:",refName,"from",compName),xlab=paste("(",refName,"-",compName,")/",compName))
lines(density((lk_total$refReconciled - lk_total$compReconciled)/lk_total$compReconciled,adjust=4),col="red")
mtext(paste("Left of 0 means",refName,"is better than",compName))


## LK GAIN plot

# TOTAL LLk

barplot(lk_total$refReconciled-lk_total$compReconciled,main=paste("TOTAL LLk:",refName))
mtext(paste("Under Zero means",refName,"iis better than",compName))

# SEQUENCE LLk

barplot(lk_total$refRecSeq-lk_total$compRecSeq,main=paste("SEQUENCE LLk:",refName))
mtext(paste("Under Zero means",refName,"iis better than",compName))

barplot(lk_total$refRecScen-lk_total$compRecScen,main=paste("SCENARIO LLk",refName))
mtext(paste("Under Zero means",refName,"iis better than",compName))

## return of the HIST & DENSITY

#REF vs COMP: sequence
hist(((lk_total$refRecSeq - lk_total$compRecSeq)/lk_total$compRecSeq),nclass=20,proba=T,main="SEQUENCE Lk Gain: REF from COMP",xlab=paste("(",refName,"-",compName,")/",compName))
lines(density((lk_total$refRecSeq - lk_total$compRecSeq)/lk_total$compRecSeq,adjust=4),col="red")
mtext(paste("Left of 0 means",refName,"is better than",compName))

#REF vs COMP: scenario
hist(((lk_total$refRecScen - lk_total$compRecScen)/lk_total$compRecScen),nclass=20,proba=T,main="SCENARIO Lk Gain: REF from COMP",xlab=paste("(",refName,"-",compName,")/",compName))
lines(density((lk_total$refRecScen - lk_total$compRecScen)/lk_total$compRecScen,adjust=4),col="red")
mtext(paste("Left of 0 means",refName,"is better than",compName))

# # lk_ensembl = read.csv("../lk_ensembl.csv",h=F)
# # lk_total = merge(lk_before_after,lk_ensembl, by="V1")
# colnames(lk_total) <- c("name","before","after")
# 
# pdf("summary_lk.pdf")
# plot(lk_total$before,lk_total$after,xlab="Before",ylab="After")
# #plot(lk_before_after$V3,,xlab="Before",ylab="After")
# hist((lk_total$after-lk_total$before)/lk_total$before,main="Normalized likelihood gain",xlab="(LKafter-LKbefore)/LKbefore",breaks=30)
# 
# #comparaison avec ensembl
# #negatif = bien
# # hist((lk_total$after-lk_total$ensemblAfter)/lk_total$ensemblAfter,main="Normalized likelihood gain",xlab="(LKafter-LKensemblAfter)/LKensemblAfter",breaks=30)
# 



dev.off()
