cextant = read.table("summary_content_extant")
dextant = read.table("summary_degres_extant")



censembl = read.table("summary_content_deco_ensemblTrees")
cstart = read.table("summary_content_deco_startingTrees")
dstart = read.table("summary_degres_deco_startingTrees")
crec = read.table("summary_content_deco_reconciledTrees")
densembl = read.table("summary_degres_deco_ensemblTrees")
drec = read.table("summary_degres_deco_reconciledTrees")

pdf("summary_degre.pdf")
plot(dextant$V1,dextant$V2/sum(cextant$V2),col="black",xlab="Number of neighbors",ylab="Proportion of genes",lwd=1,type="o",cex=1,xlim=c(0,4))
legend("topright", c("reconciled","ensembl","extant","starting"),col = c("red","blue","black","green"),lwd=3,cex=1)
points(densembl$V1,densembl$V2/sum(censembl$V2),lwd=3,col="blue",type="o")
points(dstart$V1,dstart$V2/sum(cstart$V2),lwd=3,col="green",type="o")

points(drec$V1,drec$V2/sum(crec$V2),lwd=3,col="red",type="o")

dev.off()

pdf("summary_content.pdf")
boxplot(cextant$V2,censembl$V2,crec$V2,cstart$V2,names = c("extant","ensembl","reconciled","starting"),ylab="Number of genes")
dev.off()

lk_before_after = read.csv("../lk_before_and_after.csv",h=F)
lk_ensembl = read.csv("../lk_ensembl.csv",h=F)
lk_total = merge(lk_before_after,lk_ensembl, by="V1")
colnames(lk_total) <- c("name","before","after","ensembl","ensemblAfter")

pdf("summary_lk.pdf")
#plot(lk_before_after$V2,lk_before_after$V3,xlab="Before",ylab="After")
#plot(lk_before_after$V3,,xlab="Before",ylab="After")
hist((lk_total$after-lk_total$before)/lk_total$before,main="Normalized likelihood gain",xlab="(LKafter-LKbefore)/LKbefore",breaks=30)

#comparaison avec ensembl
#negatif = bien
hist((lk_total$after-lk_total$ensemblAfter)/lk_total$ensemblAfter,main="Normalized likelihood gain",xlab="(LKafter-LKensemblAfter)/LKensemblAfter",breaks=30)

dev.off()

