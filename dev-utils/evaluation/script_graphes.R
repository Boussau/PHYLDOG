cextant = read.table("summary_content_extant")
dextant = read.table("summary_degres_extant")



cstart = read.table("summary_content_deco_startingtrees")
crec = read.table("summary_content_deco_reconciledtrees")
dstart = read.table("summary_degres_deco_startingtrees")
drec = read.table("summary_degres_deco_reconciledtrees")

jpeg("summary_degre.jpg")
# plot(dextant$V1,dextant$V2/sum(cextant$V1),col="black",cex=1,xlim=c(0,7),xlab="Number of neighbors",ylab="Proportion of genes",lwd=1,type="l")
plot(dextant$V1,dextant$V2/sum(cextant$V2),col="black",xlab="Number of neighbors",ylab="Proportion of genes",lwd=1,type="l")
legend(5,0.3, c("reconciled","starting"),col = c("red","blue"),lwd=3,cex=1)
points(dstart$V1,dstart$V2/sum(cstart$V2),lwd=3,col="blue",type="l")
points(drec$V1,drec$V2/sum(crec$V2),lwd=3,col="red",type="l")
# points(dproens$V1,dproens$V2/sum(cproens$V2),lwd=3,col="red",type="l")
# points(dprofile95$V1,dprofile95$V2/sum(cprofile95$V1),lwd=3,col="red",type="l")
# points(dphyldiag$V1,dphyldiag$V2/sum(cphyldiag$V1),lwd=3,col="green",type="l")
dev.off()

# jpeg("summary_degre_random.jpg")
# plot(dextant$V1,dextant$V2/sum(cextant$V1),col="black",cex=1,xlim=c(0,7),xlab="Number of neighbors",ylab="Proportion of genes",lwd=1,type="l")
# legend(5,0.3, c("ProfileNJ","ProfileRandom","Extant"),col = c("red","blue","black"),lwd=3,cex=1)
# points(dprofile95$V1,dprofile95$V2/sum(cprofile95$V1),lwd=3,col="red",type="l")
# points(dprrandom$V1,dprrandom$V2/sum(cprrandom$V1),lwd=3,col="green",type="l")
# dev.off()


jpeg("summary_content.jpg")
boxplot(cextant$V2,cstart$V2,crec$V2,names = c("extant","starting","reconciled"),ylab="Number of genes")
dev.off()

