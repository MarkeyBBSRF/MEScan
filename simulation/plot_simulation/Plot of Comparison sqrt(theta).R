pdf("figure1a_3genes.pdf")
mypar <- function(a=1,b=1){
  par(mar=c(2.5,2.5,1.6,1.1),mgp=c(1.5,.5,0))
  par(mfrow=c(a,b))
}
mypar(2,2)
####################### balanced
frac = c(0.1,0.2,0.3,0.4)
rTg = c(80,100,100,100)/100
rTimx = c(8,90,100,100)/100
re_test = c(22,  87 , 96,100)/100   ###unweighted method
rMESGA =c(8,62,70,87)/100  
rDen = c(8,85,94,100)/100
rCoMEt = c(21,89,99,100)/100   #weighted method
########
plot(frac,rTg,type="b",cex=1.5,pch=15, col="red",
     ylim=c(0,1.05),xlab="coverage",ylab="power",
     cex.lab=1.5, cex.axis=1, cex.main=1.5, cex.sub=1.5)
points(frac,rTimx,type="b",cex=1.5,pch=16,col="blue")
points(frac,rCoMEt,type="b",cex=1.5,pch=17,col="green")
points(frac,rDen,type="b",cex=2,pch=18,col="purple")
points(frac,rMESGA,type="b",cex=1.5,pch=7,col="brown")
points(frac,re_test,type="b",cex=1.5,pch=10,col="orange")

############################imbalanced 
rTg = c(64,97,100,100)/100
rTimx = c(7,84,97,100)/100  #55
rCoMEt = c(9,84,97,100)/100
rMESGA =c(6,73,94,100)/100
rDen = c(7,59,88,100)/100
re_test = c( 9,  83 , 94 ,97)/100
########
plot(frac,rTg,type="b",cex=1.5,pch=15, col="red",
     ylim=c(0,1.05),xlab="coverage",ylab="power",
     cex.lab=1.5, cex.axis=1, cex.main=1.5, cex.sub=1.5)
points(frac,rTimx,type="b",cex=1.5,pch=16,col="blue")
points(frac,rCoMEt,type="b",cex=1.5,pch=17,col="green")
points(frac,rDen,type="b",cex=2,pch=18,col="purple")
points(frac,rMESGA,type="b",cex=1.5,pch=7,col="brown")
points(frac,re_test,type="b",cex=1.5,pch=10,col="orange")

####
legend(x=0.25,y=0.4,legend =
         c("MEScan","MEGSA","Dendrix","TiMEx","WExT","CoMEt"),col=c("red","brown","purple","blue","green","orange"),
       lty=1,pch = c(15,7,18,16,17,10))



rTg = c(83,100,100,100)/100
rTimx = c(0,54,98,100)/100
rCoMEt = c(0,29,41,43)/100
rMESGA =c(5,54,70,87)/100
rDen = c(0,0,0,0)/100
re_test = c( 4,  30,  33  ,34)/100
########
plot(frac,rTg,type="b",cex=1.5,pch=15, col="red",
     ylim=c(0,1.05),xlab="coverage",ylab="power",
     cex.lab=1.5, cex.axis=1, cex.main=1.5, cex.sub=1.5)
points(frac,rTimx,type="b",cex=1.5,pch=16,col="blue")
points(frac,rCoMEt,type="b",cex=1.5,pch=17,col="green")
points(frac,rDen,type="b",cex=2,pch=18,col="purple")
points(frac,rMESGA,type="b",cex=1.5,pch=7,col="brown")
points(frac,re_test,type="b",cex=1.5,pch=10,col="orange")
########################### imblanced, TP53
rTg = c(65,97,100,100)/100
rTimx = c(0,50,94,100)/100
rCoMEt = c(1,27,38,42)/100
rMESGA =c(2,67,94,100)/100
rDen = c(0,0,0,0)/100
re_test = c( 2,  32 , 32  ,34)/100
########
plot(frac,rTg,type="b",cex=1.5,pch=15, col="red",
     ylim=c(0,1.05),xlab="coverage",ylab="power",
     cex.lab=1.5, cex.axis=1, cex.main=1.5, cex.sub=1.5)
points(frac,rTimx,type="b",cex=1.5,pch=16,col="blue")
points(frac,rCoMEt,type="b",cex=1.5,pch=17,col="green")
points(frac,rDen,type="b",cex=2,pch=18,col="purple")
points(frac,rMESGA,type="b",cex=1.5,pch=7,col="brown")
points(frac,re_test,type="b",cex=1.5,pch=10,col="orange")
dev.off()

