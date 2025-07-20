link <- "C:\\Users\\matklab\\Dropbox\\kiddivax\\non_bracketing\\matlab\\"


y1_ph1_1 <- read.csv(paste(link,"y1_ph1_509\\mcmc_result_1.csv",sep=""),header=F)
y1_ph1_2 <- read.csv(paste(link,"y1_ph1_509\\mcmc_result_2.csv",sep=""),header=F)
y1_ph1_3 <- read.csv(paste(link,"y1_ph1_509\\mcmc_result_3.csv",sep=""),header=F)
y1_ph1_4 <- read.csv(paste(link,"y1_ph1_509\\mcmc_result_4.csv",sep=""),header=F)

y1_ph1 <- read.csv(paste(link,"y1_ph1_0\\y1_ph1.csv",sep=""),header=F)

y1_ph1[y1_ph1 == -1] <- NA

ILI_ph1 <- read.csv(paste(link,"y1_ph1_0\\ILILAB_ph1.csv",sep=""),header=F)




inc <- 1001:11000

##############################################################################################


plotfunction <- function(data,mcmc,panel,title){


par(mar = c(4,4,2,1) )

pd <- matrix(NA,8,4)
for ( u in 1:8 ){
pd[u,1] <- sum(data[,10]/data[,8] == 2^(u+1) & !is.na(data[,8]) & !is.na(data[,10]) ) + sum(data[,12]/data[,10] == 2^(u+1) & !is.na(data[,12]) &  !is.na(data[,10]) )
pd[u,2:4] <- quantile(mcmc[,u],c(0.5,0.025,0.975))
}
pd[,1] <- pd[,1] / sum(pd[,1])


plot(NA,xlim=c(0,8),ylim=c(0,0.4),axes=FALSE,xlab="",ylab="",main="")

#points(0:7+0.4,pd[1:8,1],cex=0.7,pch=2)
points(0:7+0.5,pd[1:8,2],cex=0.7,pch=16)
for ( u in 1:8 ){
lines(rep(u-0.5,2),pd[u,3:4])
}

axis(1, at=0:7+0.5, labels=2^(2:9),cex=0.9,cex.axis=0.9)
axis(2, at=0:2*0.2,labels=c("0%","20%","40%"), las=1, pos=0,cex.axis=0.9)

mtext("Proportion",side = 2, line = 2.5)
mtext("Fold rise",side = 1, line = 2.5)
mtext(title,side=3,line=0)

if ( panel == "A"){
legend(5,0.9, legend=c("Observed","Estimated" ),pch=c(2,1))
}
title(main=panel, adj=0)

}


##############################################################################################


plotfunction2 <- function(data,mcmc,panel,title){


par(mar = c(4,4,2,1) )

pd <- matrix(NA,11,4)
for ( u in 1:11 ){
pd[u,1] <- sum(data[,10]/data[,8] == 2^(u-10),na.rm=T )  + sum(data[,12]/data[,10] == 2^(u-10),na.rm=T )  + sum(data[,14]/data[,12] == 2^(u-10),na.rm=T )
pd[u,2:4] <- quantile(mcmc[,u],c(0.5,0.025,0.975))
}
pd[,1] <- pd[,1] / sum(pd[,1])


plot(NA,xlim=c(0,11),ylim=c(0,0.4),axes=FALSE,xlab="",ylab="",main="")

#points(0:10+0.4,pd[1:11,1],cex=0.7,pch=2)
points(0:10+0.5,pd[1:11,2],cex=0.7,pch=16)
for ( u in 1:11 ){
lines(rep(u-0.5,2),pd[u,3:4])
}

axis(1, at=0:10+0.5, labels=2^c(9,NA,7,NA,5,NA,3,NA,1,NA,-1),cex=0.9,cex=0.6)
axis(2, at=0:2*0.2,labels=c("0%","20%","40%"), las=1, pos=0,cex=0.6)

mtext("Proportion",side = 2, line = 2.5)
mtext("Fold drop",side = 1, line = 2.5)
mtext(title,side=3,line=0)

#if ( panel == "A"){
#legend(5,0.9, legend=c("Observed","Estimated" ),pch=c(2,1))
#}
title(main=panel, adj=0)

}


#########################################################
pdf("C:\\Users\\matklab\\Dropbox\\kiddivax\\non_bracketing\\summary\\figure3.pdf",width=9, height=3.5)
layout(matrix( 1:2, nrow=1,byrow=T))


plotfunction(y1_ph1[y1_ph1[,5] <= 18 & !is.na(y1_ph1[,5]),],y1_ph1_3[inc,1:8],"A","Children")
plotfunction(y1_ph1[y1_ph1[,5] > 18 & !is.na(y1_ph1[,5]),],y1_ph1_3[inc,9:16],"B","Adults")
#plotfunction2(y1_ph1[y1_ph1[,5] <= 18 & !is.na(y1_ph1[,5]),],y1_ph1_4[inc,1:11],"C","H1N1pdm09, children")
#plotfunction2(y1_ph1[y1_ph1[,5] > 18 & !is.na(y1_ph1[,5]),],y1_ph1_4[inc,12:22],"D","H1N1pdm09, adults")

dev.off()


#######################################

boosting <- matrix(NA,16,3)
waning <- matrix(NA,22,3)

for ( i in 1:16){
boosting[i,] <- quantile(y1_ph1_3[inc,i],c(0.5,0.025,0.975))
}

1-quantile(y1_ph1_3[inc,1]+y1_ph1_3[inc,2],c(0.5,0.025,0.975))
1-quantile(y1_ph1_3[inc,9]+y1_ph1_3[inc,10],c(0.5,0.025,0.975))

for ( i in 1:22){
waning[i,] <- quantile(y1_ph1_4[inc,i],c(0.5,0.025,0.975))
}


1-quantile(y1_ph1_4[inc,10]+y1_ph1_4[inc,11],c(0.5,0.025,0.975))
1-quantile(y1_ph1_4[inc,21]+y1_ph1_4[inc,22],c(0.5,0.025,0.975))


### expected fold-rise and fold-drop
exp_rise_drop <- matrix(NA,length(inc),4)

for ( i in inc){
  exp_rise_drop[i-1000,1] <- sum(y1_ph1_3[i,1:8]*(2:9))
  exp_rise_drop[i-1000,2] <- sum(y1_ph1_3[i,9:16]*(2:9))
  exp_rise_drop[i-1000,3] <- sum(y1_ph1_4[i,1:11]*(-9:1))
  exp_rise_drop[i-1000,4] <- sum(y1_ph1_4[i,12:22]*(-9:1))
}

2^quantile(exp_rise_drop[,1],c(0.5,0.025,0.975))
2^quantile(exp_rise_drop[,2],c(0.5,0.025,0.975))
1/2^quantile(exp_rise_drop[,3],c(0.5,0.025,0.975))
1/2^quantile(exp_rise_drop[,4],c(0.5,0.025,0.975))

sum(exp_rise_drop[,1]<exp_rise_drop[,2])
sum(exp_rise_drop[,3]>exp_rise_drop[,4])

#################

out <- matrix(NA,11,2)
wan <- round(waning*100,1)
out[,1] <- paste(wan[1:11,1],"% (",wan[1:11,2],"%, ",wan[1:11,3],"%)",sep="")
out[,2] <- paste(wan[12:22,1],"% (",wan[12:22,2],"%, ",wan[12:22,3],"%)",sep="")

write.csv(out[11:1,],"C:\\Users\\matklab\\Dropbox\\kiddivax\\non_bracketing\\summary\\table_s5.csv")