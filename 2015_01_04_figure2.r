library(chron)

## input raw data


ph1 <- read.csv("C:\\Users\\matklab\\Google Drive\\kiddivax\\non_bracketing\\data\\y1_ph1.csv",header=F)



######################################################################################################
link <- "C:\\Users\\matklab\\Google Drive\\kiddivax\\non_bracketing\\matlab\\y1_ph1_509\\"

y1_ph1_1 <- read.csv(paste(link,"mcmc_result_1.csv",sep=""),header=F)
y1_ph1_2 <- read.csv(paste(link,"mcmc_result_2.csv",sep=""),header=F)
y1_ph1_3 <- read.csv(paste(link,"mcmc_result_3.csv",sep=""),header=F)
y1_ph1_4 <- read.csv(paste(link,"mcmc_result_4.csv",sep=""),header=F)
y1_ph1_5 <- read.csv(paste(link,"mcmc_result_5.csv",sep=""),header=F)

y1_ph1 <- read.csv(paste(link,"y1_ph1.csv",sep=""),header=F)

ILI_ph1 <- read.csv(paste(link,"ILILAB_ph1.csv",sep=""),header=F)

inc <- 1001:11000



########################################################################################################################################################################


model_estimate <- function(mcmc1,mcmc2,p1_pre,p2_pre,p1,p2){

mcmc_result_1 <- matrix(NA,length(inc),10)
for ( j in 1:10 ){
mcmc_result_1[,j] <- (1-exp(-(mcmc1[,1]*p1+mcmc1[,4]*p2)*exp(mcmc1[,7]*(j-1))))*mcmc2[,j]
}

mcmc_result_2 <- matrix(NA,length(inc),10)
for ( j in 1:10 ){
mcmc_result_2[,j] <- (1-exp(-(mcmc1[,2]*p1+mcmc1[,5]*p2)*exp(mcmc1[,7]*(j-1))))*mcmc2[,j+10]
}

mcmc_result_3 <- matrix(NA,length(inc),10)
for ( j in 1:10 ){
mcmc_result_3[,j] <- (1-exp(-(mcmc1[,3]*p1+mcmc1[,6]*p2)*exp(mcmc1[,7]*(j-1))))*mcmc2[,j+10]
}


mcmc_result_1_pre <- matrix(NA,length(inc),10)
for ( j in 1:10 ){
mcmc_result_1_pre[,j] <- (exp(-(mcmc1[,1]*p1_pre+mcmc1[,4]*p2_pre)*exp(mcmc1[,7]*(j-1))))
}

mcmc_result_2_pre <- matrix(NA,length(inc),10)
for ( j in 1:10 ){
mcmc_result_2_pre[,j] <- (exp(-(mcmc1[,2]*p1_pre+mcmc1[,5]*p2_pre)*exp(mcmc1[,7]*(j-1))))
}

mcmc_result_3_pre <- matrix(NA,length(inc),10)
for ( j in 1:10 ){
mcmc_result_3_pre[,j] <- (exp(-(mcmc1[,3]*p1_pre+mcmc1[,6]*p2_pre)*exp(mcmc1[,7]*(j-1))))
}


mat <- matrix(NA,3,3)
mat[1,] <- quantile(rowSums(mcmc_result_1_pre*mcmc_result_1),c(0.5,0.025,0.975))
mat[2,] <- quantile(rowSums(mcmc_result_2_pre*mcmc_result_2),c(0.5,0.025,0.975))
mat[3,] <- quantile(rowSums(mcmc_result_3_pre*mcmc_result_3),c(0.5,0.025,0.975))


return(c(mat[1,],mat[2,],mat[3,]))
}






########################################################################################################################################################################
## for H1N1pdm09

chpt <- 509

## from 2009/7/5 to 2010/1/16

p1_pre <- sum(ILI_ph1[288:355,])/sum(ILI_ph1)
p2_pre <- 0

p1 <- sum(ILI_ph1[356:chpt,])/sum(ILI_ph1)
p2 <- sum(ILI_ph1[(chpt+1):565,])/sum(ILI_ph1)

mcmc1 <- y1_ph1_1[inc,]
mcmc2 <- y1_ph1_2[inc,]

overall_risk <- model_estimate(mcmc1,mcmc2,p1_pre,p2_pre,p1,p2)



########

table <- matrix(0,10,6)
table[,1] <- dates( paste(c(6:12,1:3), "/01/",c(rep(2009,7),rep(2010,3)),sep="") ) - dates("06/30/2008") - 14
table[,2] <- dates( paste(c(6:12,1:3), "/",c(30,31,31,30,31,30,31,31,28,31),"/",c(rep(2009,7),rep(2010,3)),sep="") ) - dates("06/30/2008") - 14
for ( i in 1:10 ){
table[i,3] <- sum(ILI_ph1[288:min(chpt,table[i,1]),])/sum(ILI_ph1)
if ( table[i,1] > chpt ){
table[i,4] <- sum(ILI_ph1[(chpt+1):table[i,1],])/sum(ILI_ph1)
}
if ( table[i,2] <= chpt ){
table[i,5] <- sum(ILI_ph1[table[i,1]:table[i,2],])/sum(ILI_ph1)
}
if ( table[i,1] > chpt ){
table[i,6] <- sum(ILI_ph1[table[i,1]:table[i,2],])/sum(ILI_ph1)
}
if ( table[i,1] < chpt & table[i,2] > chpt ){
table[i,5] <- sum(ILI_ph1[table[i,1]:chpt,])/sum(ILI_ph1)
table[i,6] <- sum(ILI_ph1[(chpt+1):table[i,2],])/sum(ILI_ph1)
}
}
table2 <- matrix(NA,10,9)
for ( i in 1:10){
table2[i,] <- model_estimate(mcmc1,mcmc2,table[i,3],table[i,4],table[i,5],table[i,6])
}



#########################################
## plot figure 2

pdf("C:\\Users\\matklab\\Google Drive\\kiddivax\\non_bracketing\\summary\\figure2_1.pdf",width=9,height=3)

layout(matrix( 1:2, nrow=1),widths=c(5,3))

par(mar=c(2,4,1,0))



plot(0,0,xlab="",ylab="", main="", axes=F, xlim=c(0,10), ylim=c(0,0.16),type="n")
axis(1,at=c(1:9)-0.5,labels=c("Jul 09",NA,"Sep 09",NA,"Nov 09",NA,"Jan 10",NA,"Mar 10"),cex.axis=0.95)
axis(2,at=0.04*0:4, las=1, pos=0)

points(1:9-0.65,table2[2:10,1],cex=0.8,pch=16,col="black")
points(1:9-0.5,table2[2:10,4],cex=0.8,pch=17,col="black")
points(1:9-0.35,table2[2:10,7],cex=0.8,pch=18,col="black")

for ( i in 1:9){
lines(rep(i-0.65,2),table2[i+1,2:3],col="black")
lines(rep(i-0.5,2),table2[i+1,5:6],col="black")
lines(rep(i-0.35,2),table2[i+1,8:9],col="black")
}
mtext("Infection risk", side=2, line=2.5)

legend(5,0.16, cex=0.9, legend=c("Children","Adults","Older adults"), pch=c(16,17,18), lty=c(1,1,1),col=c("black","black","black"))


title(main="A", adj=0)




par(mar = c(2,1,1,1) )


plot(NA,xlim=c(0,3),ylim=c(0,0.5),axes=FALSE,xlab="",ylab="",main="")
axis(1,at=c(-1,0:2+0.5,4),labels=c(NA,"Children","Adults","Older adults",NA))
axis(2,at=0:5*0.1, las=1)

points(0:2+0.5,overall_risk[c(1,4,7)],cex=0.6,pch=16)

lines(rep(0.5,2),overall_risk[c(2,3)])
lines(rep(1.5,2),overall_risk[c(5,6)])
lines(rep(2.5,2),overall_risk[c(8,9)])

mtext("Cumulative incidence", side=2, line=2.5)

title(main="B", adj=0)

dev.off()




