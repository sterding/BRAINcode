#Rscript to fit distribution of random trancriptional signal
# Author: Xianjun Dong
# Date: Oct 21, 2015
# Version: 2.0

# input as (N, value)
#df=read.table("transcriptional.noise.rpm.txt", comment.char = "", nrows = 1000000)  # read 1 million lines
#df=log10(as.numeric(do.call('c',apply(df, 1, function(x) rep(x[2],x[1])))))

# # input as value at single nt
df=read.table("transcriptional.noise.rpm.txt", comment.char = "")  # read 1 million lines
df=log10(df[,1])

library(fitdistrplus) # install.packages('fitdistrplus')
fitn=fitdist(df,'norm')
pdf("transcriptional.noise.distribution.pdf", width=8, height=6)
hist(df, breaks=100, prob=TRUE, xlab='log10(RPM)', main='Distribution of transcriptional noise')
lines(density(df, bw=0.15))
m=round(as.numeric(fitn$estimate[1]),digits=3)
sd=round(as.numeric(fitn$estimate[2]),digits=3)
lines(density(rnorm(n=2000000, mean=m, sd=sd),bw=0.25), col='blue',lty=2)
p=round(qnorm(.05, mean=m, sd=sd, lower.tail = F), digits=3)
lines(y=c(0,0.3),x=c(p,p),col='red')
text(p,0.2,paste0("P(X>",p,") = 0.05\nRPM=10**",p,"=",round(10**p,digits=3)), adj=c(0,0))
legend("topright", c("empirical density curve", paste0("fitted normal distribution \n(mean=",m,", sd=",sd,")")), col=c('black','blue'), lty=c(1,2), bty='n')
dev.off()
write.table(10**p, "transcriptional.noise.rpm.pvalues.txt", quote=F, row.names=F, col.names=F)
