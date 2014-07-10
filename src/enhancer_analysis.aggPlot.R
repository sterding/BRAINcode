## script to analysis enhancer regions


#============================================================
# draw aggregation plot (to run in R console)
#============================================================
output_dir="~/neurogen/rnaseq_PD/results/eRNA/externalData";
pdf("aggregation.enhancers.SN.pdf", width=18, height=16)
par(mfcol=c(8,6))
for(mark in c("CAGE","Histone","TFBS","Conservation","DNase","RNAseq")){
    print(mark);
    
    n=nrow(read.table(paste(output_dir, "CAGE", paste("CAGE.fwd.bigwig",mark,"bins",sep="."), sep="/"), header=F))
    
    #CAGE
    df1=read.table(paste(output_dir, "CAGE", paste("CAGE.fwd.bigwig",mark,"bins",sep="."), sep="/"), header=F)[,-1]
    df2=read.table(paste(output_dir, "CAGE", paste("CAGE.rev.bigwig",mark,"bins",sep="."), sep="/"), header=F)[,-1]
    df=cbind(df1, df2);
    df=df[apply(df,1,sum)>0,];
    m=nrow(df);
    df=apply(df, 2, function(x) mean(x, trim=0.001)); N=length(df)/2
    o=par(mar = c(2,4,4,1))
    plot(df[1:N], type="l", col="green", ylab=NA, xlab=NA, xaxt="n", main=paste(mark, "-defined enhancers\n(N=",formatC(n,format="d",big.mark=","),")",sep=""), ylim=range(df)) 
    points(df[-c(1:N)], type="l", col="blue")
    if(mark=="CAGE")  {
        legend("topright", c("CAGE +","CAGE -"), col=c("green","blue"), lty=1, bty='n')
        mtext(side = 2, "Mean signal of CAGE \nin Substantia Nigra", line = 2, cex=0.6)
    }
    axis(1, at=50*c(0:4), labels=c("-500 bp","","0","","500 bp"))
    legend("topleft", paste("n=", formatC(m,format="d",big.mark=","), sep=''), bty='n')
    
    #Histone of SN
    df1=read.table(paste(output_dir, "Histone", paste("Histone.SN.H3K27ac.bigwig",mark,"bins",sep="."), sep="/"), header=F)[,-1]
    df2=read.table(paste(output_dir, "Histone", paste("Histone.SN.H3K4me1.bigwig",mark,"bins",sep="."), sep="/"), header=F)[,-1]
    df3=read.table(paste(output_dir, "Histone", paste("Histone.SN.H3K4me3.bigwig",mark,"bins",sep="."), sep="/"), header=F)[,-1]
    df4=read.table(paste(output_dir, "Histone", paste("Histone.SN.H3K27me3.bigwig",mark,"bins",sep="."), sep="/"), header=F)[,-1]
    df=cbind(df1, df2, df3, df4)
    df=df[apply(df,1,sum)>0,];
    m=nrow(df);
    df=apply(df, 2, function(x) mean(x, trim=0.01)); N=length(df)/4
    par(mar = c(2,4,1,1))
    plot(df[1:N], type="l", col="green", ylab=NA, xlab=NA, xaxt="n", ylim=range(df)) 
    points(df[N+c(1:N)], type="l", col="blue")
    points(df[2*N+c(1:N)], type="l", col="gold")
    points(df[3*N+c(1:N)], type="l", col="black")
    if(mark=="CAGE") {
        legend("topright", c("H3K27ac","H3K4me1", "H3K4me3","H3K27me3"), col=c("green","blue", "gold","black"), lty=1, bty='n');
        mtext(side = 2, "Mean signal of histone marks \nin Substantia Nigra", line = 2, cex=0.6)
    }
    axis(1, at=50*c(0:4), labels=c("-500 bp","","0","","500 bp"))
    legend("topleft", paste("n=", formatC(m,format="d",big.mark=","), sep=''), bty='n')

    #Histone of HeLaS3 cell
    df1=read.table(paste(output_dir, "Histone", paste("Histone.Helas3.H3K27ac.bigwig",mark,"bins",sep="."), sep="/"), header=F)[,-1]
    df2=read.table(paste(output_dir, "Histone", paste("Histone.Helas3.H3K4me1.bigwig",mark,"bins",sep="."), sep="/"), header=F)[,-1]
    df3=read.table(paste(output_dir, "Histone", paste("Histone.Helas3.H3K4me3.bigwig",mark,"bins",sep="."), sep="/"), header=F)[,-1]
    df4=read.table(paste(output_dir, "Histone", paste("Histone.Helas3.H3K27me3.bigwig",mark,"bins",sep="."), sep="/"), header=F)[,-1]
    df=cbind(df1, df2, df3, df4)
    df=df[apply(df,1,sum)>0,];
    m=nrow(df);
    df=apply(df, 2, function(x) mean(x, trim=0.01)); N=length(df)/4
    plot(df[1:N], type="l", col="green", ylab=NA, xlab=NA, xaxt="n", ylim=range(df, 10)) 
    points(df[N+c(1:N)], type="l", col="blue")
    points(df[2*N+c(1:N)], type="l", col="gold")
    points(df[3*N+c(1:N)], type="l", col="black")
    if(mark=="CAGE") {
        legend("topright", c("H3K27ac","H3K4me1", "H3K4me3","H3K27me3"), col=c("green","blue", "gold","black"), lty=1, bty='n');
        mtext(side = 2, "Mean signal of histone marks \nin HelaS3 cell", line = 2, cex=0.6)
    }
    axis(1, at=50*c(0:4), labels=c("-500 bp","","0","","500 bp"))
    legend("topleft", paste("n=", formatC(m,format="d",big.mark=","), sep=''), bty='n')
    
    #TFBS
    df=read.table(paste(output_dir, "TFBS", paste("TFBS.bigwig",mark,"bins",sep="."), sep="/"), header=F)[,-1]
    df=df[apply(df,1,sum)>0,];
    m=nrow(df);
    df=apply(df, 2, function(x) mean(x, trim=0.01));
    plot(df, type="l", col="blue", ylab=NA, xlab=NA, xaxt="n") 
    axis(1, at=50*c(0:4), labels=c("-500 bp","","0","","500 bp"))
    legend("topleft", paste("n=", formatC(m,format="d",big.mark=","), sep=''), bty='n')
    if(mark=="CAGE") mtext(side = 2, "Mean TF binding occurances \nin 161 ENCODE TFs", line = 2, cex=0.6)
    
    #Conservation
    df=read.table(paste(output_dir, "Conservation", paste("Conservation.bigwig",mark,"bins",sep="."), sep="/"), header=F)[,-1]
    df=df[apply(df,1,sum)>0,];
    m=nrow(df);
    df=apply(df, 2, function(x) mean(x, trim=0.01));
    plot(df, type="l", col="blue", ylab=NA, xlab=NA, xaxt="n", ylim=range(df,0.2)) 
    axis(1, at=50*c(0:4), labels=c("-500 bp","","0","","500 bp"))
    legend("topleft", paste("n=", formatC(m,format="d",big.mark=","), sep=''), bty='n')
    if(mark=="CAGE") mtext(side = 2, "Mean phyloP score \nin 46-way conservation", line = 2, cex=0.6)

    
    #DNase
    df=read.table(paste(output_dir, "DNase", paste("DNase.bigwig",mark,"bins",sep="."), sep="/"), header=F)[,-1]
    df=df[apply(df,1,sum)>0,];
    m=nrow(df);
    df=apply(df, 2, function(x) mean(x, trim=0.01));
    plot(df, type="l", col="blue", ylab=NA, xlab=NA, xaxt="n") 
    axis(1, at=50*c(0:4), labels=c("-500 bp","","0","","500 bp"))
    legend("topleft", paste("n=", formatC(m,format="d",big.mark=","), sep=''), bty='n')
    if(mark=="CAGE") mtext(side = 2, "Mean DNase signal \nin fetal brain (Roadmap Epigenomics)", line = 2, cex=0.6)

    ##Methylation
    #df=read.table(paste(output_dir, "Methylation", paste("Methylation.bigwig",mark,"bins",sep="."), sep="/"), header=F)[,-1]
    #df=df[apply(df,1,sum)>0,];
    #m=nrow(df);
    #df=apply(df, 2, function(x) mean(x, trim=0.001));
    #plot(df, type="l", col="blue",ylab=NA, xlab=NA, xaxt="n") 
    #axis(1, at=50*c(0:4), labels=c("-500 bp","","0","","500 bp"))
    #legend("topleft", paste("n=", formatC(m,format="d",big.mark=","), sep=''), bty='n')
    #if(mark=="CAGE") mtext(side = 2, "Mean signal of Methylation \nin Substantia Nigra", line = 2, cex=0.6)
    
    #RNAseq
    df=read.table(paste(output_dir, "RNAseq", paste("RNAseq.bigwig",mark,"bins",sep="."), sep="/"), header=F)[,-1]
    df=df[apply(df,1,sum)>0,];
    m=nrow(df);
    df=apply(df, 2, function(x) mean(x, trim=0.01));
    plot(df, type="l", col="blue", ylab=NA, xlab=NA, xaxt="n") 
    axis(1, at=50*c(0:4), labels=c("-500 bp","","0","","500 bp"))
    legend("topleft", paste("n=", formatC(m,format="d",big.mark=","), sep=''), bty='n')
    if(mark=="CAGE") mtext(side = 2, "Mean signal of RNAseq \nin SNDA", line = 2, cex=0.6)

    # intersect with RNAseq
    df=read.table(paste(output_dir, "RNAseq", paste("RNAseq",paste(mark,"_n_RNAseq",sep=""),"bins",sep="."), sep="/"), header=F)[,-1]
    df=df[apply(df,1,sum)>0,];
    m=nrow(df);
    df=apply(df, 2, function(x) mean(x, trim=0.01));
    par(mar = c(4,4,2,1))
    plot(df, type="l", col="red", ylab=NA, xlab="Center position of enhancers", xaxt="n", main=paste(mark,"- & RNAseq-defined enhancers",sep="")) 
    axis(1, at=50*c(0:4), labels=c("-500 bp","","0","","500 bp"))
    legend("topleft", paste("n=", formatC(m,format="d",big.mark=","), sep=''), bty='n')
    if(mark=="CAGE") mtext(side = 2, "Mean signal of RNAseq \nin SNDA", line = 2, cex=0.6)

}
dev.off()

quit('no');

#============================================================
# draw plot for DoD grant  (run in R console)
#============================================================
output_dir="~/neurogen/rnaseq_PD/results/eRNA/externalData";
pdf("aggregation.enhancers.forDoD.pdf", width=6, height=12)
par(mfcol=c(6,2), mar = c(5,4,4,2))
for(mark in c("CAGE","ENCODE")){
    print(mark);
    
    n=nrow(read.table(paste(output_dir, "CAGE", paste("CAGE.fwd.bigwig",mark,"bins",sep="."), sep="/"), header=F))
    
    #CAGE
    df1=apply(read.table(paste(output_dir, "CAGE", paste("CAGE.fwd.bigwig",mark,"bins",sep="."), sep="/"), header=F)[,-1], 2, function(x) mean(x, trim=10/n))
    df2=apply(read.table(paste(output_dir, "CAGE", paste("CAGE.rev.bigwig",mark,"bins",sep="."), sep="/"), header=F)[,-1], 2, function(x) mean(x, trim=10/n))
    plot(df1, type="l", col="green", ylab="Franction of CAGE tag", xlab="", xaxt="n", main=paste(mark, "-defined enhancers\n(N=",formatC(n,format="d",big.mark=","),")",sep=""), ylim=range(df1,df2)) #c(-max(df1,df2),max(df1,df2)))
    points(df2, type="l", col="blue")
    legend("topright", c("CAGE +","CAGE -"), col=c("green","blue"), lty=1, bty='n')
    axis(1, at=50*c(0:4), labels=c("-500","","0","","500"))
    
    #ENCODE
    df1=apply(read.table(paste(output_dir, "ENCODE", paste("ENCODE.H3k27ac.bigwig",mark,"bins",sep="."), sep="/"), header=F)[,-1], 2, function(x) mean(x, trim=0.01))
    df2=apply(read.table(paste(output_dir, "ENCODE", paste("ENCODE.H3k4me1.bigwig",mark,"bins",sep="."), sep="/"), header=F)[,-1], 2, function(x) mean(x, trim=0.01))
    df3=apply(read.table(paste(output_dir, "ENCODE", paste("ENCODE.H3k4me3.bigwig",mark,"bins",sep="."), sep="/"), header=F)[,-1], 2, function(x) mean(x, trim=0.01))
    df4=apply(read.table(paste(output_dir, "ENCODE", paste("ENCODE.H3k27me3.bigwig",mark,"bins",sep="."), sep="/"), header=F)[,-1], 2, function(x) mean(x, trim=0.01))
    plot(df1, type="l", col="green", ylab="Mean signal of histone marks", xlab="", xaxt="n", ylim=range(df1,df2, df3, df4)) 
    points(df2, type="l", col="blue")
    points(df3, type="l", col="gold")
    points(df4, type="l", col="black")
    legend("topright", c("H3K27ac","H3K4me1", "H3K4me3","H3K27me3"), col=c("green","blue", "gold","black"), lty=1, bty='n')
    axis(1, at=50*c(0:4), labels=c("-500","","0","","500"))
    
    #HOT
    df=apply(read.table(paste(output_dir, "HOT", paste("HOT.bigwig",mark,"bins",sep="."), sep="/"), header=F)[,-1], 2, function(x) mean(x, trim=0.01))
    plot(df, type="l", col="green", ylab="Mean TFBS occurances", xlab="", xaxt="n") 
    axis(1, at=50*c(0:4), labels=c("-500","","0","","500"))
    
    #Conservation 
    df=apply(read.table(paste(output_dir, "Conservation", paste("Conservation.bigwig",mark,"bins",sep="."), sep="/"), header=F)[,-1], 2, function(x) mean(x, trim=0.01))
    plot(df, type="l", col="green", ylab="Mean phyloP score", xlab="", xaxt="n") 
    axis(1, at=50*c(0:4), labels=c("-500","","0","","500"))
    
    #DNase
    df=apply(read.table(paste(output_dir, "DNase", paste("DNase.bigwig",mark,"bins",sep="."), sep="/"), header=F)[,-1], 2, function(x) mean(x, trim=0.01))
    plot(df, type="l", col="green", ylab="Mean DNase peak occurances", xlab=paste("Center position of\n",mark, "-defined enhancers",sep=""), xaxt="n") 
    axis(1, at=50*c(0:4), labels=c("-500","","0","","500"))
    
    #RNAseq
    df=read.table(paste(output_dir, "RNAseq", paste("RNAseq.bigwig",paste(mark,"_n_RNAseq",sep=""),"bins",sep="."), sep="/"), header=F)
    n=nrow(df)
    df=apply(df[,-1], 2, function(x) mean(x, trim=0.01))
    plot(df, type="l", col="red", ylab="Mean RNAseq signal", xlab=paste("Center position of\n",mark, "-defined enhancers",sep=""), xaxt="n", main=paste(mark, " enhancers expressed in SNDA\n(N=",formatC(n,format="d",big.mark=","),")",sep="")) 
    axis(1, at=50*c(0:4), labels=c("-500","","0","","500"))
}
        
dev.off()
