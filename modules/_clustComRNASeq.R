##########################################################################################
# This module is to perform sample hierarchical clustering analysis using cufflinks' output
# i.e. genes.fpkm_tracking and isoforms.fpkm_tracking
# Input:  1. the directory of the output of Cufflinks
#         2. the path where the output of cluster analysis will be 
# Output: 1. the excel sheets contain the FPKM values of the genes and transcripts are found in all samples
#         2. the pdf files contains a plot of the hierachical clustering 
# Usage: Rscript clustComRNASeq.R inputPath outputPath
#
##########################################################################################
args <- commandArgs(TRUE)
DsinPath=args[1]
DsoutPath=args[2]
a = list.files(DsinPath)
commGene = as.character(read.table(paste(DsinPath,a[1],"genes.fpkm_tracking",sep="/"), header=T, row.names=NULL)[,1])
commIso = as.character(read.table(paste(DsinPath,a[1],"isoforms.fpkm_tracking",sep="/"), header=T, row.names=NULL)[,1])
for (i in 1:length(a)) {
	temGene = as.character(read.table(paste(DsinPath, a[i],"genes.fpkm_tracking",sep="/"), header=T, row.names=NULL)[,1])
	temIso = as.character(read.table(paste(DsinPath, a[i],"isoforms.fpkm_tracking",sep="/"), header=T, row.names=NULL)[,1])
	commGene = intersect(commGene,temGene)
	commIso = intersect(commIso,temIso)
}
if (length(commGene)==0 || length(commIso)==0)
stop("Error: No common genes or transcripts, the program will abort...")

commGene=sort(commGene)
commIso=sort(commIso)
## extract common genes from 
xGenes=matrix(nrow=length(commGene),ncol=length(a))
xIsos=matrix(nrow=length(commIso),ncol=length(a))

for (i in 1:length(a)) {
	temGenes = read.table(paste(DsinPath, a[i],"genes.fpkm_tracking",sep="/"), header=T, row.names=NULL)
	temGenes = temGenes[which(!duplicated(temGenes[,1])),]
	whcom.gene = which(as.character(temGenes[,1]) %in% commGene)
	temGenes = temGenes[whcom.gene,]
	temGenes = temGenes[order(temGenes[,1]),]
	# cat('identical ', identical(commGene,as.character(temGenes[,1])),'\n')
	xGenes[,i]=temGenes[,10]
	temIsos = read.table(paste(DsinPath, a[i],"isoforms.fpkm_tracking",sep="/"), header=T, row.names=NULL)
	temIsos = temIsos[which(!duplicated(temIsos[,1])),]
	whcom.Iso = which(as.character(temIsos[,1]) %in% commIso)
	temIsos = temIsos[whcom.Iso,]
	temIsos = temIsos[order(temIsos[,1]),]
	xIsos[,i]=temIsos[,10]
}
## gene level
rownames(xGenes)=commGene
colnames(xGenes)=a
scaledxGenes = scale(xGenes)
hcGene <- hclust(dist(t(scaledxGenes)), "ave")
## isoform level
rownames(xIsos)=commIso
colnames(xIsos)=a
scaledxIsos = scale(xIsos)
hcIso <- hclust(dist(t(scaledxIsos)), "ave")

pdf(paste(DsoutPath,"FPKM cluster on the gene level.pdf",sep="/"), width=8,height=8)
plot(hcGene, hang=-1, main="Gene Level", xlab="Samples", sub="", y="Height", cex=0.5)
dev.off()
pdf(paste(DsoutPath, "FPKM cluster on the isoform level.pdf",sep="/"), width=8,height=8)
plot(hcIso, hang=-1, main="Isoform level", xlab="Samples", sub="", y="Height", cex=0.5)
dev.off()

write.table(xGenes,file=paste(DsoutPath, paste("Common ", nrow(xGenes), " Genes FPKM.xls", sep=""), sep="/"), sep="\t",row.names=T)
write.table(xIsos,file=paste(DsoutPath, paste("Common ", nrow(xIsos), " Isoforms FPKM.xls", sep=""), sep="/"), sep="\t",row.names=T)


