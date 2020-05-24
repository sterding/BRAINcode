###########################################
# A general framework to run WGCNA in R
# author: Xianjun Dong
# email: xdong@rics.bwh.harvard.edu
# date: 9/27/2018
# version: 1.0
# Usage: Rscript $0 -i genes.htseqcount.cufflinks.allSamples.uniq.xls -c covaraite.txt
# Note: the input data is better normalized, adjusted, and quantile normalized (see point #4 of https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/faq.html)
###########################################
library("optparse")

option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, 
              help="Input expression dataset file name", metavar="character"),
  make_option(c("-c", "--covariate"), type="character", default=NULL, 
              help="Covariate file name", metavar="character"),
  make_option(c("-a", "--annotation"), type="character", default=NULL, 
              help="Gene annotation file name", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="./", 
              help="output directory path [default=%default]", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$input) || is.null(opt$covariate)){
  print_help(opt_parser)
  stop("Both input expression matrix and covariate file must be supplied", call.=FALSE)
}

input_expression_filename=opt$input
input_covariance_filename=opt$covariate
input_annotation_filename=opt$annotation
output_dir=opt$out

# debug
# setwd("~/neurogen/rnaseq_Rot/results/merged"); input_expression_filename="genes.htseqcount.cufflinks.allSamples.uniq.xls"; input_covariance_filename="covariate.txt"; output_dir="WGCNA"; input_annotation_filename="~/neurogen/referenceGenome/Rattus_norvegicus/UCSC/rn6/Annotation/Genes/gencode.annotation.genes.bed+2"

# # Only run the following commands once to install WGCNA and flashClust on your computer
# source("http://bioconductor.org/biocLite.R")
# biocLite("WGCNA")
# install.packages("flashClust") 

# Load WGCNA and flashClust libraries every time you open R
library(WGCNA)
library(flashClust)

# ----------------------------------------------------------------
# Uploading data into R and formatting it for WGCNA --------------
# ----------------------------------------------------------------

if(tolower(tools::file_ext(input_expression_filename)) == "rds") {
  datExpr=readRDS(file.path(input_expression_filename))
}else datExpr=read.delim(file.path(input_expression_filename), header = T, row.names = 1,check.names =F)

head(datExpr); dim(datExpr);

datExpr = t(datExpr) # now samples are rows and genes are columns

if(nrow(datExpr)<15) warning("We do not recommend attempting WGCNA on a data set consisting of fewer than 15 samples.")

# Run this to check if there are gene outliers. 
gsg = goodSamplesGenes(datExpr, verbose = 3)
# remove the offending genes and samples from the data with the following:
if(!gsg$allOK){
  if(sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr)[!gsg$goodGenes], collapse= ", ")));
	if(sum(!gsg$goodSamples)>0)
		printFlush(paste("Removing samples:", paste(rownames(datExpr)[!gsg$goodSamples], collapse=", ")))
	datExpr= datExpr[gsg$goodSamples, gsg$goodGenes]
}

dim(datExpr)

#Create an object called "datTraits" that contains your trait data
datTraits = read.table(input_covariance_filename, header = T, sep="\t", row.names = 1, check.names=F)
head(datTraits); dim(datTraits);

#form a data frame analogous to expression data that will hold the clinical traits.
all(rownames(datTraits) %in% rownames(datExpr))
dim(datExpr); datExpr = datExpr[intersect(rownames(datExpr), rownames(datTraits)), , drop=F]; dim(datExpr)
dim(datTraits); datTraits = datTraits[rownames(datExpr), , drop=F]; dim(datTraits);
datTraits[] <- lapply(datTraits, function(x) if(is.factor(x)) factor(x) else x) # drop levels  after subsetting
all(rownames(datTraits) == rownames(datExpr))

if(tolower(tools::file_ext(input_annotation_filename)) == "rds") {
  datAnnotation=readRDS(file.path(input_annotation_filename))
}else datAnnotation=read.delim(file.path(input_annotation_filename), header = F, stringsAsFactors =F, col.names = c("chr","start","end","ID","score",'strand','geneName','geneType'))

head(datAnnotation)

# You have finished uploading and formatting expression and trait data
# Expression data is in datExpr, corresponding traits are datTraits

dir.create(file.path(output_dir), recursive =T, showWarnings = FALSE)
setwd(output_dir)
save(datExpr, datTraits, file="SamplesAndTraits.RData")
#load("SamplesAndTraits.RData")

pdf("WGCNA.diagnosis.plots.pdf", paper = 'a4r')
# ----------------------------------------------------------------
# Cluster samples by expression 
# ----------------------------------------------------------------
A = adjacency(t(datExpr),type="signed") # this calculates the whole network connectivity
k = as.numeric(apply(A,2,sum))-1 # standardized connectivity
Z.k = scale(k)
thresholdZ.k = -2.5 # often -2.5
outlierColor = ifelse(Z.k<thresholdZ.k,"red","black")
sampleTree = flashClust(as.dist(1-A), method = "average")
# Convert traits to a color representation where red indicates high values
datTraits2number=datTraits;
datTraits2number[] <- lapply(datTraits2number, function(x) if(is.factor(x)) as.numeric(x)-1 else x) 
traitColors = numbers2colors(datTraits2number) 

datColors = data.frame(outlier = outlierColor,traitColors)

plotDendroAndColors(sampleTree,groupLabels=names(datColors),
                    colors=datColors,main="Sample Dendrogram and Trait Heatmap")

# Now outliers have been identified. You could exclude these samples with the code below, but this scientists had a biological reason to NOT exclude these samples. It's up to you. Justify whatever decision you make.

# Remove outlying samples 
#remove.samples = Z.k<thresholdZ.k | is.na(Z.k)
#datExprOut = datExpr[!remove.samples,]
#datTraitsOut = datTraits[!remove.samples,]
#save(datExprOut, datTraitsOut, file="SamplesAndTraits_OutliersRemoved.RData")


# ----------------------------------------------------------------
# Cluster genes by expression 
# ----------------------------------------------------------------
# Choose a soft threshold power- 
allowWGCNAThreads()  # Not enableWGCNAThreads in RStudio, see note https://www.biostars.org/p/122349/
powers = c(c(1:10), seq(from =11, to=30, by=1)) #choosing a set of soft-thresholding powers
sft = pickSoftThreshold(datExpr, powerVector=powers, verbose =5, networkType="signed") #call network topology analysis function
SFT_r2=-sign(sft$fitIndices[,3])*sft$fitIndices[,2];

#sizeGrWindow(9,5)
par(mfrow= c(1,2))
cex1=0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], xlab= "Soft Threshold (power)", ylab="Scale Free Topology Model Fit, signed R^2", type= "n", main= paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels=powers, cex=cex1, col="red")
abline(h=0.90, col="red")
#from this plot, we would choose the lowest power for which the scale free topology index reaches 0.90
softPower = min(which(SFT_r2>=0.9))
text(softPower, 0.9, labels=softPower, cex=cex1, col="black")

plot(sft$fitIndices[,1], sft$fitIndices[,5], xlab= "Soft Threshold (power)", ylab="Mean Connectivity", type="n", main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1, col="red")

# ----------------------------------------------------------------
# Construct Networks                ------------------------------
# ----------------------------------------------------------------
allowWGCNAThreads()
adjacency = adjacency(datExpr, power = softPower, type = "signed") #specify network type

#translate the adjacency into TOM (topological overlap matrix) and calculate the corresponding dissimilarity:
TOM = TOMsimilarity(adjacency, TOMType="signed") # specify network type
dissTOM = 1-TOM

# ----------------------------------------------------------------
# Generate Modules
# ----------------------------------------------------------------

# Generate a clustered gene tree
geneTree = flashClust(as.dist(dissTOM), method="average")

plot(geneTree, xlab="", sub="", main= "Gene Clustering on TOM-based dissimilarity", labels= FALSE, hang=0.04)

#This sets the minimum number of genes to cluster into a module
minModuleSize = 30 
dynamicMods = cutreeDynamic(dendro= geneTree, distM= dissTOM, deepSplit=2, pamRespectsDendro= FALSE, minClusterSize = minModuleSize)

dynamicColors= labels2colors(dynamicMods)
MEList= moduleEigengenes(datExpr, colors= dynamicColors,softPower = softPower)
MEs= MEList$eigengenes
MEDiss= 1-cor(MEs)
METree= flashClust(as.dist(MEDiss), method= "average")

save(dynamicMods, MEList, MEs, MEDiss, METree, file= "Network_allSamples_signed_RLDfiltered.RData")

#plots tree showing how the eigengenes cluster together
plot(METree, main= "Clustering of module eigengenes", xlab= "", sub= "")

#set a threhold for merging modules. 
MEDissThres = 0.4 # meaning that modules with at least .6 correlation are merged
merge = mergeCloseModules(datExpr, dynamicColors, cutHeight= MEDissThres, verbose =3)
mergedColors = merge$colors
mergedMEs = merge$newMEs

#plot dendrogram with module colors below it
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors), c("Dynamic Tree Cut", "Merged dynamic"), dendroLabels= FALSE, hang=0.03, addGuide= TRUE, guideHang=0.05)

moduleColors = mergedColors
colorOrder = c("grey", standardColors())
moduleLabels = match(moduleColors, colorOrder)-1
MEs = mergedMEs

### this include eigengene for each module, module labels, color labels for the genes correponding to the modules, and gene tree
save(MEs, moduleLabels, moduleColors, geneTree, file= "Network_allSamples_signed_nomerge_RLDfiltered.RData")

# ----------------------------------------------------------------
## save genes in each module into individual files
# ----------------------------------------------------------------
genenames = colnames(datExpr)

library(foreach)
library(doParallel)
cores=detectCores()
cl <- makeCluster(cores[1]-1) #not to overload your computer
registerDoParallel(cl)

# for(module in unique(moduleColors)) {
system.time({foreach(module=unique(moduleColors)) %dopar% {
  message(paste("# save data for module = ",module))
  # Select module probes
  inModule = (moduleColors==module);
  modGenes = genenames[inModule];
  
  # save expression data
  write.table(datExpr[,inModule], file = paste("module", module, "datExpr.tab", sep="."), sep="\t", quote=F, row.names=F, col.names=T)
  
  # save gene symbol 
  write.table(data.frame(ID=modGenes, geneName=datAnnotation$geneName[match(modGenes,datAnnotation$ID)]), 
              file = paste("module", module, "anno.tab", sep="."), sep="\t", quote=F, row.names=F, col.names=T)
  
  # Select the corresponding TOM
  modTOM = TOM[inModule, inModule];
  dimnames(modTOM) = list(modGenes, modGenes)
  # Export the network into an edge list file VisANT can read
  vis = exportNetworkToVisANT(modTOM,
                              file = paste("module", module, "VisANTInput.txt", sep="."),
                              weighted = TRUE,
                              threshold = 0,
                              probeToGene = data.frame(datAnnotation$ID, datAnnotation$geneName) )
  
  ## or only the top 30 hub genes
  # nTop = 30;
  # IMConn = softConnectivity(datExpr[, modGenes]);
  # top = (rank(-IMConn) <= nTop)
  # vis = exportNetworkToVisANT(modTOM[top, top],
  #                             file = paste("module.", module, ".VisANTInputTop30.txt", sep=""),
  #                             weighted = TRUE,
  #                             threshold = 0,
  #                             probeToGene = data.frame(datAnnotation$ID, datAnnotation$geneName) )
  
  
}
})
# ----------------------------------------------------------------
# Correlate traits
# ----------------------------------------------------------------
#Define number of genes and samples
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

#Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits2number, use= "p")  # <=== FIX THIS https://www.biostars.org/p/293281/
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

#Print correlation heatmap between modules and traits
textMatrix= paste(signif(moduleTraitCor, 2), "\n(", 
						signif(moduleTraitPvalue, 1), ")", sep= "")
dim(textMatrix)= dim(moduleTraitCor)
#display the corelation values with a heatmap plot
labeledHeatmap(Matrix= moduleTraitCor, 
			xLabels= names(datTraits), 
			yLabels= names(MEs), 
			ySymbols= names(MEs), 
			colorLabels= FALSE, 
			colors= blueWhiteRed(50), 
			textMatrix= textMatrix, 
			setStdMargins= TRUE, 
			cex.text= 0.5, 
			cex.lab = 0.5,
			zlim= c(-1,1), 
			main= paste("Module-trait relationships"))

# turn off plotting
dev.off()
