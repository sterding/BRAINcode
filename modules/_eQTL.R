###########################################
# Rscript to run eQTL analysis using Matrix eQTL
# Matrix eQTL by Andrey A. Shabalin
# http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/
# Usage: Rscript $PIPELINE_PATH/_eQTL.R snp.txt expression.txt cov.txt output.txt geneloc.txt snploc.txt
# Author: Xianjun Dong
# Version: 1.0
# Date: 2014-May-19
###########################################

## TODO: incooprate the genotype-covairance interaction!!
## TODO: re-run without 
## Note: PD could potentially introduce noise (e.g. those survived neurons could be more resisitent than other died)

require('MatrixEQTL') || install.package('MatrixEQTL', repo='http://cran.revolutionanalytics.com');

args<-commandArgs(TRUE)

SNP_file_name=args[1]
expression_file_name=args[2]
covariates_file_name=args[3]
output_file_name=args[4]
gene_location_file_name = args[5]
snp_location_file_name = args[6]

useModel = modelLINEAR  # modelANOVA, modelLINEAR, or modelLINEAR_CROSS

## Settings
# Only associations significant at this level will be saved
pvOutputThreshold = 5e-3;

# Only associations significant at this level will be saved
pvOutputThreshold_cis = 2e-2;
pvOutputThreshold_tra = 1e-2;

# Distance for local gene-SNP pairs
cisDist = 1e6;

# Error covariance matrix
# Set to numeric() for identity.
errorCovariance = numeric();
# errorCovariance = read.table("Sample_Data/errorCovariance.txt");

message("## Load genotype data...")

snps = SlicedData$new();
snps$fileDelimiter = "\t";      # the TAB character
snps$fileOmitCharacters = "NA"; # denote missing values;
snps$fileSkipRows = 1;          # one row of column labels
snps$fileSkipColumns = 1;       # one column of row labels
snps$fileSliceSize = 2000;      # read file in slices of 2,000 rows
snps$LoadFile(SNP_file_name);

message("## Load gene expression data ...")

gene = SlicedData$new();
gene$fileDelimiter = "\t";      # the TAB character
gene$fileOmitCharacters = "NA"; # denote missing values;
gene$fileSkipRows = 1;          # one row of column labels
gene$fileSkipColumns = 1;       # one column of row labels
gene$fileSliceSize = 2000;      # read file in slices of 2,000 rows
gene$LoadFile(expression_file_name);

message("## Load covariates...")

cvrt = SlicedData$new();
cvrt$fileDelimiter = "\t";      # the TAB character
cvrt$fileOmitCharacters = "NA"; # denote missing values;
cvrt$fileSkipRows = 1;          # one row of column labels
cvrt$fileSkipColumns = 1;       # one column of row labels
cvrt$fileSliceSize = 2000;       # read file in slices of 2,000 rows
if(length(covariates_file_name)>0) {
cvrt$LoadFile(covariates_file_name);
}

if(gene_location_file_name != "" && snp_location_file_name!="")
{
    message("## Load gene/SNP location data...")
    snpspos = read.table(snp_location_file_name, header = TRUE, stringsAsFactors = FALSE);
    genepos = read.table(gene_location_file_name, header = TRUE, stringsAsFactors = FALSE);

    message("## Run the cis-/trans-eQTL analysis...")
    
    me = Matrix_eQTL_main(
    snps = snps, 
    gene = gene, 
    cvrt = cvrt,
    output_file_name     = paste(output_file_name,"trans.txt", sep="."),
    pvOutputThreshold     = pvOutputThreshold_tra,
    useModel = useModel, 
    errorCovariance = errorCovariance, 
    verbose = TRUE, 
    output_file_name.cis = paste(output_file_name,"cis.txt", sep="."),
    pvOutputThreshold.cis = pvOutputThreshold_cis,
    snpspos = snpspos, 
    genepos = genepos,
    cisDist = cisDist,
    pvalue.hist = "qqplot",
    min.pv.by.genesnp = FALSE,
    noFDRsaveMemory = FALSE);
    
    cat('Analysis done in: ', me$time.in.sec, ' seconds', '\n');
    cat('Detected local eQTLs:', '\n');
    show(me$cis$eqtls)
    cat('Detected distant eQTLs:', '\n');
    show(me$trans$eqtls);
} else {
    message("## Run the eQTL analysis...")

    me = Matrix_eQTL_engine(
    snps = snps,
    gene = gene,
    cvrt = cvrt,
    output_file_name = paste(output_file_name,"txt", sep="."),
    pvOutputThreshold = pvOutputThreshold,
    useModel = useModel, 
    errorCovariance = errorCovariance, 
    verbose = TRUE,
    pvalue.hist = FALSE,
    min.pv.by.genesnp = FALSE,
    noFDRsaveMemory = FALSE);
    
    
    message("## Getting results ...")

    cat('Analysis done in: ', me$time.in.sec, ' seconds', '\n');
    cat('Detected eQTLs:', '\n');
    show(me$all$eqtls);

}
pdf(paste(output_file_name, "pdf", sep="."));
plot(me);
dev.off();