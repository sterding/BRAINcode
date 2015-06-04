###########################################
## R script to calculate the RTC (Regulatory Trait Concordance) for eQTL SNPs
# Usage: Rscript $pipeline_path/modules/_RTC.R SNP.txt expression.txt cis_eQTL_output.txt
# Author: Xianjun Dong
# Version: 0.0
# Date: 2015-Apr-16
###########################################
args<-commandArgs(TRUE)

snps_file_name=args[1]  # in format of http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/Sample_Data/SNP.txt
expr_file_name=args[2]  # in format of http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/Sample_Data/GE.txt
eqtl_file_name=args[3]  # a tab-delimited file with haeder of "SNP     gene    beta    t-stat  p-value FDR"
ld_file_name = args[4]  # a file to tell SNPs in LD interval, in format of "chr	LDstart	LDend	LDID	gwas	SNPs"

# Note: it's possible to add covariate when calculating the residuals. See example code here:
# http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/CompareResults.html
cvrt_file_name=args[5] # in format of http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/Sample_Data/Covariates.txt

## function
get_SNPs_in_LD <- function(tagSNP, LD)
{
    index=grep(paste0("\\b",tagSNP,"\\b"),LD$genotype)[1]  # only the first LD is taken, if the query SNP is found (with exact word match) in multiple LDs
    return(unlist(strsplit(as.character(LD[index,6]),";")))
}

get_GWAS_in_LD <- function(tagSNP, LD)
{
    index=grep(paste0("\\b",tagSNP,"\\b"),LD$genotype)[1]  # only the first LD is taken, if the query SNP is found (with exact word match) in multiple LDs
    return(unlist(strsplit(as.character(LD[index,5]),";")))
}

# ---------------
# snps is a data frame for SNPs in a LD, with each row a SNP and each column a sample
# expr is a vecter of values for a single gene's expression
# id_GWAS is a character string for the ID of GWAS SNP
# id_eQTL is a character string for the ID of eQTL SNP
# ---------------

get_RTC_score <- function(SNPS, EXPR, id_GWAS, id_eQTL)
{ 
    if(id_GWAS == id_eQTL) return(1);
    
    eQTL_SNP = SNPS[id_eQTL, ]
    
    # SNP set except the eQTL: put the GWAS SNP at the first
    SNPS2 = rbind(SNPS[id_GWAS, ], SNPS[!(rownames(SNPS) %in% c(id_GWAS, id_eQTL)), ])
    
    # get pvalue for each SNP in SNPS2
    pvalues=c()
    
    for(i in 1:nrow(SNPS2)){
        df=data.frame(y=as.numeric(EXPR), x=as.numeric(SNPS2[i,]), x0=as.numeric(eQTL_SNP))
        df=na.exclude(df);
        res=resid(lm(y~x, df));
        pvalues=c(pvalues, cor.test(res, df$x0, method="spearman")$p.value);  # Spearman Rank Correlation
    }
    
    # Rank of GWAS SNP in all tested SNPs
    N=nrow(SNPS2);
    Rank_GWAS = which(unique(sort(pvalues)) %in% pvalues[1])  # take the first hit in unique order if multiple SNPs have the same pvalue
    RTC = (N-Rank_GWAS)/N
    return(RTC);
}


# test using Ganqiang's data
snps_file_name="/PHShome/gl871/neorogen/PDBPNeuroX/HBS_Feb21st2015/Matrix-eQTL/RTC/HC/HBS-HCSub.Matrix.uniq.genotype"
expr_file_name="/PHShome/gl871/neorogen/PDBPNeuroX/HBS_Feb21st2015/Matrix-eQTL/RTC/HC/HBS-HCSubTX.Expression.txt"
eqtl_file_name="/PHShome/gl871/neorogen/PDBPNeuroX/HBS_Feb21st2015/Matrix-eQTL/RTC/HC/HBS-HC.eQTL.trans.txt"
ld_file_name="/PHShome/gl871/neorogen/PDBPNeuroX/HBS_Feb21st2015/Matrix-eQTL/RTC/HC/HBS-HC.LD.GWAS.Genotype.txt"
gwas_file_name='/PHShome/gl871/neorogen/PDBPNeuroX/HBS_Feb21st2015/Matrix-eQTL/RTC/HC/HBS-HC.GWAS.Genotype.txt'

message("## loading input files ...")
snps = read.table(snps_file_name, header=T, stringsAsFactors =F); rownames(snps) = snps[,1]; snps=snps[,-1];
expr = read.table(expr_file_name, header=T, stringsAsFactors =F); rownames(expr) = expr[,1]; expr=expr[,-1];
eqtl = read.table(eqtl_file_name, header=T, stringsAsFactors =F); 
LD   = read.table(ld_file_name, header=T, stringsAsFactors =F); 
GWAS = read.delim(gwas_file_name, header=T, stringsAsFactors =F);

result = c();

message("## processing the eQTL gene-snp pairs ...")

for(i in 1:nrow(eqtl)){
    id_eQTL = eqtl[i,1];
    id_gene = eqtl[i,2];
    FDR = eqtl[i,6];
    
    message(paste0(" # record ",i,": ", id_eQTL," -- ", id_gene,"  ..."))
    
    # extreme case: eQTL itself is a GWAS 
    if(id_eQTL %in% GWAS$id){ 
        rtc_score = 1;
        id_GWAS = id_eQTL;
        trait = GWAS[GWAS$genotype==id_GWAS, 'Disease']
        
        result = rbind(result, c(id_eQTL, id_gene, FDR, id_GWAS, trait, rtc_score))
        message(paste("   ", id_GWAS, trait, id_gene, rtc_score))
        
    } else {
        snps_in_LD = get_SNPs_in_LD(tagSNP=id_eQTL, LD=LD)
        index = snps_in_LD %in% GWAS$genotype;  # GWAS_in_LD = get_GWAS_in_LD(tagSNP=id_eQTL, LD=LD)
        
        if(is.na(snps_in_LD) || sum(index)==0) next;  # skip if no LD found or no GWAS in LD
            
        # If >1 GWAS in the LD, take them one by one
        if(sum(index)>1) message("Warning: More than one GWAS found in the LD, take them one by one for downstream analysis!")
        for(id_GWAS in snps_in_LD[which(index)]){
            trait = GWAS[GWAS$genotype==id_GWAS, 'Disease']
            rtc_score = get_RTC_score(SNPS=snps[snps_in_LD,], EXPR=expr[id_gene,], id_GWAS=id_GWAS, id_eQTL=id_eQTL);
            result = rbind(result, c(id_eQTL, id_gene, FDR, id_GWAS, trait, rtc_score))
            message(paste("   ", id_GWAS, trait, id_gene, rtc_score))
        }
    }
}

colnames(result) = c("eQTL_SNP", "Gene", "FDR_eQTL", "GWAS_SNP","Complex_Trait", "RTC");

message("## Detected RTC results ...")

write.table(result, paste(eqtl_file_name, "RTC", "xls", sep="."), quote =F, sep="\t", row.names=F)

# Since the eQTL result might contain redudancey (e.g. if one SNP is reported as eQTL, the other SNPs in a strong LD are also reported), so we need to clean up in the output table
# for one GWAS-gene pair, only take the one with best RTC score
#require(dplyr)
#result <- summarise(group_by(), )


