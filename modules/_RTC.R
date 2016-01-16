###########################################
## R script to calculate the RTC (Regulatory Trait Concordance) for eQTL SNPs
# Usage: Rscript $pipeline_path/modules/_RTC.R SNP.txt expression.txt cis_eQTL_output.txt
# Author: Xianjun Dong
# Version: 0.0
# Date: 2015-Apr-16
###########################################
args<-commandArgs(TRUE)

snps_file_name = args[1]  # in format of http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/Sample_Data/SNP.txt
expr_file_name = args[2]  # in format of http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/Sample_Data/GE.txt
eqtl_file_name = args[3]  # a tab-delimited file with haeder of "SNP     gene    beta    t-stat  p-value FDR"
ld_file_name   = args[4]  # a file to tell SNPs in LD interval, in format of "chr17	43714849	43714850 hg38_chr17_45637484_rs2942168 	1	.	Parkinson's_disease	rs113155081;rs62055661;rs558738552;|43752078;43752039;43751598;|0.997467;0.997467;0.99", where the 4th column has GWAS SNP id at the end.
gwas_file_name = args[5]  # a file of GWAS SNPs, in format of "chr	LDstart	LDend	LDID	gwas	SNPs"

# Note: it's possible to add covariate when calculating the residuals. See example code here:
# http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/CompareResults.html
cvrt_file_name = args[6] # in format of http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/Sample_Data/Covariates.txt

# test using Ganqiang's data
snps_file_name="/PHShome/gl871/neorogen/PDBPNeuroX/HBS_Feb21st2015/Matrix-eQTL/RTC/HC/HBS-HCSub.Matrix.uniq.genotype"
expr_file_name="/PHShome/gl871/neorogen/PDBPNeuroX/HBS_Feb21st2015/Matrix-eQTL/RTC/HC/HBS-HCSubTX.Expression.txt"
eqtl_file_name="/PHShome/gl871/neorogen/PDBPNeuroX/HBS_Feb21st2015/Matrix-eQTL/RTC/HC/HBS-HC.eQTL.trans.txt"
ld_file_name="/PHShome/gl871/neorogen/PDBPNeuroX/HBS_Feb21st2015/Matrix-eQTL/RTC/HC/HBS-HC.LD.GWAS.Genotype.txt"
gwas_file_name='/PHShome/gl871/neorogen/PDBPNeuroX/HBS_Feb21st2015/Matrix-eQTL/RTC/HC/HBS-HC.GWAS.Genotype.txt'

# Note: make sure all SNP id are dbSNP based. So, convert SNP id based on Illumina chip to dbSNP v144, for both eQTL and SNP genotype
# cat ~/eRNAseq/HCILB_SNDA/final.cis.eQTL.new.d1e6.p1e-2.xls | awk '$6<=0.05' > ~/eRNAseq/HCILB_SNDA/final.cis.eQTL.new.d1e6.p1e-2.FDRpt5.xls

# assuming the All.Matrix.SNP.ID and All.Matrix.txt are sorted already and have the same order
#awk '{OFS="\t"; if(NR>1) print $2,$3-1,$3,$1,$4}' ~/neurogen/genotyping_PDBrainMap/eQTLMatrixBatch123/All.Matrix.SNP.ID | LC_ALL=C sort --parallel=8 --buffer-size=5G -k1,1 -k2,2n | intersectBed -a - -b $GENOME/Annotation/Variation/dbSNP144.hg19.bed.groupped.SNP.unstranded -wo -sorted | cut -f1-5,9 > ~/neurogen/genotyping_PDBrainMap/eQTLMatrixBatch123/All.Matrix.SNP.ID.dbSNP144

snps_file_name="~/neurogen/genotyping_PDBrainMap/eQTLMatrixBatch123/All.Matrix.txt"
expr_file_name="~/eRNAseq/HCILB_SNDA/expression.postSVA.xls"
eqtl_file_name="~/eRNAseq/HCILB_SNDA/final.cis.eQTL.new.d1e6.p1e-2.FDRpt5.xls"
ld_file_name="/data/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Annotation/GWASCatalog/gwas_catalog_v1.0-downloaded.hg19.snps_in_LD.SNAP.LD_w250.r2_0.8.bed"
gwas_file_name="/data/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Annotation/GWASCatalog/gwas_catalog_v1.0-downloaded.hg19.bed"
snpid_file_name="/data/neurogen/genotyping_PDBrainMap/eQTLMatrixBatch123/All.Matrix.SNP.ID.dbSNP144"

## function
get_SNPs_in_LD <- function(tagSNP, LD)
{
    index=grep(paste0(tagSNP,";"),apply(LD[,c('LDID','genotype')],1,paste,collapse=";"))[1]  # only the first LD is taken, if the query SNP is found (with exact word match) in multiple LDs
    snpsinld = unlist(strsplit(as.character(LD[index,8]),"|", fixed=T))[1]
    gwasinld = unlist(strsplit(as.character(LD[index,4]),"_", fixed=T))[4]
    return(c(gwasinld, unlist(strsplit(as.character(snpsinld),";", fixed=T))))
}

get_GWAS_in_LD <- function(tagSNP, LD)
{
    index=grep(paste0("\\b",tagSNP,"\\b"),LD$genotype)[1]  # only the first LD is taken, if the query SNP is found (with exact word match) in multiple LDs
    return(unlist(strsplit(as.character(LD[index,4]),"_", fixed=T))[4])
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
        if(nrow(df)<5) next;
        res=resid(lm(y~x, df));
        # correlation test between continuous and categorical varialbles
        pvalue = anova(lm(res~df$x0))$`Pr(>F)`[1]
        #pvalues=c(pvalues, cor.test(res, df$x0, method="spearman")$p.value);  # Spearman Rank Correlation
        pvalues=c(pvalues, pvalue)
    }
    
    # Rank of GWAS SNP in all tested SNPs
    N=length(pvalues);
    Rank_GWAS = which(unique(sort(pvalues)) %in% pvalues[1])  # take the first hit in unique order if multiple SNPs have the same pvalue
    RTC = (N-Rank_GWAS)/N
    return(RTC);
}

message("## loading input files ...")
SNPID = read.table(snpid_file_name, header=F, stringsAsFactors =F);   # only those with dbSNP and IlluminaSNP ID; could be multiple dbSNPID for one IlluminaSNP
colnames(SNPID) = c('chr',	'start', 'end', 'SNPid_illumina', 'imputed', 'SNPid_dbSNP')

snps = read.table(snps_file_name, header=T, stringsAsFactors =F); rownames(snps) = snps[,1]; snps=snps[,-1];
expr = read.table(expr_file_name, header=T, stringsAsFactors =F); 
com=intersect(colnames(expr),colnames(snps))
snps=snps[,com]; expr=expr[,com]
#rownames(snps) = SNPID$SNPid_dbSNP[match(rownames(snps),SNPID$SNPid_illumina)]


LD   = read.delim(ld_file_name, header=F, stringsAsFactors =F);
colnames(LD) = c('chr',	'LDstart', 'LDend', 'LDID', 'score', 'strand', 'disease', 'SNPs_in_LD')
LD$genotype=do.call(rbind,strsplit(LD$SNPs_in_LD,"|", fixed=T))[,1]

GWAS = read.delim(gwas_file_name, header=F, stringsAsFactors =F);
colnames(GWAS) = c('chr',	'start', 'end', 'gwasID', 'score', 'strand', 'disease')
GWAS$id=do.call(rbind,strsplit(GWAS$gwasID,"_", fixed=T))[,4]

eqtl = read.table(eqtl_file_name, header=F, stringsAsFactors =F); 
colnames(eqtl) = c('SNPid_illumina', 'gene', 'beta', 't.stat', 'p.value', 'FDR', 'SNP.pos')
eqtl$SNPid_dbSNP = SNPID$SNPid_dbSNP[match(eqtl$SNPid_illumina,SNPID$SNPid_illumina)]

result = c();

message("## processing the eQTL gene-snp pairs ...")

for(i in 1:nrow(eqtl)){
    id_eQTL = eqtl[i,'SNPid_dbSNP'];
    if(is.na(id_eQTL)) next; # without corresponding dbSNP id
    
    id_eQTL = unlist(strsplit(id_eQTL,","))[1] # if multiple SNPs callapsed, only pick the first one
    id_gene = eqtl[i,'gene'];
    FDR = eqtl[i,'FDR'];
    
    message(paste0(" # record ",i,": ", id_eQTL," -- ", id_gene,"  ..."))
    
    # extreme case: eQTL itself is a GWAS 
    if(id_eQTL %in% GWAS$id){ 
        rtc_score = 1;
        id_GWAS = id_eQTL;
        trait = paste(GWAS[GWAS$id==id_GWAS, 'disease'], collapse = ';')
        
        result = rbind(result, c(id_eQTL, id_gene, FDR, id_GWAS, trait, rtc_score))
        message(paste("   --- ", id_GWAS, trait, id_gene, rtc_score))
        
    } else {
        snps_in_LD = get_SNPs_in_LD(tagSNP=id_eQTL, LD=LD)
        
        if(sum(is.na(snps_in_LD))>0) {
        	message("Oh no, I don't find any LD for this SNP!"); 
        	next; # skip if no LD found
        }
        
        index = snps_in_LD %in% GWAS$id;  # GWAS_in_LD = get_GWAS_in_LD(tagSNP=id_eQTL, LD=LD)
        if(sum(index)==0) {
        	message("No way! There is no GWAS SNP found in the LD!"); 
        	next; # skip if no GWAS in LD
        }  
        
        # convert dbSNP to illuminaID
        snps_in_LD_illumina = SNPID$SNPid_illumina[pmatch(snps_in_LD,SNPID$SNPid_dbSNP)]  #NOTE: can be NA if some snps in LD are not on Illumina ChiP
        
        # If >1 GWAS in the LD, take them one by one
        if(sum(index)>1) warning("More than one GWAS found in the LD, take them one by one for downstream analysis!")
        for(id_GWAS in snps_in_LD[index]){
            trait = GWAS[GWAS$id==id_GWAS, 'disease']
            SNPs=snps[snps_in_LD_illumina,]
            rownames(SNPs) = snps_in_LD
            
            rtc_score = get_RTC_score(SNPS=SNPs, EXPR=expr[id_gene,], id_GWAS=id_GWAS, id_eQTL=id_eQTL);
            result = rbind(result, c(id_eQTL, id_gene, FDR, id_GWAS, trait, rtc_score))
            message(paste("   === ", id_GWAS, trait, id_gene, rtc_score, nrow(result)))
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


