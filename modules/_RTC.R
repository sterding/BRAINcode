###########################################
## R script to calculate the RTC (Regulatory Trait Concordance) for eQTL SNPs
# Usage: 
#    Rscript $pipeline_path/modules/_RTC.R ~/neurogen/genotyping_PDBrainMap/eQTLMatrixBatch123/All.Matrix.txt expression.postSVA.xls final.cis.eQTL.d1e6.p1e-2.FDRpt5.xls HTNE
#    Rscript $pipeline_path/modules/_RTC.R ~/neurogen/genotyping_PDBrainMap/eQTLMatrixBatch123/All.Matrix.txt expression.postSVA.xls final.cis.eQTL.d1e6.p1e-2.FDRpt5.xls gene
# bsub -q big-multi -n 4 -M 10000 Rscript $pipeline_path/modules/_RTC.R SNP.txt expression.txt cis_eQTL_output.txt
# Author: Xianjun Dong
# Version: 1.0
# Date: 2017-Mar-03
# Reference: http://biorxiv.org/content/biorxiv/suppl/2016/09/11/074682.DC1/074682-1.pdf
# Note: there are specifal cases including
# 1. GWAS SNP is same as eQTL SNP, so its RTC score is also 1
# 2. GWAS SNP is not eQTL SNP, but they have same genotypes, so its RTC score is also 1. In that case, we report both
###########################################
args<-commandArgs(TRUE)

snps_file_name = args[1]  # in format of http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/Sample_Data/SNP.txt
expr_file_name = args[2]  # in format of http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/Sample_Data/GE.txt
eqtl_file_name = args[3]  # a tab-delimited file with columns of "SNP     gene    beta    t-stat  p-value FDR", but no header line
gene_or_HTNE   = args[4]  # <gene|HTNE>

ld_file_name="~/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Annotation/GWASCatalog/gwas_catalog_v1.0-downloaded.hg19.snps_in_LD.SNAP.LD_w250.r2_0.8.bed"
# a file to tell SNPs in LD interval, in format of "chr17	43714849	43714850 hg38_chr17_45637484_rs2942168 	1	.	Parkinson's_disease	rs113155081;rs62055661;rs558738552;|43752078;43752039;43751598;|0.997467;0.997467;0.99", where the 4th column has GWAS SNP id at the end.

# # HTNEs
snps_file_name="~/neurogen/genotyping_PDBrainMap/eQTLMatrixBatch123/All.Matrix.txt"
expr_file_name="~/eRNAseq/HCILB_SNDA/expression.postSVA.xls"
eqtl_file_name="~/eRNAseq/HCILB_SNDA/final.cis.eQTL.d1e6.p1e-2.FDRpt5.xls"
# 
# # Genes
# snps_file_name="~/neurogen/genotyping_PDBrainMap/eQTLMatrixBatch123/All.Matrix.txt"
# expr_file_name="~/neurogen/rnaseq_PD/results/eQTL/HCILB_SNDA/expression.postSVA.xls"
# eqtl_file_name="~/neurogen/rnaseq_PD/results/eQTL/HCILB_SNDA/final.cis.eQTL.d1e6.p1e-2.FDRpt5.xls"

# Note: make sure all SNP id are dbSNP based. So, convert SNP id based on Illumina chip to dbSNP v144, for both eQTL and SNP genotype
# cat ~/eRNAseq/HCILB_SNDA/final.cis.eQTL.d1e6.p1e-2.xls | awk '$6<=0.05' > ~/eRNAseq/HCILB_SNDA/final.cis.eQTL.d1e6.p1e-2.FDRpt5.xls

# assuming the All.Matrix.SNP.ID and All.Matrix.txt are sorted already and have the same order
#awk '{OFS="\t"; if(NR>1) print $2,$3-1,$3,$1,$4}' ~/neurogen/genotyping_PDBrainMap/eQTLMatrixBatch123/All.Matrix.SNP.ID | LC_ALL=C sort --parallel=8 --buffer-size=5G -k1,1 -k2,2n | intersectBed -a - -b $GENOME/Annotation/Variation/dbSNP144.hg19.bed.groupped.SNP.unstranded -wo -sorted | cut -f1-5,9 > ~/neurogen/genotyping_PDBrainMap/eQTLMatrixBatch123/All.Matrix.SNP.ID.dbSNP144
snpid_file_name="~/neurogen/genotyping_PDBrainMap/eQTLMatrixBatch123/All.Matrix.SNP.ID.dbSNP144"

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
    
    # get pvalues for all SNPs in SNPS
    pvalues=c()
    gwas_pvalue = 0;
    
    for(i in 1:nrow(SNPS)){
        df=data.frame(y=as.numeric(EXPR), x=as.numeric(SNPS[i,]), x0=as.numeric(eQTL_SNP))
        df=na.exclude(df);
        if(nrow(df)<5) next;
        res=resid(lm(y~x, df));
        # Using ANOVA to test correlation between continuous and categorical varialbles
        pvalue = anova(lm(res~df$x0))$`Pr(>F)`[1]
        #pvalues=c(pvalues, cor.test(res, df$x0, method="spearman")$p.value);  # Spearman Rank Correlation
        pvalues=c(pvalues, pvalue)
        
        if(rownames(SNPS)[i] == id_GWAS) gwas_pvalue = pvalue
    }
    
    # Rank of GWAS SNP in all tested SNPs
    N=length(pvalues);
    
    ## XD: code before 2017/03/03, where the sorting is in wrong order
    #Rank_GWAS = which(unique(sort(pvalues)) %in% pvalues[1])  # take the first hit in unique order if multiple SNPs have the same pvalue
    ## XD: bug fixed on 2017/03/03
    Rank_GWAS = which(sort(pvalues, decreasing =T) %in% gwas_pvalue)[1] - 1  # take the first order if multiple SNPs have the same pvalue
    
    RTC = (N-Rank_GWAS)/N
    return(RTC);
}

message("## loading input files ...")
SNPID = read.table(snpid_file_name, header=F, stringsAsFactors =F);   # only those with dbSNP and IlluminaSNP ID; could be multiple dbSNPID for one IlluminaSNP
colnames(SNPID) = c('chr',	'start', 'end', 'SNPid_illumina', 'imputed', 'SNPid_dbSNP')
SNPID$SNPid_dbSNP = sub("([^,]*),.*","\\1",SNPID$SNPid_dbSNP) # only the first if multiple dbSNP id for one IlluminaSNP

snps = read.table(snps_file_name, header=T, stringsAsFactors =F); rownames(snps) = snps[,1]; snps=snps[,-1];
expr = read.table(expr_file_name, header=T, stringsAsFactors =F);
com=intersect(colnames(expr),colnames(snps))
snps=snps[,com]; expr=expr[,com]

LD = read.delim(ld_file_name, header=F, stringsAsFactors =F);
colnames(LD) = c('chr',	'LDstart', 'LDend', 'LDID', 'score', 'strand', 'disease', 'SNPs_in_LD')
LD$genotype=do.call(rbind,strsplit(LD$SNPs_in_LD,"|", fixed=T))[,1]

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
    pvalue = eqtl[i,'p.value'];
    
    message(paste0(" # record ",i,": ", id_eQTL," -- ", id_gene,"  ..."))
    
    index=grep(paste0(id_eQTL,";"),apply(LD[,c('LDID','genotype')],1,paste,collapse=";"))
    if(sum(index)==0) {
    	#message("It's not from any LD!"); 
    	next; # skip if not in LD
    }
    
    for(ld in index){
    	snpsinld = unlist(strsplit(as.character(LD[ld,8]),"|", fixed=T))[1]
    	gwasinld = unlist(strsplit(as.character(LD[ld,4]),"_", fixed=T))[4]
    	trait = LD$disease[ld]
    	
    	if(id_eQTL == gwasinld){  # eQTL SNP itself is a GWAS 
    		rtc_score = 1;
		 		result = rbind(result, c(id_eQTL, id_gene, pvalue, FDR, id_eQTL, trait, rtc_score))
    		message(paste("   --- ", gwasinld, trait, id_gene, rtc_score, nrow(result)))
    	} else {  # eQTL is not a GWAS 
    		# get all distinct SNPs in LD, incl. the GWAS SNPs
    		snps_in_LD = unique(c(gwasinld, unlist(strsplit(as.character(snpsinld),";", fixed=T))));
    		# convert dbSNP to illuminaID
    		snps_in_LD_illumina = SNPID$SNPid_illumina[match(snps_in_LD,SNPID$SNPid_dbSNP)]
    		# if the GWAS SNP is not genotyped
    		if(is.na(snps_in_LD_illumina[1])) next;
    		# remove other NA (if any snps in LD are not on Illumina ChiP)
    		snps_in_LD = snps_in_LD[!is.na(snps_in_LD_illumina)]
    		snps_in_LD_illumina = snps_in_LD_illumina[!is.na(snps_in_LD_illumina)]
    		
    		SNPs=snps[snps_in_LD_illumina,]; rownames(SNPs) = snps_in_LD
    	
    		rtc_score = get_RTC_score(SNPS=SNPs, EXPR=expr[id_gene,], id_GWAS=gwasinld, id_eQTL=id_eQTL);
    		result = rbind(result, c(id_eQTL, id_gene, pvalue, FDR, gwasinld, trait, rtc_score))
    		message(paste("   === ", gwasinld, trait, id_gene, rtc_score, nrow(result)))
    	}
    }
}

colnames(result) = c("eQTL_SNP", "assocaitedRNA", "pvalue_eQTL", "FDR_eQTL", "GWAS_SNP","Complex_Trait", "RTC");
result=as.data.frame(result)

message("## Detected RTC results ...")

# Since the eQTL result might contain redudancey (e.g. if one SNP is reported as eQTL, the other SNPs in a strong LD are usually also eQTL), so we need to prune the output table, e.g. for each GWAS-gene pair by only taking the one with best RTC score 
#result=read.table(paste(eqtl_file_name, "RTC", "xls", sep="."), sep="\t", header=T)

if(grepl("eRNA|HTNE", gene_or_HTNE)){
  annotation=read.table("~/eRNAseq/HCILB_SNDA/eRNA.characterize.xls", stringsAsFactors =F)
  result = cbind(result, assocaitedRNA_hostgene=do.call(rbind,strsplit(annotation$f19.Hostgene[match(result$assocaitedRNA, rownames(annotation))],"___"))[,1])
  result$assocaitedRNA_coord=result$assocaitedRNA
}
if(grepl("gene", gene_or_HTNE)){
  annotation=read.table("~/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.bed", stringsAsFactors =F, header=F); 
  result = cbind(result, assocaitedRNA_hostgene=annotation[match(result$assocaitedRNA, annotation$V4),7])
  # replace gene coordinate with gene EnsID
  result$assocaitedRNA_coord=apply(annotation[match(result$assocaitedRNA, annotation$V4),1:3],1,paste0,collapse="_")
}
# add postion of SNPs
result=cbind(result, pos_GWAS_SNP=SNPID$end[match(result$GWAS_SNP, SNPID$SNPid_dbSNP)])
result=cbind(result, pos_eQTL_SNP=SNPID$end[match(result$eQTL_SNP, SNPID$SNPid_dbSNP)])

write.table(result, paste(eqtl_file_name, "RTC", "xls", sep="."), quote =F, sep="\t", na='NA', row.names=F)

## annotate: add hostgene_GWAS_SNP and hostgene_eQTL_SNP
# sed 's/ /___/g' final.cis.eQTL.d1e6.p1e-2.FDRpt5.xls.RTC.xls | awk '{OFS="\t"; split($2,a,"_"); if(NR>1) print a[1],$9-1,$9,$0}'  | intersectBed -a - -b <(cat /data/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.bed | awk '{OFS="\t"; print $1,$2,$3,$7,$5,$6}') -wao | cut -f4-13,17 | groupBy -g 1-10 -c 11 -o distinct | awk '{OFS="\t"; split($2,a,"_"); print a[1],$10-1,$10,$0}' | intersectBed -a - -b <(cat /data/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.bed | awk '{OFS="\t"; print $1,$2,$3,$7,$5,$6}') -wao | cut -f4-14,18 | groupBy -g 1-11 -c 12 -o distinct | sed 's/___/ /g' > final.cis.eQTL.d1e6.p1e-2.FDRpt5.xls.RTC.annotated.xls

## filter: for each GWAS SNP, only take the one with best RTC score (if there are multiple eSNPs in LD)
# sed 's/ /___/g' final.cis.eQTL.d1e6.p1e-2.FDRpt5.xls.RTC.xls | sort -k2,2r -k5,5 -k7,7gr | awk '{if(id!=$5) print; id=$5;}' | awk '{OFS="\t"; split($2,a,"_"); if(NR>1) print a[1],$9-1,$9,$0}'  | intersectBed -a - -b <(cat /data/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.bed | awk '{OFS="\t"; print $1,$2,$3,$7,$5,$6}') -wao | cut -f4-13,17 | groupBy -g 1-10 -c 11 -o distinct | awk '{OFS="\t"; split($2,a,"_"); print a[1],$10-1,$10,$0}' | intersectBed -a - -b <(cat /data/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.bed | awk '{OFS="\t"; print $1,$2,$3,$7,$5,$6}') -wao | cut -f4-14,18 | groupBy -g 1-11 -c 12 -o distinct | sed 's/___/ /g' > final.cis.eQTL.d1e6.p1e-2.FDRpt5.xls.RTC.filtered.annotated.xls