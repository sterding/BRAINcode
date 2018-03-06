## Rscript to summerize the eQTL results
## Output: $eqtl_file_name.group_by_gene.annotated.xls

args<-commandArgs(TRUE)

eqtl_file_name = args[1]
gene_or_HTNE   = args[2]

# # debug
# eqtl_file_name = '~/eRNAseq/HCILB_SNDA/final.cis.eQTL.d1e6.p1e-2.FDRpt5.xls'
# gene_or_HTNE   = 'TNE'

library(biomaRt); #source("http://bioconductor.org/biocLite.R"); biocLite('biomaRt')
library(dplyr)
library(reshape2)
options(dplyr.width = Inf)

eqtl = read.table(eqtl_file_name, header = F, stringsAsFactors=F)
colnames(eqtl) = c("SNP", "gene", "beta", "t.stat", "p.value", "FDR", "SNP.pos", "chr")

# convert Illumina ID to dbSNP ID [optional]
snpid_file_name="~/neurogen/genotyping_PDBrainMap/eQTLMatrixBatch123/All.Matrix.SNP.ID.dbSNP144" # only those with dbSNP and IlluminaSNP ID
SNPID = read.table(snpid_file_name, header=F, stringsAsFactors =F);   
colnames(SNPID) = c('chr',	'start', 'end', 'SNPid_illumina', 'imputed', 'SNPid_dbSNP')
SNPID$SNPid_dbSNP = sub("([^,]*),.*","\\1",SNPID$SNPid_dbSNP) # only the first if multiple dbSNP id mapped to one IlluminaSNP
eqtl$SNPid_dbSNP = SNPID$SNPid_dbSNP[match(eqtl$SNP,SNPID$SNPid_illumina)]
  
# add host gene of eSNP
annotation=read.table("~/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.bed", stringsAsFactors =F, header=F); 
colnames(annotation)=c('chr','gene.start','gene.end','gene.ID','score','strand','gene.symbol','gene.type')
eSNP.hostgene = inner_join(eqtl, annotation, by="chr") %>% 
  filter(SNP.pos >= gene.start & SNP.pos <= gene.end) %>% 
  group_by(SNP) %>%
  summarise(eSNP_host_gene=paste(unique(gene.symbol),collapse=";")) 
eqtl =  left_join(eqtl, eSNP.hostgene, by = 'SNP')

# group by eGene
eqtl = eqtl %>% group_by(gene) %>%
  summarize(chr=first(chr),
            n_SNPs = n(),
            prox_pos_in_Mb = paste(sort(unique(round(SNP.pos/1e6))),collapse="-"),
            list_of_SNP_pos=paste(SNP.pos,collapse=";"),
            list_of_SNP_ID=paste(SNPid_dbSNP,collapse=";"),
            list_of_SNP_hostgene=paste(unique(eSNP_host_gene),collapse=";"),
            minP=min(p.value), minFDR=min(FDR))
eqtl = mutate(eqtl, list_of_SNP_hostgene=sapply(strsplit(gsub("NA;|;NA","",list_of_SNP_hostgene),';'), function(x) paste(sort(unique(x)),collapse = ';')))

# add eGene type and name (host gene if for TNE)
if(grepl("eRNA|HTNE|TNE", gene_or_HTNE)){
  annotation=read.table("~/eRNAseq/HCILB_SNDA/eRNA.characterize.xls", stringsAsFactors =F)
  eqtl$assocaitedRNA_hostgene_symbol=do.call(rbind,strsplit(annotation$f19.Hostgene[match(eqtl$gene, rownames(annotation))],"___"))[,1]
  eqtl$assocaitedRNA_hostgene_EnsID=do.call(rbind,strsplit(annotation$f19.Hostgene[match(eqtl$gene, rownames(annotation))],"___"))[,2]
  eqtl$type='TNE'
}
if(grepl("gene", gene_or_HTNE)){
  annotation=read.table("~/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.bed", stringsAsFactors =F, header=F); 
  eqtl$assocaitedRNA_hostgene_symbol = annotation[match(eqtl$gene, annotation$V4),7]
  eqtl$assocaitedRNA_hostgene_EnsID = eqtl$gene
  eqtl$type = annotation[match(eqtl$gene, annotation$V4),8]
  eqtl$type = ifelse(eqtl$type=='protein_coding','mRNA','ncRNA') # [opitonal]
}

# get description, MIM, GO for the assocaitedRNA_hostgene
Ens=unique(sub("\\..*","",as.character(eqtl$assocaitedRNA_hostgene_EnsID)))

myattributes <- c("ensembl_gene_id",
                  "entrezgene",
                  "hgnc_symbol",
                  "description",
                  # 'mim_morbid_description', # OMIM 
                  'namespace_1003',
                  'name_1006')  # GO term
dd = useMart("ENSEMBL_MART_ENSEMBL",host="dec2013.archive.ensembl.org") %>% 
  useDataset(mart=., dataset="hsapiens_gene_ensembl") %>% 
  getBM(mart=., attributes=myattributes, filter='ensembl_gene_id', values=Ens) 
dd = dd %>% distinct %>% 
  rename(ensgene=ensembl_gene_id,
         entrez=entrezgene,
         symbol=hgnc_symbol,
         #MIM=mim_morbid_description,
         GOdomain=namespace_1003,
         GO=name_1006) %>%
  filter(GOdomain!="") %>%  # only include no empty ones
  group_by(ensgene, entrez, symbol, description, GOdomain) %>%
  summarise(GO=paste(GO, collapse ="; "), n=n()) %>%
  dcast(ensgene + entrez + symbol + description ~ GOdomain, value.var = 'GO', fill='NA') # convert GO to 3 columns: BP/CC/MF, according to the GOdomain
dd$description=gsub(" \\[.*","",dd$description)

## add OMIM phynotype (download from OMIM http://omim.org/downloads/ijetP6U0QkuItLrgN39O3Q)
omim=read.delim("http://data.omim.org/downloads/ijetP6U0QkuItLrgN39O3Q/genemap2.txt", header = F, sep="\t", stringsAsFactors = F, comment.char = "#")
colnames(omim)=c("Chromosome", "Genomic_Position_Start", "Genomic_Position_End", "Cyto_Location", "Computed_Cyto_Location", "Mim_Number", "Gene_Symbols", "Gene_Name", "Approved_Symbol", "Entrez_Gene_ID", "Ensembl_Gene_ID", "Comments", "Phenotypes", "Mouse_Gene_Symbol")
dd$OMIM=omim[match(dd$entrez, omim$Entrez_Gene_ID),"Phenotypes"]
dd$OMIM[dd$OMIM==""]=NA

# add to eqtl table
eqtl=cbind(eqtl, dd[match(sub("\\..*","",as.character(eqtl$assocaitedRNA_hostgene_EnsID)), dd$ensgene),c('description','OMIM')])

write.table(eqtl, file = paste0(eqtl_file_name, '.group_by_gene.annotated.xls'), sep="\t",row.names = F)

quit('no')

## trashbin


# source("https://bioconductor.org/biocLite.R"); biocLite('gageData')
library(gageData)
data(kegg.sets.hs)
data(sigmet.idx.hs)
kegg.sets.hs = kegg.sets.hs[sigmet.idx.hs]
g=rep(seq_along(kegg.sets.hs), sapply(kegg.sets.hs, length))
index=g[match(dd$entrez, unlist(kegg.sets.hs))]
dd$KEGG=names(kegg.sets.hs)[index]

## add GDA (gene-disease assocaition) from DisGeNet (http://www.disgenet.org/)
library(disgenet2r) # library(devtools); install_bitbucket("albags/disgenet2r")
gq <- disgenetGene(gene = dd$entrez, 
                   database = "ALL", 
                   score = c(">", 0.1)
)
gqlist=extract(gq)
do.call(rbind, strsplit(as.character(gqlist$c1.diseaseClassName),"; "))
write.xlsx(gqlist, file="TableS7.topgene.disgenet.xls")
# gqlist=Reduce(rbind, by(gqlist,gqlist$c2.name,head, n=3)) # only top 3 per group
dd$disgenet=gqlist[match(dd$entrez, gqlist$c2.geneId),"Phenotypes"]

df=cbind(df, dd[match(sub("\\..*","",as.character(df$associatedTxhostEns)), dd$ensgene),])

write.xlsx(df, file='TableS7.eQTL.RTC.xls', sheetName = "all_eSNP2", append=F)