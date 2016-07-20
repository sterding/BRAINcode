require(xlsx) # install.packages('xlsx')
df=read.xlsx("~/Dropbox/PDBrainMap/manuscript/Nature/TableS7.eQTL.RTC.xls", sheetName = "top_eSNP.long", stringsAsFactors=F)
head(df)
cn=colnames(df)
df=cbind(df, do.call(rbind, strsplit(df$Associated.transcript.host.locus, "___")))
colnames(df)=c(cn, 'associatedTxhostSymbol','associatedTxhostEns','associatedTxhostType')

df[] <- lapply(df, as.character)  #http://stackoverflow.com/questions/2851015/convert-data-frame-columns-from-factors-to-characters

library(biomaRt); #source("http://bioconductor.org/biocLite.R"); biocLite('biomaRt')
library(dplyr)
options(dplyr.width = Inf)

myattributes <- c("ensembl_gene_id",
                  "entrezgene",
                  "hgnc_symbol",
                  "description",
                  # 'mim_morbid_description', # OMIM 
                  'namespace_1003',
                  'name_1006')  # GO term

fix_genes <- . %>% 
  #tbl_df %>% 
  distinct %>% 
  rename(ensgene=ensembl_gene_id,
         entrez=entrezgene,
         symbol=hgnc_symbol,
         #MIM=mim_morbid_description,
         GOdomain=namespace_1003,
         GO=name_1006) %>%
  filter(GOdomain!="") %>%  # only include no empty ones
  group_by(ensgene, entrez, symbol, description, GOdomain) %>%
  summarise(GO=paste(GO, collapse ="; "), n=n()) 

Ens=unique(sub("\\..*","",as.character(df$associatedTxhostEns)))
dd = useMart("ENSEMBL_MART_ENSEMBL",host="dec2013.archive.ensembl.org") %>% 
  useDataset(mart=., dataset="hsapiens_gene_ensembl") %>% 
  getBM(mart=., attributes=myattributes, filter='ensembl_gene_id', values=Ens) %>% 
  fix_genes

# convert GO to 3 columns: BP/CC/MF, according to the GOdomain
library(reshape2)
dd = dcast(dd, ensgene + entrez + symbol + description ~ GOdomain, value.var = 'GO', fill='NA')

# ## use org.Hs.egPATH or gageData
# #source("http://bioconductor.org/biocLite.R"); biocLite(c("BiocUpgrade","AnnotationDbi","org.Hs.eg.db"))
# library("AnnotationDbi")
# library("org.Hs.eg.db")
# x <- org.Hs.egPATH
# mapped_genes <- mappedkeys(x)
# xx <- as.list(x[mapped_genes])

# source("https://bioconductor.org/biocLite.R"); biocLite('gageData')
library(gageData)
data(kegg.sets.hs)
data(sigmet.idx.hs)
kegg.sets.hs = kegg.sets.hs[sigmet.idx.hs]
g=rep(seq_along(kegg.sets.hs), sapply(kegg.sets.hs, length))
index=g[match(dd$entrez, unlist(kegg.sets.hs))]
dd$KEGG=names(kegg.sets.hs)[index]

## add OMIM phynotype (download from OMIM http://omim.org/downloads/ijetP6U0QkuItLrgN39O3Q)
omim=read.delim("http://data.omim.org/downloads/ijetP6U0QkuItLrgN39O3Q/genemap2.txt", header = F, sep="\t", stringsAsFactors = F, comment.char = "#")
head(omim)
colnames(omim)=c("Chromosome", "Genomic_Position_Start", "Genomic_Position_End", "Cyto_Location", "Computed_Cyto_Location", "Mim_Number", "Gene_Symbols", "Gene_Name", "Approved_Symbol", "Entrez_Gene_ID", "Ensembl_Gene_ID", "Comments", "Phenotypes", "Mouse_Gene_Symbol")
dd$OMIM=omim[match(dd$entrez, omim$Entrez_Gene_ID),"Phenotypes"]

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

dim(dd)

df=cbind(df, dd[match(sub("\\..*","",as.character(df$associatedTxhostEns)), dd$ensgene),])

write.xlsx(df, file='TableS7.eQTL.RTC.xls', sheetName = "all_eSNP2", append=F)
