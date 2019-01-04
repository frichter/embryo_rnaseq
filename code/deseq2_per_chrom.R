# Felix Richter
# felix.richter@icahn.mssm.edu
# 10/5/2018
# description: DESeq per chromosome
##############################################################

# export TMP=/sc/orga/projects/chdiTrios/Felix/embryology/tmp_dir
# module load R/3.5.1 ##  R/3.4.1 ## 3.5.1 does not support colorout
# R

## set the home directory
setwd("D:/Dropbox/PhD/")
setwd("/Users/frichter/Dropbox (Personal)/PhD/")
setwd("/Users/felixrichter/Dropbox/PhD/")
setwd("/hpc/users/richtf01/")
options(stringsAsFactors=FALSE)

## installing custom packages in ~/rLocalPackages
# module load R/3.5.1
# R CMD INSTALL -l . colorout_1.1-1.tar.gz
# R CMD INSTALL -l . R_peer_source_1.3.tgz

library(peer,lib="/hpc/users/richtf01/rLocalPackages_v3_5_1/")
library(colorout,lib="/hpc/users/richtf01/rLocalPackages_v3_5_1/")

## load external libraries (order matters)
p = c("WGCNA", "DESeq2", "limma", "variancePartition",
      ## generic dataprocessing packages:
      "readr", "magrittr", "purrr", "dplyr", "ggplot2", "tidyr")
lapply(p, require, character.only = TRUE)

parent_dir = "/sc/orga/projects/chdiTrios/Felix/embryology/"
# parent_dir = "embryo_rnaseq/expression_data/"

dds = readRDS(paste0(parent_dir, 'dds_reproseq.RDS'))


# filtering lowly expressed genes
keep = rowSums(fpkm(dds) >= 1) >= 1
sum(keep)
# keep = keep & low_enough
sum(keep)
dds = dds[keep,]

# design(dds) = ~Batch + Embryo_Stage + Arrested + Euploid + chr1

# get gene info from some database
# AnnotationDbi::columns(org.Hs.eg.db::org.Hs.eg.db)
# ensembldb::listColumns(EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86)
# ensembldb::transcripts(EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86)
# ensembldb::listColumns(org.Hs.eg.db::org.Hs.eg.db)
# AnnotationDbi::columns(EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86)
info_cols = c("GENEBIOTYPE", "SEQNAME", "SEQSTRAND", "GENESEQSTART", "GENESEQEND")

gene_info = AnnotationDbi::select(
  EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86,
  keys=rownames(dds), columns=info_cols, keytype="SYMBOL")

per_chrom_genes = gene_info %>% filter(SEQNAME == '22') %>% select(SYMBOL) %>% unlist %>% as.character %>% 
  unique

dds = dds[per_chrom_genes, ]
design(dds) = ~Batch + Embryo_Stage + Arrested + chr22
dds = DESeq(dds)

## which contrasts are we running
resultsNames(dds)
res = results(dds, name = "chr22", independentFiltering = T)
summary(res)
## print RPKMs to file (only for genes with non-NA p-values)
res %>% as.data.frame %>% 
  mutate(gene = row.names(res)) %>%
  filter(!is.na(padj)) %>% 
  arrange(padj) %>% head

dds
vsd = varianceStabilizingTransformation(dds, blind=TRUE)

assay(vsd)
expr_df = assay(vsd)[c('GSTT2' , 'EMID1'), ] %>% t
expr_df %<>% as.data.frame %>% 
  mutate(Specimen.ID = row.names(expr_df))

## Batch + Embryo_Stage + Arrested + 
colData(dds) %>% 
  as.data.frame %>% 
  mutate(Specimen.ID = row.names(colData(dds))) %>% 
  inner_join(expr_df) %>% 
  ggplot(aes(x = chr22, y = GSTT2, col = Embryo_Stage)) + 
  geom_point() +
  geom_smooth(method = 'lm', se = F) + 
  scale_color_manual(values = c(''))
  theme_classic()







