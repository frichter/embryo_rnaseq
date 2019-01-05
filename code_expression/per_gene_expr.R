# Felix Richter
# felix.richter@icahn.mssm.edu
# 9/23/2018
# description: Plot PCAs
##############################################################

## set the home directory
setwd("D:/Dropbox/PhD/")
setwd("/Users/frichter/Dropbox (Personal)/PhD/")
setwd("/Users/felixrichter/Dropbox/PhD/")
setwd("/hpc/users/richtf01/")
options(stringsAsFactors=FALSE)

p = c("DESeq2", "limma", "gplots",
      ## generic dataprocessing packages:
      "readr", "magrittr", "purrr", "dplyr", "ggplot2", "tidyr")
lapply(p, require, character.only = TRUE)

################################################
# Plot expression of specific genes
################################################

stage_i = 'All_n72' ## Blastocyst Cleavage Compacted_Morula All_n72 Blastocyst_d7
# parent_dir = paste0("/sc/orga/projects/chdiTrios/Felix/embryo_rnaseq/expression_data/",
#                     stage_i, '/')
# locally:
parent_dir = paste0("embryo_rnaseq/expression_data/", stage_i, '/')

info = readRDS(paste0(parent_dir, 'info_2018_10_08.RDS'))
dds = readRDS(paste0(parent_dir, 'dds_2018_10_08.RDS'))
vsd = readRDS(paste0(parent_dir, 'vsd_assay_2018_10_08.RDS'))

## see how expression changes if you regress out batch
# datExpr = as.matrix(t(vsd))
# colnames(datExpr) = rownames(vsd)
# rownames(datExpr) = colnames(vsd)
# fit = lm(datExpr ~ as.factor(info$Batch))
# datExpr_res = residuals(fit)
# datExpr = datExpr_res
# vsd = datExpr %>% t


hla_neg_sup_path = c("HFE", "HLA-DOA", "HLA-DOB", "TAPBPL", "THBS1", "FCGR2B", "LILRB2")

gene_i = hla_neg_sup_path # c("POU5F1", "SOX2", "KLF4", "MYC", "NANOG") ## hla_neg_sup_path
gene_i = gene_i[gene_i %in% rownames(vsd)]

vsd[gene_i, ] %>% t %>% cor
expr_long = vsd[gene_i, ] %>% as.data.frame %>% 
  mutate(gene = gene_i) %>% 
  gather(key = Specimen.ID, value = expr_norm, -gene) %>% 
  inner_join(info)

p = expr_long %>% 
  # filter(grepl("HLA", gene)) %>% 
  # mutate(Day.of.Development = factor(Day.of.Development)) %>% 
  ## Number.of.Cells Day.of.Development
  ggplot(aes(x = Day.of.Development, y = expr_norm, col = gene)) +
  geom_point() +
  geom_smooth(method = "lm", se = F) +
  # geom_smooth(se = F) +
  # facet_wrap(~Arrested) +
  theme_classic()

p


fit = lm(expr_norm ~ Day.of.Development, expr_long %>% filter(gene == "TAPBPL"))
summary(fit)


