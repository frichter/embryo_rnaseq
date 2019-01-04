# Felix Richter
# felix.richter@icahn.mssm.edu
# 12/20/2018
# description: run Variance Partition and dream for embryos
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

## load packages
library(peer,lib="/hpc/users/richtf01/rLocalPackages_v3_5_1/")
library(colorout,lib="/hpc/users/richtf01/rLocalPackages_v3_5_1/")
## load external libraries (order matters)
p = c("WGCNA", "DESeq2", "limma", "variancePartition", "gplots",
      ## generic dataprocessing packages:
      "readr", "magrittr", "purrr", "dplyr", "ggplot2", "tidyr")
lapply(p, require, character.only = TRUE)

################################################
# load and clean data
################################################

## Blastocyst Cleavage Compacted_Morula All_n72 Blastocyst_d7 Blastocyst_icm_te
stage_i = 'All_n72'

parent_dir = paste0("/sc/orga/projects/chdiTrios/Felix/embryo_rnaseq/expression_data/",
                    stage_i, '/')
# locally:
# parent_dir = paste0("embryo_rnaseq/expression_data/", stage_i, '/')

info = readRDS(paste0(parent_dir, 'info_2018_10_08.RDS'))
dds = readRDS(paste0(parent_dir, 'dds_2018_10_08.RDS'))
expr = readRDS(paste0(parent_dir, 'vsd_assay_2018_10_08.RDS'))

## convert to factors
info %<>% mutate(Patient = factor(Patient),
                 # Age_binned = cut(Age, 5) %>% factor) %>% 
                 Age_binned = ifelse(Age >= 35, '>=35', '<35') %>% factor) %>% 
  mutate_at(vars(Patient, Batch, Embryo_Stage, Arrested, Euploid), as.factor)


# info %>% filter(is.na(Expansion) | is.na(ICM) | is.na(TE))
# info %>% group_by(Day.of.Development, Expansion, ICM, TE) %>% tally
# info %>% group_by(Day.of.Development, Number.of.Cells) %>% tally

################################################
# Compute Canonical Correlation Analysis (CCA)
# between all pairs of variables
################################################

# returns absolute correlation value
form = ~ Batch + Embryo_Stage + Arrested + Patient + Age + Euploid + 
  Day.of.Development + Dx + Age_binned

C = canCorPairs( form, info)
# Plot correlation matrix
plotCorrMatrix( C )

################################################
# variance partition
################################################

info$Specimen.ID == colnames(expr)
rownames(info) = info$Specimen.ID

## for cleavage only, use Number.of.Cells
## For blastocyst, use Expansion, ICM, and TE

## ICM, TE, 
info %>% group_by(Day.of.Development, Expansion, Age_binned, Euploid) %>% tally
info %>% group_by(Day.of.Development, Euploid, Arrested, Age_binned) %>% tally

## (1|Age_binned) (1|Patient) + (1|Batch) + (1|Embryo_Stage) + 
form = ~ (1|Patient) + Day.of.Development +
  # (1|Euploid) ## All_n72
  # (1|Arrested) + ## All_n72
  # (1|ICM) + (1|TE) ## Blastocyst_wICM_TE
  Expansion ## Blastocyst
  # Number.of.Cells ## Cleavage
  # (1|Age_binned) ## not Blastocyst_wICM_TE

res = fitVarPartModel(expr[1:4,], form, info )

# evaluating for collinearity
colinearityScore(res[[3]])

## actually run variance partition
varPart = fitExtractVarPartModel(expr, form, info)


## check form and incl form in name
form
saveRDS(varPart, paste0(parent_dir, 'vp_pt_dod_ploid_2018_12_20.RDS'))
# varpart_2018_12_20.RDS vp_pt_dod_exp_icm_te_2018_12_20.RDS
# vp_pt_dod_exp_2018_12_20.RDS vp_pt_dod_Ncells_2018_12_20.RDS
# vp_pt_dod_ploid_2018_12_20.RDS

# getting the following error:
# Error in { : task 14985 failed - "Downdated VtV is not positive definite"
## error is due to info/form variable correlations

################################################
# plot variance partition results
################################################

## gene sets of interest
hla_neg_sup_path = c("HFE", "HLA-DOA", "HLA-DOB", "TAPBPL", "THBS1", "FCGR2B", "LILRB2")
ipsc_path = c("POU5F1", "SOX2", "KLF4", "MYC", "NANOG")

vp = readRDS(paste0(parent_dir, 'varpart_2018_12_20.RDS'))

vp = sortCols( vp )
plotPercentBars( vp[hla_neg_sup_path,] )
plotPercentBars( vp[ipsc_path,] )
plotVarPart( vp )

## problem: no directionality. But can definitely use as a sanity check that
## variance is not explained by Patient
vp %>% 
  mutate(gene = row.names(vp)) %>% 
  # filter(gene %in% hla_neg_sup_path)
  # filter(Age_binned > 0.5)
  # filter(Arrested > 0.5)
  arrange(desc(Arrested)) %>% ## Day.of.Development
  head(10) %>% select(gene) %>% unlist %>% paste(collapse = ', ')
  # filter(gene %in% c('FGA', 'FGB'))
  # filter(gene %in% c('KCNJ8', 'KCNJ1', 'KCNJ15', 'GNGT2'))

################################################
# dream
################################################

info$Specimen.ID == colnames(expr)
rownames(info) = info$Specimen.ID

## use the subset of variables important for variance
## (1|Age_binned) (1|Patient) + (1|Batch) + (1|Embryo_Stage) + 
form = ~ (1|Patient) + 
  Euploid +
  # (1|Arrested) +
  # (1|Age_binned) +
  Day.of.Development
# (1|Embryo_Stage) + (1|Arrested) alone works

# The variable to be tested should be a fixed effect
# Get the contrast matrix for the hypothesis test
L = getContrast( expr, form, info, "Euploideuploid")
L

# Fit the dream model on each gene
# Apply the contrast matrix L for the hypothesis test  
# By default, uses the Satterthwaite approximation for the hypothesis test
fitmm = dream( expr, form, info, L)
fiteb = eBayes( fitmm )

saveRDS(fiteb, paste0(parent_dir, 'dream_pt_dod_ploid_2018_12_20.RDS'))



################################################
# analyze dream results
################################################

fiteb = readRDS(paste0(parent_dir, 'dream_pt_dod_2018_12_20.RDS'))

top_genes = topTable( fiteb, number = 4e5)
top_genes %<>% mutate(gene = row.names(top_genes)) 

top_genes %>% select(gene, everything()) %>% 
  # filter(logFC < -2, adj.P.Val < 0.05)
  # filter(gene %in% hla_neg_sup_path)
  # filter(gene %in% ipsc_path)
  # filter(gene %in% c('FGA', 'FGB'))
  filter(gene %in% c('KCNJ8', 'KCNJ1', 'KCNJ15', 'GNGT2'))

top_genes %>% 
  select(gene, logFC, adj.P.Val, AveExpr) %>% 
  # filter(log2FoldChange < 0) %>% 
  write_tsv('embryo_rnaseq/expression_data/All_n72/dream_pt_dod_all_genes.txt')



