# Felix Richter
# felix.richter@icahn.mssm.edu
# 9/23/2018
# description: run DESeq2 for embryo data
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
p = c("WGCNA", "DESeq2", "limma", "variancePartition", "gplots",
      ## generic dataprocessing packages:
      "readr", "magrittr", "purrr", "dplyr", "ggplot2", "tidyr")
lapply(p, require, character.only = TRUE)

################################################
# Split cleavage and blastocyst data, export
# VST results
################################################

## Blastocyst Cleavage Compacted_Morula All_n72 Blastocyst_d7 Blastocyst_icm_te
stage_i = 'Blastocyst_icm_te'
## morula has Day.of.Development but no ICM or TE so cannot be considered blastocyst
## no Expansion, Number.of.Cells so cannot be considered cleavage

info = read_tsv('embryo_rnaseq/metadata/info_n72_repro_rna.txt')
txlength_per_id = readRDS('embryo_rnaseq/expression_data/txlength_n72_repro_rna_2018_10_08.RDS')
repro_cts = readRDS('embryo_rnaseq/expression_data/cts_n72_repro_rna_2018_10_08.RDS')

## subset for the desired stage
if(stage_i %in% c('Blastocyst', 'Cleavage')) {
  info %<>% filter(Embryo_Stage == stage_i)
}
## filter for Day.of.Development being only 7
# info %<>% filter(Day.of.Development == 7)
# stage_i %<>% paste0("_d7")
info %<>% filter(Embryo_Stage == 'Blastocyst') %>% filter(!is.na(TE), !is.na(ICM))

repro_cts = repro_cts[, info$Specimen.ID]
txlength_per_id = txlength_per_id[, info$Specimen.ID]

## DESeq2 wants factors
info %<>% mutate_at(vars(Batch, ICM, TE, Arrested, Euploid, Embryo_Stage), as.factor)

dds = DESeqDataSetFromMatrix(countData = repro_cts,
                             colData = info,
                             design = ~Batch)
assays(dds)[["avgTxLength"]] = txlength_per_id

## update to appropriate design matrix
if(stage_i == 'Cleavage') {
  ## all are Arrested
  ## Day.of.Development is 3 or 5
  design(dds) = ~ Batch + Day.of.Development + Number.of.Cells + Euploid
} else if(stage_i == 'Blastocyst') {
  design(dds) = ~ Batch + Day.of.Development + Expansion + ICM + TE + Arrested + Euploid
} else if(stage_i == 'All_n72') {
  design(dds) = ~ Embryo_Stage + Batch + Day.of.Development + Expansion +
    Arrested + Euploid
  ## ICM + TE + Number.of.Cells + 
} else if(stage_i == 'Blastocyst_icm_te') {
  design(dds) = ~ Batch + Day.of.Development + Expansion + ICM + TE + Arrested + Euploid
}

# filtering lowly expressed genes
# keep = rowSums(fpkm(dds) > 1/10) >= 1
## form GO:0002578
hla_neg_sup_path = c("HFE", "HLA-DOA", "HLA-DOB", "TAPBPL", "THBS1", "FCGR2B", "LILRB2")
keep = rowSums(fpkm(dds) >= 1) >= 1
sum(keep)

hla_neg_sup_path %in% rownames(dds[keep, ])
# keep = keep & low_enough
sum(keep)
dds = dds[keep,]

## saving
saveRDS(dds, paste0('embryo_rnaseq/expression_data/', stage_i, '/dds_2018_10_08.RDS'))
saveRDS(info, paste0('embryo_rnaseq/expression_data/', stage_i, '/info_2018_10_08.RDS'))

## performing VST
vsd = varianceStabilizingTransformation(dds, blind=TRUE)
dim(assay(vsd))

saveRDS(assay(vsd), paste0('embryo_rnaseq/expression_data/', stage_i, '/vsd_assay_2018_10_08.RDS'))
saveRDS(vsd, paste0('embryo_rnaseq/expression_data/', stage_i, '/vsd_2018_10_08.RDS'))
