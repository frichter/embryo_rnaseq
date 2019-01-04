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
# cd ~/rLocalPackages
# R CMD INSTALL -l . colorout_1.1-1.tar.gz
# R CMD INSTALL -l . R_peer_source_1.3.tgz

library(peer,lib="/hpc/users/richtf01/rLocalPackages_v3_5_1/")
library(colorout,lib="/hpc/users/richtf01/rLocalPackages_v3_5_1/")

## load external libraries (order matters)
p = c("WGCNA", "DESeq2", "limma", "variancePartition", "gplots",
      ## generic dataprocessing packages:
      "readr", "magrittr", "purrr", "dplyr", "ggplot2", "tidyr")
lapply(p, require, character.only = TRUE)

parent_dir = "/sc/orga/projects/chdiTrios/Felix/embryo_rnaseq/expression_data/"
# parent_dir = "embryo_rnaseq/expression_data/"

## load Mondale's data
# load(paste0(parent_dir, "DESeq2_LRT.RDATA"))
# ls()
# saveRDS(dds_full, paste0(parent_dir, "dds_full.RDS"))
dds_full = readRDS(paste0(parent_dir, "dds_full.RDS"))
info = read_tsv(paste0(parent_dir, '../metadata/info_n72_repro_rna.txt'))

# keep the expr where you have metadata
info %<>% filter(Specimen.ID %in% colnames(dds_full))
dds_full = dds_full[, info$Specimen.ID]

# colData(dds_full) %>% as.data.frame %>% 
#   mutate(sample_id = row.names(colData(dds_full))) %>% 
#   write_tsv(paste0(parent_dir, 'info_embryo.txt'))

## keep only 

## keep a subset of genes
# low_enough = rowMax(fpkm(dds_full)) < 1e5 
# sum(low_enough)
# keep = rowSums(fpm(dds_full) >= 1) >= 2
keep = rowSums(fpkm(dds_full) >= 1) >= 1
sum(keep)
# keep = keep & low_enough
sum(keep)
dds_full = dds_full[keep,]

# fpkm(dds_full)["DNMT1", ]
# row.names(dds_full[!low_enough, ])
#   rowMax(fpkm(dds_full[!low_enough, ]))

## VST
vsd = varianceStabilizingTransformation(dds_full, blind=TRUE)
dim(assay(vsd))

# saveRDS(vsd, paste0(parent_dir, "min1fpkm_1sample/vsd_full.RDS"))
vsd = readRDS(paste0(parent_dir, "min1fpkm_1sample/vsd_full.RDS"))
# saveRDS(assay(vsd), paste0(parent_dir, "min1fpkm_1sample/vsd_assay.RDS"))

## PCA (get a sense of the data)
colData(vsd) %>% as.data.frame %>% group_by(Embryo_Stage, Arrested) %>% tally
colData(vsd) %>% as.data.frame %>% group_by(Batch) %>% tally
colData(vsd) %>% as.data.frame %>% group_by(Patient) %>% tally ## 81 samples from 23 patients
colData(vsd) %>% head

# Patient Age Batch Embryo_Stage Arrested
group_i = "Embryo_Stage"
# DESeq2 is obnoxious and makes you modify the plotPCA function to get PCs beyond 1 and 2
p_data = plotPCA(vsd, intgroup=group_i, returnData=TRUE)
percentVar = round(100 * attr(p_data, "percentVar"), digits = 1)
p = p_data %>%
  ggplot(., aes(PC1, PC2, color=group)) +
  geom_point(size=1.5) + ## 
  xlab(paste0("PC1: ",percentVar[1],"% variance")) + ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  # coord_fixed() +
  ### order is: nuclear, ribo, wc
  scale_color_manual(values = c("royalblue", "red", "black")) +
  # scale_color_manual(values = c("gold", "blue")) + ##c("blue", "grey")) + ##
  # geom_text(aes(label = group), col = "black", show.legend = FALSE, check_overlap = F, hjust = "inward") +
  theme_classic()
p

# filename = paste0("d1_d2_rnaseq/figures/pca_2018_08_28/", data_subset, 
#                   "_from_", data_source, "_pc3_pc4_", group_i, "_2018_08_28.pdf")
# filename = paste0(parent_dir, "../figures/pca_2018_09_28/minfpkm1_pc1_pc2", group_i, "_2018_09_28.png")
# ggsave(filename, p, width = 4, height = 2.5, units = "in")


################################################
# Combine metadata 10/8/2018
################################################

## load count data
parent_dir = "/sc/orga/projects/chdiTrios/Felix/embryology/"
parent_dir = "embryo_rnaseq/expression_data/"
dds_full = readRDS(paste0(parent_dir, "dds_full.RDS"))

# combine and clean metadata
rnaseq_sample_ids = row.names(colData(dds_full))
rna_info = colData(dds_full) %>% as.data.frame %>% mutate(Specimen.ID = rnaseq_sample_ids)
reproseq_data = read_tsv('embryo_rnaseq/reproseq/reproseq_data.txt')
names(reproseq_data) = names(reproseq_data) %>% make.names
reproseq_data %<>% mutate(Specimen.ID = paste0('Sample_', Specimen.ID))

# 74/106 reproseq samples have RNAseq, 74/81 rnaseq samples have reproseq data
sum(reproseq_data$Specimen.ID %in% rnaseq_sample_ids)
# rnaseq_sample_ids[!(rnaseq_sample_ids %in% reproseq_data$Specimen.ID)] %>% paste(collapse = ", ")
## 72/74 reproseq samples with RNAseq passed QC :)
reproseq_data %>% filter(passed.QC. == 1) %>% filter(Specimen.ID %in% rnaseq_sample_ids)

repro_rna_info = reproseq_data %>% select(-Age, -Arrested, -Euploid) %>%
  inner_join(rna_info, by = c("Specimen.ID")) %>% 
  filter(passed.QC. == 1)

## clean the data
# AMH, Partner.Age, Male.factor is NA for everything
## Number.of.Cells is available for every Cleavage stage embryo
## Expansion, ICM, and TE are only for blastocyst

repro_rna_info %<>% select(-AMH, -Partner.Age, -Male.factor, -Embryo.Stage,
                           - Blastocyst, -Compacted_Morula, -Cleavage,
                           ## removing Patients bc it is same as Patient
                           -Patients)

write_tsv(repro_rna_info, 'embryo_rnaseq/metadata/info_n72_repro_rna.txt')

## export counts and gene length data 
repro_cts = counts(dds_full)[, repro_rna_info$Specimen.ID]

saveRDS(repro_cts, 'embryo_rnaseq/expression_data/cts_n72_repro_rna_2018_10_08.RDS')

txlength_per_id = assays(dds_full)$avgTxLength[, repro_rna_info$Specimen.ID]
saveRDS(txlength_per_id, 'embryo_rnaseq/expression_data/txlength_n72_repro_rna_2018_10_08.RDS')

