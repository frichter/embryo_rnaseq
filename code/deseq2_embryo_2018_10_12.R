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

## load external libraries (order matters)
p = c("WGCNA", "DESeq2", "limma", "variancePartition", "gplots", "edgeR",
      "annotate", "org.Hs.eg.db", "topGO", "goseq",
      ## generic dataprocessing packages:
      "readr", "magrittr", "purrr", "dplyr", "ggplot2", "tidyr")
lapply(p, require, character.only = TRUE)


hla_neg_sup_path = c("HFE", "HLA-DOA", "HLA-DOB", "TAPBPL", "THBS1", "FCGR2B", "LILRB2")

################################################
# trying DESeq2 and per-gene plotting
################################################

## Blastocyst Cleavage Compacted_Morula All_n72 Blastocyst_d7 Blastocyst_icm_te
stage_i = 'All_n72'
# parent_dir = paste0("/sc/orga/projects/chdiTrios/Felix/embryo_rnaseq/expression_data/",
#                     stage_i, '/')
# locally:
parent_dir = paste0("embryo_rnaseq/expression_data/", stage_i, '/')

info = readRDS(paste0(parent_dir, 'info_2018_10_08.RDS'))
dds = readRDS(paste0(parent_dir, 'dds_2018_10_08.RDS'))


## trying DESeq2 on Day.of.Development Embryo_Stage + Batch + 
# design(dds) = ~ Day.of.Development + Arrested + Euploid
design(dds) = ~ Day.of.Development + Patient
dds = DESeq(dds)

# dds_deseq_pt_2018_12_20.RDS dds_deseq_noBatch_2018_10_10.RDS
saveRDS(dds, paste0(parent_dir, 'dds_deseq_pt_2018_12_20.RDS'))
# dds = readRDS(paste0(parent_dir, 'dds_deseq_noBatch_2018_10_10.RDS'))

resultsNames(dds)
## Arrested_ongoing_vs_arrested Day.of.Development Euploid_euploid_vs_aneuploid
res = results(dds, name = "Day.of.Development", independentFiltering = T)
summary(res)
# resOrdered = resLFC[order(resLFC$padj),]
resOrdered = res[order(res$padj),]
de_df = resOrdered %>%
  as.data.frame %>% 
  mutate(gene = row.names(resOrdered))
dayofdev_non_na_genes = de_df %>% filter(!is.na(padj)) %>% 
  select(gene) %>% unlist %>% as.character
  
de_df %>% filter(!is.na(padj)) %>% 
  dplyr::select(gene, everything()) %>%
  arrange(padj) %>% mutate(de_rank = rank(padj)) %>%
  # filter(log2FoldChange < 0) %>% 
  # filter(padj < 0.05) %>% select(gene) %>% unlist %>% as.character
  # filter(gene %in% hla_neg_sup_path) %>% as.data.frame
  filter(!is.na(padj)) %>%
  dplyr::select(gene, everything()) %>%
  # filter(padj < 0.05) %>% 
  select(gene, log2FoldChange, padj, baseMean) %>%
  # filter(log2FoldChange < 0) %>% 
  write_tsv('embryo_rnaseq/expression_data/All_n72/day_of_dev_pt_covar_all_genes.txt')
  
de_genes = de_df %>% filter(!is.na(padj)) %>% 
  dplyr::select(gene, everything()) %>%
  arrange(padj) %>% mutate(de_rank = rank(padj)) %>% 
  filter(log2FoldChange < -1) %>% ## 0 -2.5
  filter(padj < 0.05) %>% # arrange(log2FoldChange) %>% head
  select(gene) %>% unlist %>% as.character
de_genes_high = de_df %>% filter(!is.na(padj)) %>% 
  dplyr::select(gene, everything()) %>%
  arrange(padj) %>% mutate(de_rank = rank(padj)) %>%
  filter(log2FoldChange > 1) %>%
  filter(padj < 0.05) %>% select(gene) %>% unlist %>% as.character
  

#### GO analysis (using functions from GO_embryo.R)
expr_genes = row.names(resOrdered)
gene_map = getgo(expr_genes, 'hg19', 'geneSymbol')
sum(is.na(names(gene_map)))
sum(!is.na(names(gene_map)))
gene_map = gene_map[!is.na(names(gene_map))]

all_genes_w_scores = as.numeric(expr_genes %in% de_genes)
sum(all_genes_w_scores)
length(all_genes_w_scores)
names(all_genes_w_scores) = expr_genes
# ontology_i = "MF"
num_go_terms = 5

go_results = map_df(c("CC", "MF", "BP"), GetTopTermsPerOntology, all_genes_w_scores,
                    ## de_genes 
                    gene_map, de_genes, num_go_terms=5)

go_results %>% select(-sig_genes) %>% 
  mutate(term_rank = rank(classic)) %>% arrange(classic) %>% head
  # filter(grepl("negative", Term)) %>% 
  # filter(GO.ID == "GO:0002584")
  # filter(GO.ID == "GO:0002578") %>% as.data.frame
  filter(classic < 0.05) %>% 
  arrange(classic) %>% as.data.frame

go_results %>% filter(GO.ID == "GO:0003677") %>% select(sig_genes)

# group_by(ontology) %>% slice(1:5) %>% ungroup  

de_df %>% 
  # filter(gene %in% hla_neg_sup_path)
  # filter(gene %in% ipsc_path)
  # filter(gene %in% c('FGA', 'FGB'))
  filter(gene %in% c('KCNJ8', 'KCNJ1', 'KCNJ15', 'GNGT2'))
### plotting significant DE genes
vsd = readRDS(paste0(parent_dir, 'vsd_assay_2018_10_08.RDS'))
## 0002578
hla_neg_sup_path = c("HFE", "HLA-DOA", "HLA-DOB", "TAPBPL", "THBS1", "FCGR2B", "LILRB2")
ipsc_path = c("POU5F1", "SOX2", "KLF4", "MYC", "NANOG")
## 0002584
# hla_neg_sup_path = c("HFE", "HLA-DOA", "HLA-DOB", "TAPBPL")

gene_i = hla_neg_sup_path # ipsc_path # 
gene_i = gene_i[gene_i %in% rownames(vsd)]
## if only considering non-na genes:
# gene_i = gene_i[gene_i %in% dayofdev_non_na_genes]


vsd[gene_i, ] %>% t %>% cor
expr_long = vsd[gene_i, ] %>% as.data.frame %>% 
  mutate(gene = gene_i) %>% 
  gather(key = Specimen.ID, value = expr_norm, -gene) %>% 
  inner_join(info)

col_order = ifelse(sort(gene_i) %in% de_genes, 'red', 'black')
sort(gene_i) %in% de_genes_high

info %>% group_by(Embryo_Stage, Day.of.Development) %>% tally

max(expr_long$expr_norm)
p = expr_long %>% 
  # filter(grepl("HLA", gene)) %>% 
  # filter(gene %in% dayofdev_non_na_genes) %>%
  mutate(DEG = ifelse(gene %in% de_genes, 'Padj\n<0.05', 'NS')) %>% 
  # mutate(Day.of.Development = factor(Day.of.Development)) %>% 
  ## Number.of.Cells Day.of.Development Expansion
  ggplot(aes(x = Day.of.Development, y = expr_norm, col = gene)) +
  geom_point(col = "grey60", size = 0.25) +
  geom_smooth(method = "lm", se = T, fill = "grey60") +
  # geom_smooth(se = F) +
  facet_wrap(~gene, ncol = 2) + 
  scale_color_manual(values = col_order) +
  ylim(1, 15) +
  theme_classic()

p
ggsave('embryo_rnaseq/hla_expr/expr_v_t_GO0002578_facet.png', p, width = 3.25, height = 4.75)
  

p = expr_long %>% 
  mutate(DEG = ifelse(gene %in% de_genes, 'Padj\n<0.05', 'NS')) %>% 
  mutate(Day.of.Development = factor(Day.of.Development, levels = c(7, 6, 5, 3))) %>%
  ## Number.of.Cells Day.of.Development Expansion
  ggplot(aes(x = gene, y = expr_norm, fill = Day.of.Development)) +
  # geom_point(col = "grey60", size = 0.25) +
  # geom_smooth(method = "lm", se = F) +
  geom_boxplot(outlier.size = 0.25) +
  # geom_smooth(se = F) +
  scale_color_manual(values = rep("black", 4)) +
  scale_fill_manual(values = rep("grey", 4)) +
  ylim(1, 15) +
  coord_flip() +
  theme_classic()

p
 
