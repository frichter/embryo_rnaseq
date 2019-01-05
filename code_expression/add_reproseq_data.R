# Felix Richter
# felix.richter@icahn.mssm.edu
# 9/23/2018
# description: run DESeq2 for embryo data
##############################################################

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


################################################
# Add reproseq per-chromosome data to DDS
################################################

## prepared with prep_chr_dataframe.py
repro_chr_df = read_tsv('embryo_rnaseq/reproseq/reproseq_chr_df.txt')
repro_chr_df %<>% mutate(Specimen.ID = paste0('Sample_', Specimen.ID))

stage_i = 'All_n72' ## Blastocyst Cleavage Compacted_Morula All_n72 Blastocyst_d7
# parent_dir = paste0("/sc/orga/projects/chdiTrios/Felix/embryo_rnaseq/expression_data/",
#                     stage_i, '/')
# locally:
parent_dir = paste0("embryo_rnaseq/expression_data/", stage_i, '/')

info = readRDS(paste0(parent_dir, 'info_2018_10_08.RDS'))
dds = readRDS(paste0(parent_dir, 'dds_2018_10_08.RDS'))

rnaseq_sample_ids = row.names(colData(dds))
rna_info = colData(dds) %>% as.data.frame %>% mutate(Specimen.ID = rnaseq_sample_ids)

# confirm that none of the chromosomes are correlated with >0.6 (there is one pair with 0.5-0.6)
# should be 24 (for 24 chromosomes) if none are correlated
sum(abs(cor(repro_chr_df %>% select(-Specimen.ID) %>% as.matrix)) > 0.6)
# 74/106 reproseq samples have RNAseq, 74/81 rnaseq samples have reproseq data
sum(repro_chr_df$Specimen.ID %in% rnaseq_sample_ids)

repro_chr_df %>% head %>% as.data.frame

## update design matrix and reproseq column data
# dds_full = readRDS(paste0(parent_dir, "dds_full.RDS"))
## possibly repeat filtering with just the samples with reproseq data
design(dds)
# design(dds_full) = ~
colData(dds)
## confirm all the colData is in repro_rna_chr_info
# names(colData(dds))[!(names(colData(dds)) %in% names(repro_chr_df))]
chr_mtx = repro_chr_df %>% select(chr1:chrY) %>% as.matrix
row.names(chr_mtx) = repro_chr_df$Specimen.ID
chr_mtx = chr_mtx[rownames(colData(dds)), ]

colData(dds) = cbind(colData(dds), chr_mtx)

saveRDS(dds, paste0(parent_dir, 'dds_per_chrom_2018_10_08.RDS'))

################################################
# Perform per-chromosome DESeq
################################################

# filtering lowly expressed genes
keep = rowSums(fpkm(dds) >= 1) >= 1
sum(keep)
# keep = keep & low_enough
sum(keep)
dds = dds[keep,]

design(dds) 

###############################################################
# Calculate fit between genes as a function of chromosome
# count, see if there is a trend
###############################################################

# get gene info from some database
info_cols = c("GENEBIOTYPE", "SEQNAME", "SEQSTRAND", "GENESEQSTART", "GENESEQEND")
gene_info = AnnotationDbi::select(
  EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86,
  keys=rownames(dds), columns=info_cols, keytype="SYMBOL")

vsd = readRDS(paste0(parent_dir, 'vsd_assay_2018_10_08.RDS'))

colData(dds) %>% as.data.frame %>% group_by(chr17) %>% tally

data_cor = colData(dds) %>% as.data.frame %>% 
  select(Day.of.Development, Number.of.Cells, X..Losses:Total...Errors,
         Age, chr1:chrY) %>% cor
dim(data_cor)
sum(abs(data_cor) > 0.8)
data_cor[, 'chr17']
data_cor['chr6', c(1:6)]

PerChromFit = function(chr_i, vsd, gene_info) {
  print(chr_i)
  per_chrom_genes = gene_info %>% filter(SEQNAME == chr_i) %>% select(SYMBOL) %>% 
    unlist %>% as.character %>% unique 
  expr_long = vsd[per_chrom_genes, ] %>% as.data.frame %>% 
    mutate(gene = per_chrom_genes) %>% 
    gather(key = Specimen.ID, value = expr_norm, -gene) %>% 
    inner_join(repro_chr_df)
  form_i = as.formula(paste0('expr_norm ~ chr', chr_i))
  fit = lm(form_i, expr_long)
  chr_coef = summary(fit)$coefficients %>% as.data.frame %>% slice(2) %>% 
    rename(lm_p = `Pr(>|t|)`) %>% select(Estimate, lm_p) %>% 
    mutate(chr_i = chr_i)
  return(chr_coef)
}

chr_vec = c(as.character(c(1:22)), "X", "Y")
chr_fit_df = map_df(chr_vec, PerChromFit, vsd, gene_info)

chr_fit_df %>% 
  arrange(lm_p) %>% head(10)

p = expr_long %>% 
  ggplot(aes(x = chr17, y = expr_norm)) +
  geom_point() +
  geom_smooth(method = "lm", se = T, col = 'red') +
  # geom_smooth(se = F) +
  # facet_wrap(~Arrested) +
  theme_classic()

p

### dang actually can't mix expr_norm for different genes into a single model...
### since not normalized to gene length so I think longer genes will be weighted more highly


###############################################################
# Run per-chrom DESeq, start w chr22 since smallest and also
# significant overall association between chr count and expr
###############################################################

stage_i = 'Cleavage' ## Blastocyst Cleavage Compacted_Morula All_n72 Blastocyst_d7
# parent_dir = paste0("/sc/orga/projects/chdiTrios/Felix/embryo_rnaseq/expression_data/",
#                     stage_i, '/')
# locally:
parent_dir = paste0("embryo_rnaseq/expression_data/", stage_i, '/')

dds = readRDS(paste0(parent_dir, 'dds_per_chrom_2018_10_08.RDS'))

chr_vec = c(as.character(c(1:22)), "X", "Y")

full_chr_form = paste0('chr', chr_vec) %>% paste(., collapse = " + ") %>% 
  paste0('~ Number.of.Cells + ', .) %>% 
  as.formula
# design(dds) = ~ Day.of.Development + Number.of.Cells + Euploid + chr22
# design(dds) = full_chr_form
design(dds) = ~Number.of.Cells + Euploid
dds = DESeq(dds)

# saveRDS(dds, paste0(parent_dir, 'dds_all_chrom_deseq_2018_10_13.RDS'))
saveRDS(dds, paste0(parent_dir, 'dds_ploidy_deseq_2018_10_13.RDS'))

## which contrasts are we running
resultsNames(dds)
res = results(dds, name = "Euploid_euploid_vs_aneuploid", independentFiltering = T)
summary(res)

resOrdered = res[order(res$padj),]
de_df = resOrdered %>%
  as.data.frame %>% 
  mutate(gene = row.names(resOrdered))

PerChromDEGFEt = function(chr_i, de_df, gene_info) {
  per_chrom_genes = gene_info %>% filter(SEQNAME == chr_i) %>% select(SYMBOL) %>% 
    unlist %>% as.character %>% unique 
  chr_fet_df = de_df %>% 
    mutate(chr_i_gene = gene %in% per_chrom_genes) %>% 
    filter(!is.na(padj)) %>% 
    mutate(p_sig = padj < 0.05) %>% 
    mutate(pos_lfc = log2FoldChange > 0) %>% 
    mutate(p_sig_pos_lfc = p_sig & pos_lfc) %>% 
    ## pos_lfc p_sig_pos_lfc p_sig
    mutate(col_interst = pos_lfc) %>% 
    group_by(chr_i_gene, col_interst) %>% tally %>% 
    ungroup %>% 
    mutate(test_col = paste0(chr_i_gene, col_interst)) %>% 
    mutate(test_col =  test_col %>% gsub("ALSE|RUE", "", .)) %>% 
    select(test_col, n) %>% 
    spread(key = test_col, value = n) %>% 
    mutate(or = (TT/TF)/(FT/FF),
           fet_p = fisher.test(cbind(c(TT, TF), c(FT, FF)))$p.value) %>% 
    mutate(chr_i = chr_i)
  return(chr_fet_df)
}

chr_vec = c(as.character(c(1:22)), "X", "Y")
chr_de_df = map_df(chr_vec, PerChromDEGFEt, de_df, gene_info)

chr_de_df %>% as.data.frame %>% arrange(fet_p)
## are the most significant chromosomes also the most correlated? No
data_cor = colData(dds) %>% as.data.frame %>% 
  select(Day.of.Development, Number.of.Cells, X..Losses:Total...Errors,
         Age, chr1:chrY) %>% cor
data_cor[, 'chr22'] %>% sort
# fisher.test(cbind(c(239, 125), c(10365, 7602)))

colData(dds) %>% as.data.frame %>% 
  group_by(chr17, chr22) %>% tally


###############################################################
# DESeq ploidy
###############################################################

design(dds) = ~ Number.of.Cells + Euploid
dds = DESeq(dds)

# saveRDS(dds, paste0(parent_dir, 'dds_all_chrom_deseq_2018_10_13.RDS'))
saveRDS(dds, paste0(parent_dir, 'dds_ploidy_deseq_2018_10_13.RDS'))

res = results(dds, name = "Euploid_euploid_vs_aneuploid", independentFiltering = T)
summary(res)

resOrdered = res[order(res$padj),]
de_df = resOrdered %>%
  as.data.frame %>% 
  mutate(gene = row.names(resOrdered))

purp_nodes = read_tsv('embryo_rnaseq/wgcna/Cleavage/cytoscape_modules/cs_nodes_purple.txt') %>% 
  rename(gene = nodeName) %>% select(gene)
salmon_nodes = read_tsv('embryo_rnaseq/wgcna/Cleavage/cytoscape_modules/cs_nodes_salmon.txt') %>% 
  rename(gene = nodeName) %>% select(gene)

de_df %>% 
  filter(padj < 0.05, log2FoldChange < 0) %>% arrange(padj) %>% 
  mutate(rank_de = rank(padj)) %>% 
  # select(gene) %>% write_tsv(paste0(parent_dir, 'deseq_sig_up_aneuploid.txt'))
  # filter(gene %in% c('NPM2', 'WEE2'))
  filter(gene %in% salmon_nodes$gene)

vsd = readRDS(paste0(parent_dir, 'vsd_assay_2018_10_08.RDS'))
info = readRDS(paste0(parent_dir, 'info_2018_10_08.RDS'))

gene_i = c('AC023085.1', 'ARHGAP31-AS1')

expr_long = vsd[gene_i, ] %>% as.data.frame %>% 
  mutate(gene = gene_i) %>% 
  gather(key = Specimen.ID, value = expr_norm, -gene) %>% 
  inner_join(info)

p = expr_long %>% 
  # filter(grepl("HLA", gene)) %>% 
  # mutate(Day.of.Development = factor(Day.of.Development)) %>% 
  ## Number.of.Cells Day.of.Development
  ggplot(aes(x = Euploid, y = expr_norm, col = gene)) +
  geom_jitter() +
  # geom_smooth(method = "lm", se = F) +
  theme_classic()

p

colData(dds) %>% as.data.frame %>% filter(Euploid == 'euploid')
  group_by(Euploid) %>% tally
