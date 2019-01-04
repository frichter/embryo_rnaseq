# Felix Richter
# felix.richter@icahn.mssm.edu
# 9/23/2018
# description: run WGCNA for embryo data
##############################################################

# export TMP=/sc/orga/projects/chdiTrios/Felix/embryo_rnaseq/tmp_dir
# module load R/3.5.1 ##  R/3.4.1 ## 3.5.1 does not support colorout
# R

## set the home directory
setwd("D:/Dropbox/PhD/")
setwd("/Users/frichter/Dropbox (Personal)/PhD/")
setwd("/Users/felixrichter/Dropbox/PhD/")
setwd("/hpc/users/richtf01/")
options(stringsAsFactors=FALSE)

## load external libraries (order matters)
p = c("WGCNA", "DESeq2", "limma", "annotate", "org.Hs.eg.db", "topGO", "goseq", "gplots",
      ## generic dataprocessing packages:
      "readr", "magrittr", "purrr", "dplyr", "ggplot2", "tidyr")
lapply(p, require, character.only = TRUE)

####################################
# Analyze modules of interest 
# locally (following up on wgcna_embryo.R)
####################################

hla_neg_sup_path = c("HFE", "HLA-DOA", "HLA-DOB", "TAPBPL", "THBS1", "FCGR2B", "LILRB2")

# stage_i = 'Cleavage' ## Blastocyst Cleavage Compacted_Morula All_n72 Blastocyst_d7

# for All_n72:
mod_i = read_tsv('embryo_rnaseq/wgcna/All_n72/cytoscape_module/cs_edges_brown.txt')
# for cleavage: cs_edges_salmon.txt and purple
# mod_i = read_tsv('embryo_rnaseq/wgcna/Cleavage/cytoscape_modules/cs_edges_purple.txt')
mod_i %<>% select(fromNode:weight)
# mod_i %<>% filter(weight >= 0.03)

uniq_nodes = c(mod_i$fromNode, mod_i$toNode) %>% unique
n_interest = sum(hla_neg_sup_path %in% uniq_nodes)

mod_i %>% 
  # filter(weight >= 0.05) %>%
  # filter(fromNode %in% de_genes, toNode %in% de_genes) %>% 
  mutate(gene_interest_edges = 
           (fromNode %in% hla_neg_sup_path) | (toNode %in% hla_neg_sup_path)) %>% 
  mutate(n_bg = ifelse(gene_interest_edges, n_interest, length(uniq_nodes) - n_interest)) %>% 
  mutate(gene_interest_edges = ifelse(gene_interest_edges, "Y", "N") %>% factor) %>% 
  # group_by(gene_interest_edges) %>% 
  # summarise(med = median(weight), tot_degree = n(), deg_cent = tot_degree/unique(n_bg))
  summarise(wst_p = wilcox.test(weight ~ gene_interest_edges)$p.value)

ggplot(aes(x = weight, col = gene_interest_edges)) +
  geom_density() +
  scale_color_manual(values = c("black", "red")) +
  theme_classic()


Calc_degree_per_gene = function(gene_i, mod_i) {
  degree_i = mod_i %>% 
    # filter(weight >= 0.04) %>%
    filter((fromNode %in% gene_i) | (toNode %in% gene_i)) %>% 
    nrow
  return(degree_i)
}

node_degree = map(uniq_nodes, Calc_degree_per_gene, mod_i)

degree_df = cbind("gene" = uniq_nodes, "degree" = node_degree) %>% as.data.frame %>% 
  mutate(degree = as.numeric(degree)) %>% 
  mutate(gene_interest = gene %in% hla_neg_sup_path)

degree_df %>% filter(gene_interest)
mod_i %>% filter(fromNode == "HFE")

degree_df %>% arrange(desc(degree)) %>% mutate(deg_rank = rank(1/degree)) %>% 
  filter(gene_interest)
degree_df %>% mutate(gene_interest = factor(gene_interest)) %>% 
  summarise(p = wilcox.test(degree ~ gene_interest)$p.value)

##
gs_mm_loc = "embryo_rnaseq/wgcna/All_n72/bicor_signed_beta9_min30_minCoreKME5neg1_mergecutheight1neg1_static995_reassignThreshold1eneg6_minKMEtoStay1neg1_pamT_GS_MM_18_10_08.csv"
# top_n_genes 
gene_mm_gs = read_csv(gs_mm_loc)
gene_mm_gs %<>% select(-X1)

gene_mm_gs %>% 
  # filter(!(geneSymbol %in% hla_neg_sup_path)) %>% filter(moduleColor != "brown") %>% dim
  # group_by(moduleColor) %>% tally
  # filter(moduleColor == "brown") %>% 
  # select(geneSymbol, MM.brown, p.MM.brown) %>% 
  # filter(geneSymbol %in% de_genes) %>% 
  # arrange(p.MM.brown) %>% 
  # mutate(gene_rank = rank(p.MM.brown)) %>% 
  ## hla_neg_sup_path
  filter(geneSymbol %in% ipsc_path) %>%
  as.data.frame
# select(geneSymbol, moduleColor)

## over-representation in brown module
fisher.test(cbind(c(3, 4), c(1841-3, 26280)))
## beta=12 brown module (has )
fisher.test(cbind(c(3, 4), c(1231-3, 28125 - (1231 + 4))))

## over-representation in brown module + DE genes
fisher.test(cbind(c(3, 4), c(608-3, 28125 - (608 + 4))))

length(de_genes)
fisher.test(cbind(c(4, 3), c(2498-4, 28125 - (2498 + 3))))


### top connections per HLA gene
top_genes = 50

mod_i %>% filter(fromNode %in% hla_neg_sup_path, toNode %in% hla_neg_sup_path)

top_edges = mod_i %>% 
  # filter(weight >= 0.04) %>%
  mutate(gene_interest_edges = 
           (fromNode %in% hla_neg_sup_path) | (toNode %in% hla_neg_sup_path)) %>% 
  filter(gene_interest_edges) %>% 
  mutate(gene_i = ifelse(fromNode %in% hla_neg_sup_path, fromNode, toNode)) %>% 
  group_by(gene_i) %>% arrange(desc(weight)) %>% slice(1:top_genes) %>% ungroup %>% 
  select(fromNode:weight)

# check if any edges have a <0 correlation (they shouldn't)
sum((vsd[top_edges$fromNode, ] %>% t %>% cor) < 0)

# weight = |(1 + cor)/2|^beta
# 2*(w^(1/9)) - 1 = cor

top_edges %>% 
  ## back calculate correlation from the weight
  mutate(weight_cor =  abs(2*(weight^(1/9)) - 1)) %>% 
  ## consider adding in the direct connection between these two
  write_tsv('embryo_rnaseq/hla_expr/hla_brownMod_t50_edges_per_hlaSup.txt')

##### Plotting modules
infile = 'embryo_rnaseq/wgcna/All_n72/bicor_signed_beta12_min30_minCoreKME5neg1_mergecutheight1neg1_static995_reassignThreshold1eneg6_minKMEtoStay1neg1_pamT_trait_module_18_10_08.txt'
# infile = 'embryo_rnaseq/wgcna/Cleavage/bicor_signed_beta16_min30_mergecutheight1neg1_static995_minKMEtoStay1neg1_pamT_trait_module_18_10_08.txt'
mod_cor = read_tsv(infile)

cor_mtx = mod_cor %>% select(contains("cor")) %>% 
  ## for All_n72:
  select(Day.of.Developmentcor, Cleavagecor, Morulacor, Blastocystcor, Euploidcor) %>%
  ## for Cleavage:
  # select(Euploidcor, Day.of.Developmentcor, Number.of.Cellscor) %>%
  as.matrix
new_mod_names = c("I", "II", "III", "IV",
                  # "V", "VI", "VII", "VIII", "IX", "X", "XI", "XII", "XIII",
                  "Unassigned")
rownames(cor_mtx) = new_mod_names
mod_cor %>% select(module) %>% mutate(new_mod_names = new_mod_names)

colnames(cor_mtx) = colnames(cor_mtx) %>% gsub("cor$", "", .) %>% 
  gsub("Number.of.Cells", "Cell Count", .) %>% 
  gsub("Day.of.Development", "Day of Dev.", .)

color_scale = colorpanel(50, "Blue", "Black", "Yellow")
# color_scale = colorpanel(100, "Blue", "Black", "Yellow")[c(1:25, 45:55, 75:100)]
color_scale = c(rep(color_scale[[1]], 15), color_scale, rep(color_scale[[length(color_scale)]], 5))

filename = gsub("_18_10_08.txt", '_18_10_14.pdf', infile)
pdf(file = filename, width = 4.5, height = 3)
par(mar = c(6, 8.5, 3, 3))
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = cor_mtx,
               xLabels = colnames(cor_mtx),
               yLabels = rownames(cor_mtx) ,
               colorLabels = FALSE,
               colors = color_scale,
               setStdMargins = FALSE,
               cex.text = 0.2,
               # cex.lab.y = 0.5, 
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()
