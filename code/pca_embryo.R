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
# Plot separate PCAs for cleavage and
# blastocyst data
################################################

stage_i = 'Cleavage' ## Blastocyst Cleavage Compacted_Morula All_n72 Blastocyst_d7
# parent_dir = paste0("/sc/orga/projects/chdiTrios/Felix/embryo_rnaseq/expression_data/",
#                     stage_i, '/')
# locally:
parent_dir = paste0("embryo_rnaseq/expression_data/", stage_i, '/')

info = readRDS(paste0(parent_dir, 'info_2018_10_08.RDS'))
dds = readRDS(paste0(parent_dir, 'dds_2018_10_08.RDS'))
vsd = readRDS(paste0(parent_dir, 'vsd_2018_10_08.RDS'))

info %>% select(Patient, Euploid) %>% unique %>% dim

# together: Patient Age Batch Embryo_Stage Arrested Day.of.Development
# for cleavage: Number.of.Cells, Day.of.Development
# for blastocyst look at: Arrested, TE, ICM, Expansion
# for all: Batch, X..Losses, X..Gains, Total...Errors, MAPD, Euploid, Age
# colData(vsd)$Expansion = factor(colData(vsd)$Expansion)
group_i = "Number.of.Cells"
# DESeq2 is obnoxious and makes you modify the plotPCA function to get PCs beyond 1 and 2
p_data = plotPCA(vsd, intgroup=group_i, returnData=TRUE)
percentVar = round(100 * attr(p_data, "percentVar"), digits = 1)

p = p_data %>%
  # mutate(group = gsub("TD.*SpoolFinal_", "", group)) %>% 
  # mutate(group = factor(group)) %>%
  ggplot(., aes(PC1, PC2, color=group)) +
  geom_point(size=2, show.legend = T) + ## 
  xlab(paste0("PC1: ",percentVar[1],"% variance")) + ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  # coord_fixed() +
  ### order is: nuclear, ribo, wc
  # scale_color_manual(values = c("royalblue", "red", "black"), name = group_i) +
  # scale_color_manual(name = group_i, values = c("gold", "blue")) + ## c("blue", "grey")) + ##
  scale_color_continuous(low="blue", high="gold", name = group_i, na.value = 'grey90') +
  # geom_text(aes(label = group), col = "black", show.legend = FALSE, check_overlap = F, hjust = "inward") +
  theme_classic()
# p

# filename = paste0(parent_dir, "../figures/pca_2018_10_08/minfpkm1_pc1_pc2", group_i, "_2018_09_28.png")
filename = paste0("embryo_rnaseq/figures/pca_2018_10_08/", stage_i,
                  "_minfpkm1_pc1_pc2_", group_i, "_2018_10_14.png")
ggsave(filename, p, width = 4, height = 2.5, units = "in")


sample_grep = p_data %>% filter(PC1 > 0, PC2 < -10) %>% select(name) %>% 
  unlist %>% paste(., collapse = "|")

info %>% #group_by(Day.of.Development) %>% tally
  # filter(Dx == "Unexplained/idiopathic") %>% tally
  filter(grepl(sample_grep, Specimen.ID)) %>% 
  # filter(grepl("93869_C4_THS_017_BxE7|91122_C1_THS_027_BxE5|91330_C6_THS_018_BxE3", Specimen.ID)) %>% 
  as.data.frame

info %>% filter(is.na(ICM)) %>% head %>% as.data.frame

# Expansion
exp_df = p_data %>% filter(!is.na(Number.of.Cells)) %>% select(PC1, Number.of.Cells)
cor.test(exp_df$PC1, exp_df$Number.of.Cells)

exp_df = p_data %>% filter(!is.na(Expansion)) %>% select(PC1, Expansion)
cor.test(exp_df$PC1, exp_df$Expansion)

