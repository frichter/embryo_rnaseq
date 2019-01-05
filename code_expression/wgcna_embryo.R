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

## installing custom packages in ~/rLocalPackages
# module load R/3.5.1
# R CMD INSTALL -l . colorout_1.1-1.tar.gz
# R CMD INSTALL -l . R_peer_source_1.3.tgz

library(peer,lib="/hpc/users/richtf01/rLocalPackages_v3_5_1/")
library(colorout,lib="/hpc/users/richtf01/rLocalPackages_v3_5_1/")

## load external libraries (order matters)
p = c("WGCNA", "DESeq2", "limma", "gplots",
      ## generic dataprocessing packages:
      "readr", "magrittr", "purrr", "dplyr", "ggplot2", "tidyr")
lapply(p, require, character.only = TRUE)

### WGCNA

###
source("d1_d2_rnaseq/code_pipeline/wgcna_plotting_functions.R")
enableWGCNAThreads(4)

## load data
stage_i = 'All_n72' ## Blastocyst Cleavage Compacted_Morula All Blastocyst_d7 All_n72
parent_dir = paste0("/sc/orga/projects/chdiTrios/Felix/embryo_rnaseq/expression_data/",
                    stage_i, '/')
# locally:
# parent_dir = paste0("embryo_rnaseq/expression_data/", stage_i, '/')

wgcna_dir = paste0(parent_dir, '../../wgcna/', stage_i, '/')
info = readRDS(paste0(parent_dir, 'info_2018_10_08.RDS'))
vsd = readRDS(paste0(parent_dir, 'vsd_assay_2018_10_08.RDS'))

# info = read_tsv(paste0(parent_dir, '../info_embryo.txt'))
# vsd = readRDS(paste0(parent_dir, "vsd_assay.RDS"))

datExpr = as.matrix(t(vsd))
colnames(datExpr) = rownames(vsd)
rownames(datExpr) = colnames(vsd)

## regressing out batch
# fit = lm(datExpr ~ as.factor(info$Batch))
# datExpr_res = residuals(fit)
# datExpr = datExpr_res

## create traits dataframe
if(stage_i == 'Cleavage') {
  datTraits = model.matrix( ~ Day.of.Development + Number.of.Cells + Euploid + Age + 0, info) %>% 
    as.data.frame %>% select(-Euploidaneuploid)
} else if(grepl('Blastocyst', stage_i)) {
  # ICM + TE + are NA, but only for 1 ID for d7 blastocysts
  info %<>% mutate_at(vars(ICM, TE), as.character)
  info %<>% mutate_at(vars(ICM, TE), function(x) ifelse(is.na(x), "0", x))
  datTraits = model.matrix( ~ Expansion + ICM + TE + # Day.of.Development + 
                              Arrested + Euploid + Age + 0, info) %>% as.data.frame
  datTraits %<>% select(-ICM0)
} else if(grepl('All_n72', stage_i)) {
  datTraits = model.matrix( ~ Embryo_Stage + Day.of.Development + 
                              Arrested + Euploid + Age + 0, info) %>% as.data.frame
}

head(datTraits)
# confirm same dimensions
dim(datTraits)
dim(vsd)

## previously used:
# datTraits = model.matrix( ~ Embryo_Stage + Age + Arrested + Euploid + 0, info) %>% as.data.frame
## clean column names:
names(datTraits) = names(datTraits) %>% 
  gsub("Embryo_Stage|euploid|Compacted_", "", .) %>% 
  gsub("Arrestedo", "O", .)

## none of the stages are highly correlated with Batch
# datTraits = model.matrix( ~ Embryo_Stage + Age + Arrested + Euploid + Batch + 0, info) %>% as.data.frame
# cor(datTraits)[1:7, ]

## soft thresholding
powers = c(c(1:10), seq(from = 12, to=30, by=2)) # , seq(from = 25, to=100, by=5)
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, networkType = "signed", verbose = 5)
filename = paste0(wgcna_dir, "soft_thresholds_signed_18_10_11.pdf")
PlotSoftThreshold(sft, filename) 

## checking out the HLA path
## form GO:0002584
hla_neg_sup_path = c("HFE", "HLA-DOA", "HLA-DOB", "TAPBPL")
## form GO:0002578
hla_neg_sup_path = c("HFE", "HLA-DOA", "HLA-DOB", "TAPBPL", "THBS1", "FCGR2B", "LILRB2") #[c(1:2, 5)]

# hla_neg_sup_path %in% colnames(datExpr)
# datExpr[, hla_neg_sup_path] %>% cor

## All_n72 beta 9 for >0.8, 12for >0.85, for 16 >0.9
## Blastocyst_d7 beta 10 for >0.8, 12 for >0.85, 14 for >0.9
beta_choice = 12
## 14 when using Blastocyst and Cleavage separately for 0.8, for 0.85 use 16
## no minimum threshold: 7
## min threshold only: plateaus at 8 with >0.8, hits 0.85 at 20..
## also tried: 9, 12
mod_size_choice = 30
wgcna_file_base = wgcna_dir %>% paste0(., "bicor_signed_beta", beta_choice, 
                                       "_min", mod_size_choice,
                                       "_minCoreKME5neg1_",
                                       "mergecutheight1neg1_static995_",
                                       "reassignThreshold1eneg6_",
                                       "minKMEtoStay1neg1_pamT_")
# max power for signed is 30, default beta for n<20 is 18 (FAQ 6)
# https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/faq.html
net = blockwiseModules(datExpr, power = beta_choice,
                       networkType = "signed",
                       TOMType = "signed", 
                       detectCutHeight = 0.995,
                       minModuleSize = mod_size_choice,
                       reassignThreshold = 1e-6,
                       minCoreKME = 0.5,
                       minKMEtoStay = 0.1, # 
                       mergeCutHeight = 0.1, # (1 - correlation btw eigengenes) for merging modules
                       corType="bicor",
                       numericLabels = TRUE,
                       pamStage = TRUE, 
                       # this also greatly determines how "clean". Geschwind uses negative (probably means false)
                       # pamRespectsDendro = TRUE,
                       maxBlockSize = 31000,
                       saveTOMs = TRUE,
                       saveTOMFileBase = wgcna_file_base,
                       verbose = 3)

# saveRDS(net, paste0(wgcna_file_base, "bwm_out_18_10_08.RDS"))
net = readRDS(paste0(wgcna_file_base, "bwm_out_18_10_08.RDS"))

moduleLabels = net$colors
moduleColors = labels2colors(moduleLabels)
modMembers = data.frame(Gene = colnames(datExpr), Module = moduleColors)

# modMembers %>% filter(Gene %in% hla_neg_sup_path)
unique(moduleColors) %>% length

# Plot the dendrogram and the module colors underneath
# pdf(file = paste0(wgcna_file_base, "dendro_genes_min20_18_04_23.pdf"), width = 12, height = 9)
pdf(file = paste0(wgcna_file_base, "dendro_genes_18_10_08.pdf"),
    width = 9, height = 4.5)
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

# get module eigengenes
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)

## label the gene expression level
# library(gplots)
color_scale = colorpanel(50, "Blue", "Black", "Yellow")
# displayColors(color_scale)

MET = orderMEs(cbind(MEs, datTraits))
## for no traits, use:
# MET = orderMEs(MEs)
filename = paste0(wgcna_file_base, "eigengene_dendro_18_10_08.pdf")
#  eigengene_dendro_18_08_30.pdf eigengene_dendro_NOTRAITS_18_08_30.pdf
# eigengene_dendro_MFonly_18_08_30.pdf eigengene_dendro_Method_18_08_30.pdf
pdf(file = filename, width = 6, height = 6)
# Plot the relationships among the eigengenes and the trait 
# Plot the dendrogram
par(cex = 1.0)
plotEigengeneNetworks(MET, "Eigengene dendrogram", marDendro = c(0,4,2,0),
                      signed = TRUE,
                      plotHeatmaps = FALSE)
# Plot the heatmap matrix (note: this plot will overwrite the dendrogram plot)
par(cex = 1.0)
plotEigengeneNetworks(MET, "Eigengene adjacency heatmap",
                      # marHeatmap = c(3,4,2,2),# previous
                      marHeatmap = c(7,7,2,2),
                      signed = TRUE,
                      heatmapColors = color_scale,
                      plotDendrograms = FALSE,
                      xLabelsAngle = 90)
dev.off()

# plot heatmap relationship between modules and traits
filename = paste0(wgcna_file_base, "trait_module_18_10_08.pdf")
moduleTraitCor = plotTraitModule(datExpr, moduleColors, datTraits, filename, color_scale)

# sending gene modules to files

trait_interest = "Number.of.Cells" # "Day.of.Development" # Expansion Number.of.Cells
filename = paste0(wgcna_file_base, "GS_MM_18_10_08.csv")
createGSMMTable(datExpr, moduleColors, datTraits, trait_interest, filename)


####################################
# Export every module to cytoscape
####################################
load(paste0(wgcna_file_base, "-block.1.RData"))

TOM = as.matrix(TOM)
# TOM = TOMsimilarityFromExpr(datExpr, power = beta_choice)

# Select modules
# modules = [[2]]
## remove vst_ from folder name
cytoprefix = gsub("bicor_signed_.*", "", wgcna_file_base)
cytoprefix = paste0(cytoprefix, "cytoscape_module/")
dir.create(cytoprefix)

## be sure to create a cytoscape_modules folder in the results folder
walk(unique(moduleColors), WrapperForCytoscapeExport, datExpr, TOM, moduleColors, cytoprefix) 

modules = "magenta"

## link for GO enrichment:
# https://labs.genetics.ucla.edu/horvath/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/FemaleLiver-04-Interfacing.R

