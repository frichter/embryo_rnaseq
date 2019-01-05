
# Felix Richter
# felix.richter@icahn.mssm.edu
# 9/23/2018
# description: GO enrichment of the embryo output
##############################################################

## set the home directory
setwd("D:/Dropbox/PhD/")
setwd("/Users/frichter/Dropbox (Personal)/PhD/")
setwd("/Users/felixrichter/Dropbox/PhD/")
setwd("/hpc/users/richtf01/")
options(stringsAsFactors=FALSE)

## load external libraries (order matters)
p = c("limma", "edgeR", "annotate", "org.Hs.eg.db", "DESeq2", "topGO", "goseq",
      ## generic dataprocessing packages:
      "readr", "magrittr", "purrr", "dplyr", "ggplot2", "tidyr")
lapply(p, require, character.only = TRUE)


# parent_dir = "/sc/orga/projects/chdiTrios/Felix/embryology/"
# parent_dir = "embryo_rnaseq/expression_data/"
# vsd = readRDS(paste0(parent_dir, "min1fpkm_1sample/vsd_assay.RDS"))


########################################################
# load module results (GS MM table)
########################################################

# gs_mm_loc = "embryo_rnaseq/wgcna/All_n72/bicor_signed_beta9_min30_minCoreKME5neg1_mergecutheight1neg1_static995_reassignThreshold1eneg6_minKMEtoStay1neg1_pamT_GS_MM_18_10_08.csv"
gs_mm_loc = "embryo_rnaseq/wgcna/Cleavage/bicor_signed_beta16_min30_minCoreKME5neg1_mergecutheight1neg1_static995_reassignThreshold1eneg6_minKMEtoStay1neg1_pamT_GS_MM_18_10_08.csv"
# top_n_genes 
gene_mm_gs = read_csv(gs_mm_loc)
gene_mm_gs %<>% select(-X1)
first_mod_col = 11 # 15 for all, 11 for cleavage
names(gene_mm_gs)[first_mod_col:ncol(gene_mm_gs)] = paste0(
  names(gene_mm_gs)[first_mod_col:ncol(gene_mm_gs)], "_")

############################
# Define the gene universe
# (i.e., genes w GO terms)
############################

# expr_genes = unique(row.names(vsd))
expr_genes = gene_mm_gs$geneSymbol %>% unique
# names(expr_genes) 
# get GO IDs from gene names
# supportedGenomes()
# supportedGeneIDs()
gene_map = goseq::getgo(expr_genes, 'hg19', 'geneSymbol')
sum(expr_genes %in% names(gene_map))
# look for NAs in the gene_map (how many genes are represented on GO?)
sum(is.na(names(gene_map)))
sum(!is.na(names(gene_map)))
gene_map = gene_map[!is.na(names(gene_map))]

########################################################
# Get interesting per-module genes
########################################################

module_color_i = 'salmon'
GetModuleGenes = function(module_color_i, gene_mm_gs, gene_map, top_n_genes, first_mod_col) {
  print(module_color_i)
  module_mm_gs = gene_mm_gs %>% 
    filter(moduleColor == module_color_i) %>% 
    select(geneSymbol, contains("GS"), contains(paste0("MM.", module_color_i, "_")))
  names(module_mm_gs) = c(names(module_mm_gs)[1:(first_mod_col-2)], "Module_membership", "MM_p")
  # top_genes = module_mm_gs %>% arrange(MM_p) %>%
  #   filter(geneSymbol %in% names(gene_map)) %>% 
  #   # slice(1:top_n_genes) %>%
  #   select(geneSymbol) %>% unlist %>% as.character
  # return(top_genes)
  ## if dataframes instead of lists are desired:
  top_gene_df = module_mm_gs %>% arrange(MM_p) %>%
    filter(geneSymbol %in% names(gene_map)) %>%
    slice(1:top_n_genes) %>%
    select(geneSymbol) %>% mutate(mod = module_color_i)
  return(top_gene_df)
}

mod_list = gene_mm_gs$moduleColor %>% unique
names(mod_list) = mod_list

## for lists
mod_genes = map(mod_list, GetModuleGenes, gene_mm_gs, gene_map, 50, first_mod_col)

## for dataframes
mod_gene_df = map_df(mod_list, GetModuleGenes, gene_mm_gs, gene_map, 10, first_mod_col)
write_tsv(mod_gene_df, 'embryo_rnaseq/wgcna/Cleavage/t10_genes_per_module.txt')
# 
# names(gene_map) %>% as.data.frame %>% 
#   write_tsv('embryo_rnaseq/wgcna/bg_genes_w_GO_terms.txt')


##################################
## GO enrichment
##################################

# function to return indices of genes
mySelGenes = function(score) {
  return (score != 0)
}

GetTermGenesPerResult = function(go_id, gene_set_for_test, module_go) {
  # go_id = results_final[3, 1]
  term_genes = genesInTerm(module_go, go_id) %>% unlist %>% as.character
  term_genes_in_gene_set = paste(term_genes[term_genes %in% gene_set_for_test], collapse = ",")
  return(term_genes_in_gene_set)
}

GetTopTermsPerOntology = function(ontology_i, all_genes_w_scores,
                                  gene_map, gene_set_for_test, num_go_terms) {
  print(ontology_i)
  module_go = new("topGOdata",
                  description = "GO terms associated with a WGCNA module",
                  ontology = ontology_i, # CC MF BP
                  allGenes = all_genes_w_scores, 
                  geneSel = mySelGenes,
                  annot = annFUN.gene2GO,
                  gene2GO = gene_map)
  total_terms = usedGO(module_go) %>% length
  test_statistic = new("classicCount", testStatistic = GOFisherTest, name = "FET")
  results_fet = getSigGroups(module_go, test_statistic)
  print(results_fet)
  ## keep all terms with total_terms, keep top x with num_go_terms
  results_final = GenTable(module_go, classic = results_fet, topNodes = total_terms)
  ## note that p-value is NOMINAL and up to user to correct
  results_final %<>% mutate(ontology = ontology_i)
  ## append a comma-sep list of genes in significant terms
  gene_set_term_gene_list = map(results_final$GO.ID, GetTermGenesPerResult,
                                gene_set_for_test, module_go) %>% unlist
  results_final %<>% mutate(sig_genes = gene_set_term_gene_list)
  return(results_final)
}

### running one-offs

## take the subset of interest
gene_set_for_test = mod_genes[['salmon']]

## annotate all genes as either being or not being in the gene set of interest
names(gene_map)
## either use all expr_genes or only those with GO terms
gene_universe = names(gene_map) # expr_genes
all_genes_w_scores = as.numeric(gene_universe %in% gene_set_for_test)
sum(all_genes_w_scores)
length(all_genes_w_scores)
names(all_genes_w_scores) = gene_universe

## keep all terms instead of significant terms so that you can correct ofr multiple hypotheses
go_results = map_df(c("CC", "MF", "BP"), GetTopTermsPerOntology, all_genes_w_scores,
                    gene_map, gene_set_for_test, num_go_terms=5)
go_results %>% select(-sig_genes) %>% 
  filter(classic < 0.01) %>% arrange(classic) %>% 
  group_by(ontology) %>% slice(1:5) %>% as.data.frame

## benchmark how p-values were calculated
fisher.test(cbind(c(20, 354-20), c(673-20, 30942 - (673 + 334))))
fisher.test(cbind(c(20, 334), c(369-20, 14891 - (369 + 334))))

##################################
# Loop over modules for GO
# enrichment
##################################

EnrichPerModule = function(gene_set_for_test, mod_name, gene_map, expr_genes) {
  print(mod_name)
  print(length(gene_set_for_test))
  ## annotate all genes as either being or not being in the gene set of interest
  all_genes_w_scores = as.numeric(expr_genes %in% gene_set_for_test)
  print(sum(all_genes_w_scores))
  print(length(all_genes_w_scores))
  names(all_genes_w_scores) = expr_genes
  go_results = map_df(c("CC", "MF", "BP"), GetTopTermsPerOntology,
                      all_genes_w_scores, gene_map, gene_set_for_test, num_go_terms=5)
  go_results %<>% mutate(module = mod_name)
  return(go_results)
}

go_df = map2_df(mod_genes, names(mod_genes), EnrichPerModule,
                gene_map, expr_genes)
go_df %>% write_tsv('embryo_rnaseq/wgcna/Cleavage/module_all_GO_beta16.txt')

## 4463 tests for MF, 15395 for BP, 1882 for CC. These are *probably* not independent
go_df %>% 
  mutate(p_value = gsub("<1e-30|< 1e-30", "1e-30", classic) %>% as.numeric) %>% 
  filter(module %in% c("salmon", "purple")) %>%
  # filter(grepl("immune", Term)) %>%
  # group_by(ontology) %>% 
  filter(module != "grey") %>%
  mutate(p_adj = 0.05/n()) %>% ungroup %>% 
  # filter(p_value < p_adj) %>% 
  # group_by(module, ontology) %>% tally
  # select(-sig_genes) %>% #arrange(classic) %>% head
  # filter(ontology == "BP") %>%
  group_by(module) %>% arrange(p_value) %>% slice(1:4) %>%
  # filter(p_value == min(p_value)) %>%
  ungroup %>% as.data.frame
  # filter(module == "brown") %>% 
  # filter(Significant < 10) %>% 
  # select(module, p_value, Term, ontology, GO.ID, Significant, Annotated, sig_genes) %>% 
  # write_tsv('embryo_rnaseq/wgcna/module_GOterms_beta8_t200_top_paths.txt')

go_df = read_tsv('embryo_rnaseq/wgcna/Cleavage/module_all_GO_beta16.txt')
go_df %>% head



