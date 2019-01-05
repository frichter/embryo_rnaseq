
###############################################################
# Run 1-vs-all DESeq, see if any can identify RNAseq samples
# that are NOT mosaic
###############################################################

stage_i = 'All_n72' ## Blastocyst Cleavage Compacted_Morula All_n72 Blastocyst_d7
# parent_dir = paste0("/sc/orga/projects/chdiTrios/Felix/embryo_rnaseq/expression_data/",
#                     stage_i, '/')
# locally:
parent_dir = paste0("embryo_rnaseq/expression_data/", stage_i, '/')

dds = readRDS(paste0(parent_dir, 'dds_per_chrom_2018_10_08.RDS'))

# filtering lowly expressed genes
keep = rowSums(fpkm(dds) >= 1) >= 1
sum(keep)
# keep = keep & low_enough
sum(keep)
dds = dds[keep,]

colData(dds)$Patient = factor(colData(dds)$Patient)

id_list = colData(dds)$Patient %>% unique %>% as.character

colData(dds)$spec_id_aneuploid = colData(dds) %>% as.data.frame %>% 
  mutate(spec_id_aneuploid = ifelse(Euploid != 'euploid', Specimen.ID, 'Euploid')) %>% 
  select(spec_id_aneuploid) %>% unlist %>% as.character %>% factor
# id_i = id_list[[19]]

# per_id_form = paste0('~ Day.of.Development + ', id_i) %>% as.formula
# design(dds) = ~ Day.of.Development + Patient
design(dds) = ~ spec_id_aneuploid
dds = DESeq(dds)
# saveRDS(dds, paste0(parent_dir, 'per_id_deseq/dds_deseq_', id_i, '_2018_10_14.RDS'))
saveRDS(dds, paste0(parent_dir, 'dds_deseq_per_Aneuploid_spec_id_2018_10_14.RDS'))

## loop over patients (only aneuploid patients!)

GetSomyDEGfraction = function(id_i, dds, gene_info, somy) {
  print(id_i)
  ## get DE genes
  # resultsNames(dds)
  # result_i = paste0('Patient_', id_i, '_vs_75888')
  result_i = paste0('spec_id_aneuploid_', id_i, '_vs_Euploid')
  res = results(dds, name = result_i, independentFiltering = T)
  summary(res)
  
  resOrdered = res[order(res$padj),]
  de_df = resOrdered %>%
    as.data.frame %>% 
    mutate(gene = row.names(resOrdered))
  
  if(somy == 'Trisomy') {
    de_df %<>% mutate(sig_lfc = log2FoldChange > 0)
  } else if(somy == 'Monosomy') {
    de_df %<>% mutate(sig_lfc = log2FoldChange < 0)
  }
  
  ## are upregulated DE genes associated with the corresponding trisomies?
  # get PT trisomies
  chr_str = colData(dds) %>% as.data.frame %>% filter(spec_id_aneuploid == id_i) %>%
    select(Trisomy) %>% unlist %>% as.character %>% unique
  chr_list = chr_str %>% strsplit(',') %>% unlist
  per_chrom_genes = gene_info %>% filter(SEQNAME %in% chr_list) %>% select(SYMBOL) %>% 
    unlist %>% as.character %>% unique
  
  chr_fet_df = de_df %>% 
    mutate(chr_i_gene = gene %in% per_chrom_genes) %>% 
    filter(!is.na(padj)) %>% 
    mutate(p_sig = padj < 0.005) %>%
    mutate(p_sig_lfc = p_sig & sig_lfc) %>% 
    ## pos_lfc p_sig_pos_lfc p_sig
    mutate(col_interst = p_sig_lfc) %>% 
    group_by(chr_i_gene, col_interst) %>% tally %>% 
    ungroup %>% 
    mutate(test_col = paste0(chr_i_gene, col_interst)) %>% 
    mutate(test_col =  test_col %>% gsub("ALSE|RUE", "", .)) %>% 
    select(test_col, n) %>% 
    spread(key = test_col, value = n) 
  if(ncol(chr_fet_df) == 4) {
    chr_fet_df %<>% 
      mutate(or = (TT/TF)/(FT/FF),
             fet_p = fisher.test(cbind(c(TT, TF), c(FT, FF)))$p.value,
             frac_degs_on_aneuploid_str = paste0(as.character(TT), '/', as.character((TT + FT))),
             frac_degs_on_aneuploid = TT/(TT + FT))
  }
  chr_fet_df %<>% mutate(Specimen.ID = id_i, chr_str = chr_str)
  print(chr_fet_df)
  return(chr_fet_df)
}

id_list = colData(dds) %>% as.data.frame %>% filter(!is.na(Trisomy)) %>% 
  select(spec_id_aneuploid) %>% unlist %>% as.character %>% unique

tri_de_df_01 = map_df(id_list, GetSomyDEGfraction, dds, gene_info, 'Trisomy')
tri_de_df_01_control = map_df(id_list, GetSomyDEGfraction, dds, gene_info, 'Monosomy')

tri_de_df_01 %>% arrange(desc(frac_degs_on_aneuploid))
tri_de_df_01 %>% arrange(fet_p) %>% filter(fet_p < 0.05/47) %>% 
  head %>% as.data.frame
## are downregulated DEGs on the monosomy chromosomes?

id_list = colData(dds) %>% as.data.frame %>% filter(!is.na(Monosomy)) %>% 
  select(spec_id_aneuploid) %>% unlist %>% as.character %>% unique
mon_de_df = map_df(id_list, GetSomyDEGfraction, dds, gene_info, 'Monosomy')
mon_de_df_control = map_df(id_list, GetSomyDEGfraction, dds, gene_info, 'Trisomy')

mon_de_df %>% arrange(desc(frac_degs_on_aneuploid))
mon_de_df %>% arrange(fet_p) %>% filter(fet_p < 0.05/47) %>% 
  head %>% as.data.frame

### Plotting the results

tri_de_df_01 %>% select(Specimen.ID, chr_str, contains('frac'), fet_p)

p_df_tri = colData(dds) %>% as.data.frame %>% 
  select(Specimen.ID) %>% 
  left_join(tri_de_df_01) %>%
  # left_join(mon_de_df) %>%
  filter(!is.na(chr_str)) %>% 
  mutate(TT = ifelse(is.na(TT), 0, TT)) %>% 
  mutate(FT = ifelse(is.na(FT), 0, FT)) %>% 
  mutate(frac_degs_on_aneuploid = TT/(TT + FT)) %>% 
  mutate(frac_degs_on_aneuploid_str = paste0(as.character(TT), '/', as.character((TT + FT)))) %>% 
  arrange(frac_degs_on_aneuploid) %>% 
  filter(!is.nan(frac_degs_on_aneuploid))

p = p_df %>% 
  mutate(chr_str = factor(chr_str, levels = unique(p_df$chr_str))) %>% 
  mutate(p_sig = ifelse((fet_p < 0.05/ncol(p_df)) & (or > 1), "P<0.05", "NS")) %>% 
  mutate(p_sig = ifelse(is.na(fet_p), "NS", p_sig)) %>% 
  ggplot(aes(x = chr_str, y = frac_degs_on_aneuploid, label = frac_degs_on_aneuploid_str, col = p_sig)) +
  geom_segment(aes(x = chr_str, y = 0, 
                   xend = chr_str, yend = frac_degs_on_aneuploid), color = "grey80") +
  geom_point(size = 3) + xlab("") + ylab("") +
  geom_text(nudge_y = 0.15, color = "grey60") + 
  scale_color_manual(values = c("black", "red")) +
  coord_flip() +
  # ylim(0, 0.55) +
  theme_classic()
p

## monosomy_deg.png trisomy_deg.png
filename = 'embryo_rnaseq/figures/ploidy_2018_10_14/monosomy_deg.png'
ggsave(filename, p, width = 6.5, height = 6.5, units = "in")


sig_id_tri = p_df_tri %>% filter((fet_p < 0.05/ncol(p_df)) & (or > 1)) %>% 
  select(Specimen.ID) %>% unlist %>% as.character

p_df_mono %>% 
  filter((fet_p < 0.05/ncol(p_df)) & (or > 1)) %>% 
  filter(Specimen.ID %in% sig_id_tri)

