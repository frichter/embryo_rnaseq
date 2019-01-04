
## 72/74 reproseq samples with RNAseq passed QC :)
reproseq_data %>% filter(passed.QC. == 1) %>% filter(Specimen.ID %in% rnaseq_sample_ids)

## confirm other columns match up
# reproseq_data %>% 
#   inner_join(rna_info, by = c("Specimen.ID")) %>% 
#   # filter(Age.x != Age.y)
#   # group_by(Arrested.x, Arrested.y, Euploid.x, Euploid.y) %>% tally

## check column overlaps, if all reproseq data is available then 
## don't bother joining and cleaning
names(reproseq_data)[!(names(reproseq_data) %in% names(rna_info))]
# repro_rna_info = reproseq_data %>% select(-Age, -Arrested, -Euploid) %>%
#   inner_join(rna_info, by = c("Specimen.ID")) %>% 
#   filter(passed.QC. == 1)

## clean the data
# AMH, Partner.Age, Male.factor is NA for everything
# repro_rna_info %<>% select(-AMH, -Partner.Age, -Male.factor, -Embryo.Stage,
#                            - Blastocyst, -Compacted_Morula, -Cleavage,
#                            ## Patients == Patient
#                            -Patients)
repro_rna_info = rna_info 

## Number.of.Cells is available for every Cleavage stage embryo
## Expansion, ICM, and TE are only for blastocyst
repro_rna_info %>% head %>% as.data.frame
repro_rna_info %>% group_by(Day.of.Development, Embryo_Stage, 
                            Number.of.Cells, Expansion, ICM, TE) %>% tally
repro_rna_info %>% group_by(Embryo_Stage, Expansion, ICM, TE) %>% tally

repro_rna_info %>% 
  filter(!is.na(Expansion)) %>% 
  group_by(Expansion, Day.of.Development, Embryo_Stage) %>% tally

# join with per-chromosome counts
repro_rna_chr_info = repro_rna_info %>% inner_join(repro_chr_df)

repro_rna_chr_info %>% head %>% as.data.frame