library("tidyverse")

### SNAKEMAKE I/O ###
ctrp_db        <- snakemake@input[['ctrp_db']]
gdsc_db        <- snakemake@input[['gdsc_db']]
prism_db      <- snakemake@input[['prism_db']]
cell_line_annotation <- snakemake@input[['cell_line_annotation']]

db_save    <- snakemake@output[['full_database']]
db_save_rd <- snakemake@output[['full_database_rd']]

db_summary_save    <- snakemake@output[['drug_summary']]
db_summary_save_rd <- snakemake@output[['drug_summary_rd']]


prism <- read_csv(prism_db) %>%
    select(-cid) %>%
    mutate(project = "PRISM") %>%
    rename(drug_id = broad_id,
           cid = smiles
           )

ctrp  <- read_csv(ctrp_db) %>%
    mutate(project = "CTRP") %>%
    rename(drug_id = broad_id,
           cid = smile
           )

gdsc  <- read_csv(gdsc_db) %>%
    mutate(project = "GDSC") %>%
    rename(cid = pubchem) %>%
    mutate(drug_id = as.character(drug_id))

## Cell line annotation
cell_line_annotation <- read_csv(cell_line_annotation) %>%
    select(ModelID,
           StrippedCellLineName,
           original_lineage,
           OncotreeLineage,
           OncotreeSubtype,
           PrimaryOrMetastasis
           ) %>%
    rename(depmap_id = ModelID,
           cell_line = StrippedCellLineName,
           origin = original_lineage,
           origin_subtype = OncotreeSubtype,
           tumor_type = PrimaryOrMetastasis
           )

## attempt to merge the data by outer join and add project 
combined_data <- prism %>%
    bind_rows(ctrp) %>%
    bind_rows(gdsc) %>%
    select(-OncotreeLineage) %>%
    left_join(cell_line_annotation, "depmap_id") %>%
    rename(celligner_origin = OncotreeLineage)


## generate summary
drug_summary <- combined_data %>%
    group_by(project, drug_id) %>%
    summarise(name = first(name),
              target = first(target),
              moa = first(moa),
              lines_tested = first(profiled_lines),
              min_auc = min(auc), 
              max_auc = max(auc),
              mean_auc = mean(auc),
              cid = first(cid)
              ) %>%
    distinct()


write_csv(combined_data, file = db_save)
write_csv(drug_summary, file = db_summary_save)

save(combined_data, file = db_save_rd)    
save(drug_summary, file = db_summary_save_rd)    
