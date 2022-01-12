library("tidyverse")


### SNAKEMAKE I/O ###
compound_data        <- snakemake@input[['compound_data']]
compound_meta        <- snakemake@input[['compound_meta']]
lines_compounds      <- snakemake@input[['compound_lines']]

comma_file <- snakemake@output[['csv_db']]
rdata      <- snakemake@output[['rdata_db']]

full_table <- compound_data %>%
    map(read_csv) %>%
    reduce(bind_rows)

## Annotate how many lines were profiled by comp.
lines_compounds <- read_csv(lines_compounds) 

## Now keep columns of interest
# broad_id, name, auc, ic50, depmap_id, lineage, moa
filtered_data <- full_table %>%
    select(master_cpd_id,
           broad_cpd_id,
           area_under_curve,
           apparent_ec50_umol,
           DepMap_ID,
           lineage
           ) %>%
    left_join(y = lines_compounds, by = "broad_cpd_id") %>%
    rename(broad_id = broad_cpd_id,
           auc = area_under_curve,
           ec50 = apparent_ec50_umol,
           depmap_id = DepMap_ID)


## We need to annotate the moa for CRTP
compound_meta <- read_tsv(compound_meta) %>%
    select(master_cpd_id,
           gene_symbol_of_protein_target,
           target_or_activity_of_compound
           ) %>%
    rename(target = gene_symbol_of_protein_target,
           moa = target_or_activity_of_compound)

filtered_annotated_data <- filtered_data %>%
    left_join(compound_meta, by = "master_cpd_id") %>%
    select(-master_cpd_id)


write_csv(x = filtered_annotated_data, file = comma_file)
save(filtered_annotated_data, file = rdata)    