library("tidyverse")


### SNAKEMAKE I/O ###
compound_data        <- snakemake@input[['compound_data']]
cell_line_annotation  <- snakemake@input[['cell_line_annotation']]
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
    select(DRUG_ID,
           DRUG_NAME,
           AUC,
           LN_IC50,
           SANGER_MODEL_ID,
           lineage,
           PUTATIVE_TARGET,
           PATHWAY_NAME
           ) %>%
    left_join(y = lines_compounds, by = "DRUG_ID") %>%
    rename(drug_id = DRUG_ID,
           name = DRUG_NAME,
           auc = AUC,
           ic50 = LN_IC50,
           cell_id = SANGER_MODEL_ID,
           moa = PATHWAY_NAME,
           target = PUTATIVE_TARGET)

## Attempt to get DepMap IDs using SANGER passport
cell_line_annotation <- read_csv(file = "../../annotation/cell_line_annotation.csv") %>%
    select(Sanger_Model_ID, DepMap_ID) %>%
    rename(cell_id = Sanger_Model_ID,
           depmap_id = DepMap_ID) %>%
    filter(!is.na(cell_id))

filtered_annotated_data <- filtered_data %>%
    left_join(cell_line_annotation, by = "cell_id") %>%
    select(-cell_id)

write_csv(x = filtered_annotated_data, file = comma_file)
save(filtered_annotated_data, file = rdata)    