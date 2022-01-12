library("tidyverse")

### SNAKEMAKE I/O ###
compound_data        <- snakemake@input[['compound_data']]
lines_compounds      <- snakemake@input[['compound_lines']]

comma_file <- snakemake@output[['csv_db']]
rdata      <- snakemake@output[['rdata_db']]

full_table <- compound_data %>%
    map(read_csv) %>%
    reduce(bind_rows)

## Annotate how many lines were profiled by comp.
lines_compounds <- read_csv(lines_compounds) 

## Now keep columns of interest
filtered_data <- full_table %>%
    select(broad_id, name, auc, ic50, depmap_id, lineage, moa, target) %>%
    left_join(y = lines_compounds, by = "broad_id")



write_csv(x = filtered_data, file = comma_file)
save(filtered_data, file = rdata)