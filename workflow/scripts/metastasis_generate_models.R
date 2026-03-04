library("tidyverse")


### SNAKEMAKE I/O ###
metastasis_penetrance        <- snakemake@input[["metastasis_data"]]
cell_lines_annotation  <- snakemake@input[["cell_lines_annotation"]]
count_matrix           <- snakemake@input[["count_matrix"]]

met_model_candidates   <- snakemake@output[["met_model_candidates"]]
met_lines_profiled      <- snakemake@output[["met_lines_profiled"]]

## load tables
metastasis_penetrance <- readr:::read_csv(metastasis_penetrance) %>%
                  as.data.frame()

cell_lines_annotation <- data.table::fread(cell_lines_annotation) %>%
                         as.data.frame()

ccle_counts <- readRDS(count_matrix)

metastasis_penetrance <- left_join(
    x = metastasis_penetrance,
    y = cell_lines_annotation,
    by = c("cell_line" = "CCLEName")
)  

## Get rid of cell lines and met data for which no RNASeq profiling was done
available_lines <- intersect(colnames(ccle_counts), metastasis_penetrance$ModelID)

selected_models <- metastasis_penetrance %>%
    filter(ModelID %in% available_lines) %>%
    select(
        met_model,
        ModelID,
        mean, 
        penetrance,
        cell_line,
        OncotreeLineage,
        original_lineage,
        PrimaryOrMetastasis
        ) %>%
    mutate(
        OncotreeLineage = as_factor(OncotreeLineage),
        met_model = as_factor(met_model)
    )

## Get the nº. lines/compound
lines_and_mets <- selected_models %>%
    group_by(met_model) %>%
    summarise(profiled_lines = n_distinct(ModelID)) %>%
    as.data.frame()

write.csv(lines_and_mets, file = met_lines_profiled, row.names=FALSE)


if(!dir.exists(file.path(met_model_candidates))){
    dir.create(met_model_candidates)}

by(selected_models,
   selected_models$met_model,
   FUN=function(i) write.csv(i, paste0(met_model_candidates, '/', i$met_model[1], ".csv"), , row.names=FALSE))


