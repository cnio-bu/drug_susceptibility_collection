
library("tidyverse")

### SNAKEMAKE I/O ###
gdsc_dose_response_curves      <- snakemake@input[["dose_response_curves"]]
cell_lines_annotation          <- snakemake@input[["cell_lines_annotation"]]
count_matrix                   <- snakemake@input[["count_matrix"]]

auc_models_candidates          <- snakemake@output[["auc_models_candidates"]]
compounds_lines_profiled       <- snakemake@output[["compounds_lines_profiled"]]

# Setup log file
log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")


## Load tables

candidate_curves <- readxl::read_xlsx(gdsc_dose_response_curves) %>%
                    filter(RMSE <= 0.3) %>%
                    as.data.frame()

cell_lines_annotation <- data.table::fread(cell_lines_annotation) %>%
                        as.data.frame()

ccle_counts <- readRDS(count_matrix)

## Filter out drug/lines associations where the line is not profiled
available_lines    <- colnames(ccle_counts)
rm(ccle_counts)

## common lines check
common_lines <- intersect(
    candidate_curves$SANGER_MODEL_ID, 
    $Sanger_Model_ID
    )

## Filter the curves based on expression availability and annotation
candidate_curves_filtered <- candidate_curves %>%
    filter(
        SANGER_MODEL_ID %in% common_lines
    ) %>%
    left_join(
        y = cell_lines_annotation[, c("Sanger_Model_ID", "DepMap_ID", "lineage")],
        by = c("SANGER_MODEL_ID" = "Sanger_Model_ID")
    ) %>%
    filter(
        DepMap_ID %in% available_lines
    )

## get the number of profiled lines/compound
lines_by_compound <- candidate_curves_filtered %>%
    group_by(DRUG_ID) %>%
    summarise(
        profiled_lines = n_distinct(SANGER_MODEL_ID),
        cv  = sd(AUC) / mean(AUC)
    ) %>%
    as.data.frame()

write.csv(lines_by_compound, file=compounds_lines_profiled, row.names=FALSE)

## Get rid of models with less than 10 profiled lines
compounds_to_test <- lines_by_compound %>%
    filter(profiled_lines >= 10 & cv >= 0.1) %>%
    arrange(desc(cv)) %>% 
    distinct(DRUG_ID, .keep_all = TRUE) %>%
    pull(DRUG_ID)

candidate_curves   <- filter(candidate_curves, DRUG_ID %in% compounds_to_test)

candidate_curves$lineage <- as.factor(candidate_curves$lineage)

if(!dir.exists(file.path(auc_models_candidates))){
    dir.create(auc_models_candidates)}

by(candidate_curves,
   candidate_curves$DRUG_ID,
   FUN=function(i) write.csv(i, 
                             paste0(auc_models_candidates, 
                                    "/", i$DRUG_ID[1], ".csv"), ,
                             row.names=FALSE)
   )