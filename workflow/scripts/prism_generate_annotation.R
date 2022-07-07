library("tidyverse")


### SNAKEMAKE I/O ###
response_curves        <- snakemake@input[["response_curves"]]
cell_lines_annotation  <- snakemake@input[["cell_lines_annotation"]]
count_matrix           <- snakemake@input[["count_matrix"]]

compounds_lines_profiled   <- snakemake@output[["compounds_lines_profiled"]]
auc_models_candidates      <- snakemake@output[["auc_models_candidates"]]

## load tables
response_curves <- data.table::fread(response_curves) %>%
                  as.data.frame()

cell_lines_annotation <- data.table::fread(cell_lines_annotation) %>%
                         as.data.frame()

ccle_counts <- readRDS(count_matrix)

## Filter out LOW QC response curves
## MTS are technical redos of some assays
## Keep HTS for now, in the future maybe mix, as it seems MTS-10 is better
response_curves <- response_curves[!is.na(response_curves$auc),]

response_curves <- response_curves %>%
    filter(screen_id == "HTS002" & !is.na(auc) & r2 >= 0.7 & lower_limit <= 1)

## Get rid of cell lines and curves for which no RNASeq profiling was done
available_lines <- intersect(colnames(ccle_counts), response_curves$depmap_id)
response_curves <- response_curves[response_curves$depmap_id %in% available_lines, ]

## Get the nº. lines/compound and the CV
lines_and_compounds <- response_curves %>%
    group_by(broad_id) %>%
    summarise(
        profiled_lines = n_distinct(depmap_id),
         cv = sd(auc) / mean(auc),
    ) %>%
    as.data.frame()

write.csv(lines_and_compounds, file = compounds_lines_profiled, row.names=FALSE)

## Keep compounds with at least 10 profiled lines and CV >= 0.1
compounds_to_test <- lines_and_compounds %>%
    filter(profiled_lines >= 10 & cv >= 0.1) %>%
    pull(broad_id)

response_curves   <- filter(response_curves, broad_id %in% compounds_to_test)

## Annotate cell line histology
response_curves <- merge(x=response_curves,
                         y=cell_lines_annotation[,c("DepMap_ID", "lineage")],
                         by.x="depmap_id",
                         by.y="DepMap_ID",
                         all.x=TRUE,
                         all.y=FALSE)

response_curves$lineage <- as.factor(response_curves$lineage)

if(!dir.exists(file.path(auc_models_candidates))){
    dir.create(auc_models_candidates)}

by(response_curves,
   response_curves$broad_id,
   FUN=function(i) write.csv(i, paste0(auc_models_candidates, "/", i$broad_id[1], ".csv"), , row.names=FALSE))
