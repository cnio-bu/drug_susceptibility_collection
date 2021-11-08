library('tidyverse')

### SNAKEMAKE I/O ###
dose_response_curves           <- snakemake@input[['curves_data']]
compound_meta                  <- snakemake@input[['compount_meta']]
experiment_meta                <- snakemake@input[['experiment_meta']]
cell_meta                      <- snakemake@input[['cell_meta']]
cell_lines_annotation          <- snakemake@input[['cell_lines_annotation']]
count_matrix                   <- snakemake@input[['rna_count_matrix']]

auc_models_candidates          <- snakemake@output[['auc_models_candidates']]
compounds_lines_profiled       <- snakemake@output[['compounds_lines_profiled']]

# Setup log file
log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

## Load tables
drug_data  <- data.table::fread(dose_response_curves) %>%
    as.data.frame()

drug_meta <- data.table::fread(compound_meta) %>%
    as.data.frame()

experiment_meta <- data.table::fread(experiment_meta) %>%
    as.data.frame()

# This table and the below one seem redundant, but CTRP data might be outdated.
# We are loading it ONLY to translate experiment_ids to cell line IDs to match
# to ccle_info
cell_meta <- data.table::fread(cell_meta) %>%
    as.data.frame()

ccle_info <- data.table::fread(cell_lines_annotation) %>%
    as.data.frame() ## I don't like tibbles


## Load count matrix
ccle_counts <- readRDS(count_matrix)

## Add broad_cpd_id to compounds id
drug_data_annotated <- merge(x=drug_data,
                             y=drug_meta[,c('master_cpd_id',  'broad_cpd_id')],
                             by.x='master_cpd_id',
                             by.y = 'master_cpd_id',
                             all.x=TRUE,
                             all.y=FALSE)


## create a named vector of master ccle ids
master_ids <- experiment_meta[!duplicated(experiment_meta$experiment_id), c('experiment_id', 'master_ccl_id')] 

drug_data_annotated$master_ccl_id <- master_ids[drug_data_annotated$experiment_id, 'master_ccl_id']

## Use master_ccl_id to  annotated drug_data with ACH depmap ids
drug_data_annotated$human_ccl_name <- cell_meta[drug_data_annotated$master_ccl_id, 'ccl_name']

# get rid of experiments without ccl_name. These are collaborators lines and
# thus no RNA-Seq data is available.
drug_data_annotated_lines <- drug_data_annotated[!is.na(drug_data_annotated$human_ccl_name),]

# now get the ACH ids for the remaining lines
annotations_fields_of_interest <- c('DepMap_ID', 'stripped_cell_line_name', 'lineage', 'is_undifferentiated')

drug_data_annotated_lines_depmap <- merge(x=drug_data_annotated_lines,
                                          y=ccle_info[,annotations_fields_of_interest],
                                          by.x='human_ccl_name',
                                          by.y='stripped_cell_line_name',
                                          all.x=TRUE)


## Get rid of lines withou ACH-id If somehow there is any missing left
drug_data_annotated_lines_depmap <- drug_data_annotated_lines_depmap[!is.na(drug_data_annotated_lines_depmap$DepMap_ID),]

## Get rid of compounds without area_under_curve
drug_data_annotated_lines_depmap <- drug_data_annotated_lines_depmap[!is.na(drug_data_annotated_lines_depmap$area_under_curve),]

## Get rid of cell lines and curves for which no RNASeq profiling was done
available_lines <- intersect(colnames(ccle_counts), drug_data_annotated_lines_depmap$DepMap_ID)
drug_data_annotated_lines_depmap <- drug_data_annotated_lines_depmap[drug_data_annotated_lines_depmap$DepMap_ID %in% available_lines,]

# If there are duplicate cell lines by compound, keep the one with the lowest auc
drug_data_annotated_lines_depmap <- drug_data_annotated_lines_depmap[order(drug_data_annotated_lines_depmap$area_under_curve),]
drug_line_duplicated             <- duplicated(drug_data_annotated_lines_depmap[,c('broad_cpd_id', 'DepMap_ID')])
drug_data_annotated_lines_depmap <- drug_data_annotated_lines_depmap[!drug_line_duplicated,]

## Get the nº. lines/compound
lines_and_compounds <- drug_data_annotated_lines_depmap %>%
    group_by(broad_cpd_id) %>%
    summarise(profiled_lines = n_distinct(DepMap_ID)) %>%
    as.data.frame()

# Save this table. It is definitely useful
write.csv(lines_and_compounds, file = compounds_lines_profiled, row.names=FALSE)

## Keep compounds with at least 10 profiled lines
compounds_to_test <- lines_and_compounds[lines_and_compounds$profiled_lines >= 10, 'broad_cpd_id']
drug_data_annotated_lines_depmap <- filter(drug_data_annotated_lines_depmap,  broad_cpd_id %in% compounds_to_test)

## set lineage if undiff
drug_data_annotated_lines_depmap[drug_data_annotated_lines_depmap$is_undifferentiated, 'lineage'] <- 'undifferentiated'

drug_data_annotated_lines_depmap$lineage          <- as.factor(drug_data_annotated_lines_depmap$lineage)
drug_data_annotated_lines_depmap$area_under_curve <- as.numeric(drug_data_annotated_lines_depmap$area_under_curve)

if(!dir.exists(file.path(auc_models_candidates))){
    dir.create(auc_models_candidates)}

by(drug_data_annotated_lines_depmap,
   drug_data_annotated_lines_depmap$broad_cpd_id,
   FUN=function(i) write.csv(i, paste0(auc_models_candidates, '/', i$broad_cpd_id[1], ".csv"), , row.names=FALSE))
