library('tidyverse')


### SNAKEMAKE I/O ###
response_curves <- snakemake@input[['response_curves']]
celligner_data  <- snakemake@input[['celligner_data']]
count_matrix    <- snakemake@input[['count_matrix']]

cell_line_annotation       <- snakemake@output[['cell_line_annotation']]
compounds_lines_profiled   <- snakemake@output[['compounds_lines_profiled']]
auc_models_candidates      <- snakemake@output[['auc_models_candidates']]


### FUNCTIONS DEFINITIONS ###
save_candidate_compound  <- function(df){
write.csv(df,paste0(auc_models_candidates, '/',unique(df$broad_id),".csv"))
return(df)
}

## load tables
response_curves <- data.table::fread(response_curves) %>%
                  as.data.frame()

celligner_data <- data.table::fread(celligner_data) %>%
                  as.data.frame()


ccle_counts <- readRDS(count_matrix)

# Keep cell lines data exclusively
is_line <- grepl(pattern='^ACH-', x=celligner_data$sampleID, fixed=FALSE, perl=TRUE)

celligner_data <- celligner_data[is_line, c('sampleID', 'sampleID_CCLE_Name', 
                                            'CL_tumor_class', 'undifferentiated_cluster')]

rownames(celligner_data) <- celligner_data$sampleID
celligner_data$sampleID_CCLE_Name <- as.factor(celligner_data$sampleID_CCLE_Name)

# Set lines from the undif. cluster as undifferentiated
celligner_data[celligner_data$undifferentiated_cluster, 'CL_tumor_class'] <- 'undifferentiated'
celligner_data$CL_tumor_class <- as.factor(celligner_data$CL_tumor_class)

## I'm not really sure we need this anymore, but It might be useful downstream so
## I'm saving it for now
write.csv(x=celligner_data, file=cell_line_annotation)

## MTS are technical redos of some assays
## Keep HTS for now, in the future maybe mix, as it seems MTS-10 is better
response_curves <- response_curves[response_curves$screen_id == 'HTS002',]

## Keep compounds for which we have AUC data
response_curves <- response_curves[!is.na(response_curves$auc),]

## Get rid of models with low performance
response_curves <- response_curves[response_curves$r2 >= 0.6,]

## Get rid of cell lines and curves for which no RNASeq profiling was done
available_lines <- intersect(colnames(ccle_counts), response_curves$depmap_id)
response_curves <- response_curves[response_curves$depmap_id %in% available_lines,]

## Get the nº. lines/compound
lines_and_compounds <- response_curves %>%
    group_by(broad_id) %>%
    summarise(profiled_lines = n_distinct(depmap_id)) %>%
    as.data.frame()

write.csv(lines_and_compounds, file = compounds_lines_profiled, row.names=FALSE)

## Keep compounds with at least 10 profiled lines
compounds_to_test <- lines_and_compounds[lines_and_compounds$profiled_lines >= 10, 'broad_id']
response_curves   <- filter(response_curves, broad_id %in% compounds_to_test)

## Annotate cell line histology
response_curves <- merge(x=response_curves,
                         y=celligner_data[,c('sampleID', 'CL_tumor_class')],
                         by.x='depmap_id',
                         by.y='sampleID',
                         all.x=TRUE,
                         all.y=FALSE)

response_curves$CL_tumor_class <- as.factor(response_curves$CL_tumor_class)

if(!dir.exists(file.path(auc_models_candidates))){
    dir.create(auc_models_candidates)}

by(response_curves,
   response_curves$broad_id,
   FUN=function(i) write.csv(i, paste0(auc_models_candidates, '/', i$broad_id[1], ".csv"), , row.names=FALSE))
