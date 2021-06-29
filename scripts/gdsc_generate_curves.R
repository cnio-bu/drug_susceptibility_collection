
library('tidyverse')

### SNAKEMAKE I/O ###
gdsc_dose_response_curves      <- snakemake@input[['dose_response_curves']]
normalized_lines_arrays        <- snakemake@input[['normalized_arrays']]
cell_lines_annotation          <- snakemake@input[['cell_lines_annotation']]

auc_models_candidates          <- snakemake@output[['auc_models_candidates']]
compounds_lines_profiled       <- snakemake@output[['compounds_lines_profiled']]

# Setup log file
log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

candidate_curves <- readxl::read_xlsx(gdsc_dose_response_curves) %>%
                    filter(RMSE <= 0.3) %>%
                    as.data.frame()

cell_lines_annotation <- data.table::fread(cell_lines_annotation) %>%
                        as.data.frame()

normalized_arrays <- readRDS(normalized_lines_arrays)


## Filter out drug/lines associations where the line is not profiled
available_lines <- colnames(normalized_arrays)
rm(normalized_arrays)

## Match SANGER and CCLE info
## I Manually tested it and only 1 line cannot be matched
sanger_lines <- unique(candidate_curves$SANGER_MODEL_ID)
ccle_lines   <- unique(cell_lines_annotation$Sanger_Model_ID)

common_lines <- intersect(sanger_lines, ccle_lines)

candidate_curves <- filter(candidate_curves, SANGER_MODEL_ID %in% common_lines)

## Annotate cell line histology
candidate_curves <- merge(x=candidate_curves,
                         y=cell_lines_annotation[,c('Sanger_Model_ID', 'lineage')],
                         by.x='SANGER_MODEL_ID',
                         by.y='Sanger_Model_ID',
                         all.x=TRUE,
                         all.y=FALSE)

## get the number of profiled lines/compound
lines_by_compound <- candidate_curves %>%
                     filter(CELL_LINE_NAME %in% available_lines) %>%
                     group_by(DRUG_ID) %>%
                     summarise(profiled_lines = n_distinct(SANGER_MODEL_ID)) %>%
                     as.data.frame()

##TODO: SAVE THIS
write.csv(lines_by_compound, file=compounds_lines_profiled, row.names=FALSE)

## Get rid of models with less than 10 profiled lines
compounds_to_test  <- lines_by_compound[lines_by_compound$profiled_lines >= 10, 'DRUG_ID']
candidate_curves   <- filter(candidate_curves, DRUG_ID %in% compounds_to_test)

candidate_curves$lineage <- as.factor(candidate_curves$lineage)

if(!dir.exists(file.path(auc_models_candidates))){
    dir.create(auc_models_candidates)}

by(candidate_curves,
   candidate_curves$DRUG_ID,
   FUN=function(i) write.csv(i, paste0(auc_models_candidates, '/', i$DRUG_ID[1], ".csv"), , row.names=FALSE))
