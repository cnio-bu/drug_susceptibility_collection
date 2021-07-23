library('limma')

### SNAKEMAKE I/O ###
compound_to_test        <- snakemake@input[['compound_to_test']]
normalized_arrays       <- snakemake@input[['normalized_arrays']]
ebayes_model            <- snakemake@output[['ebayes']]

## Load and set data types manually just in case
normalized_arrays           <- readRDS(normalized_arrays)
compound_to_test            <- read.csv(compound_to_test)

compound_to_test$AUC            <- as.numeric(compound_to_test$AUC)
compound_to_test$lineage        <- as.factor(compound_to_test$lineage)


## Subset the arrays
lines_to_test      <- compound_to_test$CELL_LINE_NAME
normalized_arrays  <- normalized_arrays[,lines_to_test]

rownames(compound_to_test) <- compound_to_test$CELL_LINE_NAME

design <- model.matrix(~lineage + AUC, data=compound_to_test)

## reorder count_matrix so that cols matches rows from design
normalized_arrays <- normalized_arrays[,rownames(design)]

fit <- lmFit(normalized_arrays, design)
fit <- eBayes(fit)

saveRDS(fit, ebayes_model)
