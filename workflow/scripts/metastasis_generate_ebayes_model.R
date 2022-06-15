suppressMessages(library("edgeR"))
suppressMessages(library("limma"))

### SNAKEMAKE I/O ###
met_type_to_test      <- snakemake@input[["mets_to_test"]]
raw_gene_counts       <- snakemake@input[["raw_gene_counts"]]

ebayes_model          <- snakemake@output[["ebayes"]]

## Logging
log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")


## Load and set data types manually just in case
ccle_counts           <- readRDS(raw_gene_counts)
met_type_to_test      <- read.csv(met_type_to_test)

## Sanity check
met_type_to_test$mean            <- as.numeric(met_type_to_test$mean)
met_type_to_test$lineage        <- as.factor(met_type_to_test$lineage)

## Subset the counts
lines_to_test <- met_type_to_test$DepMap_ID
count_matrix  <- ccle_counts[, lines_to_test]

rownames(met_type_to_test) <- met_type_to_test$DepMap_ID

## voom model
design <- model.matrix(~lineage + mean, data=compound_to_test)

## reorder count_matrix so that cols matches rows from design
count_matrix <- count_matrix[ ,rownames(design)]

dge  <- DGEList(counts=count_matrix)
keep <- filterByExpr(dge, design=design)

dge <- dge[keep, keep.lib.sizes=FALSE]
dge <- calcNormFactors(dge)

## TODO: save voom model plot
v <- voom(dge, plot=FALSE, normalize.method="quantile")

fit <- lmFit(v, design)
fit <- eBayes(fit)

saveRDS(fit, ebayes_model)





