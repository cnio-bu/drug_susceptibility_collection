library('edgeR')
library('limma')

### SNAKEMAKE I/O ###
compound_to_test      <- snakemake@input[['compound_to_test']]
raw_gene_counts       <- snakemake@input[['raw_gene_counts']]

ebayes_model          <- snakemake@output[['ebayes']]

## Load and set data types manually just in case
ccle_counts           <- readRDS(raw_gene_counts)
compound_to_test      <- read.csv(compound_to_test)

# This might be the only key difference between this script and Dep/PRISM generation
compound_to_test$area_under_curve            <- as.numeric(compound_to_test$area_under_curve)
compound_to_test$lineage                     <- as.factor(compound_to_test$lineage)

## Subset the counts
lines_to_test <- compound_to_test$DepMap_ID
count_matrix  <- ccle_counts[,lines_to_test]

rownames(compound_to_test) <- compound_to_test$DepMap_ID

## voom model
design <- model.matrix(~lineage + area_under_curve, data=compound_to_test)

## reorder count_matrix so that cols matches rows from design
count_matrix <- count_matrix[,rownames(design)]

dge  <- DGEList(counts=count_matrix)
keep <- filterByExpr(dge, design=design)

dge <- dge[keep, keep.lib.sizes=FALSE]
dge <- calcNormFactors(dge)

## TODO: save voom model plot
v <- voom(dge, plot=FALSE, normalize.method='quantile')

fit <- lmFit(v, design)
fit <- eBayes(fit)

saveRDS(fit, ebayes_model)