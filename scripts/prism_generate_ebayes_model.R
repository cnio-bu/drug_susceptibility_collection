library('edgeR')
library('limma')

### SNAKEMAKE I/O ###
compound_to_test      <- snakemake@input[['compound_to_test']]
raw_gene_counts       <- snakemake@input[['raw_gene_counts']]

ebayes_model          <- snakemake@output[['ebayes']]

## Load and set data types manually just in case
ccle_counts           <- readRDS(raw_gene_counts)
compound_to_test      <- read.csv(compound_to_test)

compound_to_test$auc            <- as.numeric(compound_to_test$auc)
compound_to_test$CL_tumor_class <- as.factor(compound_to_test$CL_tumor_class)

## Subset the counts
lines_to_test <- compound_to_test$depmap_id
count_matrix  <- ccle_counts[,lines_to_test]

rownames(compound_to_test) <- compound_to_test$depmap_id

## voom model
design <- model.matrix(~CL_tumor_class + auc, data=compound_to_test)

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





