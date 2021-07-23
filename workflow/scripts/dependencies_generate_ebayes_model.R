library('edgeR')
library('limma')

### SNAKEMAKE I/O ###
dep_to_test           <- snakemake@input[['dependency_to_test']]
raw_gene_counts       <- snakemake@input[['raw_gene_counts']]

ebayes_model          <- snakemake@output[['ebayes']]

## Load and set data types manually just in case
ccle_counts           <- readRDS(raw_gene_counts)
dep_to_test           <- read.csv(dep_to_test)

dep_to_test$probability         <- as.numeric(dep_to_test$probability)
dep_to_test$lineage             <- as.factor(dep_to_test$lineage)

## Subset the counts
lines_to_test <- dep_to_test$cell_line
count_matrix  <- ccle_counts[,lines_to_test]

rownames(dep_to_test) <- dep_to_test$cell_line

## voom model
design <- model.matrix(~lineage + probability, data=dep_to_test)

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





