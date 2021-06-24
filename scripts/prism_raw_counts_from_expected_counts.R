library('dplyr')

### SNAKEMAKE I/O ###
raw_expected_counts <- snakemake@input[['raw_expected_counts']]
raw_gene_counts     <- snakemake@output[['raw_gene_counts']]

## load CCLE raw counts
ccle_counts <- data.table::fread(raw_expected_counts) %>%
    as.data.frame()

## Convert it to a matrix where each gene is a row and each sample, a column.
## Input is RSEM expected counts (thus I have to round up to units)
rownames(ccle_counts) <- ccle_counts$V1
ccle_counts$V1        <- NULL
ccle_counts           <- t(ccle_counts)
rownames(ccle_counts) <- sub(pattern = " \\(.*", replacement = '',
                             x=rownames(ccle_counts))

ccle_counts           <- na.omit(ccle_counts)
ccle_counts           <- round(ccle_counts, digits=0)
ccle_counts           <- as.matrix(ccle_counts)

saveRDS(object=ccle_counts, file=raw_gene_counts)
