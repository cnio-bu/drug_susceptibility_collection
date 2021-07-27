library('dplyr')

### SNAKEMAKE I/O ###
raw_expected_counts <- snakemake@input[['raw_expected_counts']]
raw_gene_counts     <- snakemake@output[['raw_gene_counts']]

## load CCLE raw counts
ccle_counts <- data.table::fread(raw_expected_counts) %>%
    as.data.frame() ## I don't like tibbles

## Convert it to a matrix where each gene is a row and each sample, a column.
## Input is RSEM expected counts (thus I have to round up to units)
rownames(ccle_counts) <- ccle_counts$V1
ccle_counts$V1        <- NULL
ccle_counts           <- t(ccle_counts)

## as data.frame again because of the transposition
ccle_counts           <- as.data.frame(ccle_counts)
ccle_counts$Gene      <- rownames(ccle_counts)

# Most of the genes with NA values are ERCC-0000X spike-ins
ccle_counts           <- na.omit(ccle_counts)

# BROAD format for gene names is HUGO (ENSEMBL)
ccle_counts$Gene<- sub(pattern = " \\(.*", replacement = '', x=ccle_counts$Gene)

## Keep the most variable gene when two ENSEMBL ids point to the same HGNC
ccle_counts$gene_variance <- apply(X=ccle_counts[,colnames(ccle_counts)!='Gene'], MARGIN=1, FUN=var)

## sort by variance in descending order and keep the first of the duplicated
## entries (the one with the highest variance)
ccle_counts <- ccle_counts[order(-ccle_counts$gene_variance),]
ccle_counts <- ccle_counts[!duplicated(ccle_counts$Gene, fromLast=FALSE),]

# Get rid of the variance column
ccle_counts$gene_variance <- NULL

# Set the genes as rownames now that they are not duplicated and get rid 
# of the column
rownames(ccle_counts) <- ccle_counts$Gene
ccle_counts$Gene      <- NULL

ccle_counts           <- round(ccle_counts, digits=0)
ccle_counts           <- as.matrix(ccle_counts)

saveRDS(object=ccle_counts, file=raw_gene_counts)



