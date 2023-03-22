library('tidyverse')

### SNAKEMAKE I/O ###
raw_expected_counts <- snakemake@input[['raw_expected_counts']]
raw_gene_counts     <- snakemake@output[['raw_gene_counts']]

# load CCLE raw counts
ccle_counts <- data.table::fread(raw_expected_counts)

# Convert it to a matrix where each gene is a row and each sample, a column.
# Input is RSEM expected counts (thus I have to round up to units)
cell_lines <- colnames(ccle_counts)[-1]
ccle_counts <- ccle_counts %>%
  column_to_rownames("V1") %>%
  t() %>%
  as.data.frame()
rownames(ccle_counts) <- cell_lines

# Most of the genes with NA values are ERCC-0000X spike-ins
ccle_counts <- na.omit(ccle_counts)

# BROAD format for gene names is HUGO (ENSEMBL)
ccle_counts <- ccle_counts %>%
  mutate(gene_variance = unname(apply(ccle_counts, MARGIN = 1, FUN = var)),
         Gene = str_remove(rownames(ccle_counts), pattern = " \\(.*")) %>%
  # Get rid of remaining spikes
  filter(!grepl(Gene, pattern = "ERCC-", fixed = TRUE))

# Keep the most variable gene when two ENSEMBL ids point to the same HGNC
ccle_counts <- ccle_counts %>%
  group_by(Gene) %>%
  filter(gene_variance == max(gene_variance)) %>%
  column_to_rownames("Gene") %>%
  select(-gene_variance)

# Round counts
ccle_counts <- round(as.matrix(ccle_counts), digits = 0)

# Save object
saveRDS(ccle_counts, file = raw_gene_counts)
