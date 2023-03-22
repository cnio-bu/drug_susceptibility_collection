library('tidyverse')

### SNAKEMAKE I/O ###
raw_expected_counts <- snakemake@input[["raw_expected_counts"]]
cell_line_info <- snakemake@input[["cell_line_info"]]
protein_coding_genes <- snakemake@input[["protein_coding_genes"]]
raw_gene_counts <- snakemake@output[["raw_gene_counts"]]

### SNAKEMAKE PARAMS ###
coding_genes <- snakemake@params["coding_genes_only"]

## Load CCLE raw counts
ccle_counts <- data.table::fread(raw_expected_counts)

## Load cell line equivalences
cell_lines <- data.table::fread(cell_line_info)

## Keep RNA cell line equivalences
cell_lines <- cell_lines %>%
  filter(ProfileType == "rna") %>%
  rename(V1 = ProfileID) %>%
  select(V1, ModelID) %>%
  unique()

## Convert CCLE raw counts to a matrix where each gene is a row and each sample, 
## a column.
## Input is RSEM expected counts (thus I have to round up to units)
genes <- colnames(ccle_counts)[-1]
ccle_counts <- ccle_counts %>%
  merge(cell_lines, by = "V1") %>%
  column_to_rownames("ModelID") %>%
  select(-V1) %>%
  t() %>%
  as.data.frame()
rownames(ccle_counts) <- genes

## Most of the genes with NA values are ERCC-0000X spike-ins
ccle_counts <- na.omit(ccle_counts)

## BROAD format for gene names is HUGO (ENSEMBL)
ccle_counts <- ccle_counts %>%
  mutate(gene_variance = unname(apply(ccle_counts, MARGIN = 1, FUN = var)),
         Gene = str_remove(rownames(ccle_counts), pattern = " \\(.*"),
         Ensembl = str_extract(rownames(ccle_counts), 
                               pattern = "ENSG[0-9]+")) %>%
  ## Get rid of remaining spikes
  filter(!grepl(Gene, pattern = "ERCC-", fixed = TRUE))

## If coding_genes = TRUE, just keep protein coding genes
if (coding_genes) {
  hgnc_coding_genes <- data.table::fread(protein_coding_genes) %>%
    rename(Gene = symbol, Ensembl = ensembl_gene_id) %>%
    select(Gene, Ensembl) %>%
    unique()
  ccle_counts <- ccle_counts %>%
    select(-Gene) %>%
    merge(hgnc_coding_genes)
}

## Keep the most variable gene when two ENSEMBL ids point to the same HGNC
ccle_counts <- ccle_counts %>%
  group_by(Gene) %>%
  filter(gene_variance == max(gene_variance)) %>%
  column_to_rownames("Gene") %>%
  select(-gene_variance, -Ensembl)

## Round counts
ccle_counts <- round(as.matrix(ccle_counts), digits = 0)

## Save object
saveRDS(ccle_counts, file = raw_gene_counts)
