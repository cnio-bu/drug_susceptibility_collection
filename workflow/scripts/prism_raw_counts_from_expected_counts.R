suppressMessages(library('tidyverse'))

### SNAKEMAKE I/O ###
raw_expected_counts <- snakemake@input[["raw_expected_counts"]]
protein_coding_genes <- snakemake@input[["protein_coding_genes"]]
raw_gene_counts <- snakemake@output[["raw_gene_counts"]]

### SNAKEMAKE PARAMS ###
coding_genes <- as.logical(snakemake@params["coding_genes_only"])

## Load CCLE raw counts
ccle_counts <- data.table::fread(raw_expected_counts)

## Convert CCLE raw counts to a matrix where each gene is a row and each sample, a column.
## Input is RSEM expected counts (at gene level derived from unstranded RNA-seq for all genes in human)
# Keep default entry for model and remove unnecessary columns
ccle_counts <- ccle_counts %>%
  filter(IsDefaultEntryForModel == "Yes") %>%
  select(-1, -c(SequencingID, IsDefaultEntryForModel, ModelConditionID, IsDefaultEntryForMC)) %>%
  column_to_rownames("ModelID") %>%
  t() %>%
  as.data.frame()

## Most of the genes with NA values are ERCC-0000X spike-ins
ccle_counts <- na.omit(ccle_counts)

## BROAD format for gene names is HUGO (ENSEMBL)
ccle_counts <- ccle_counts %>%
  mutate(gene_variance = unname(apply(ccle_counts, MARGIN = 1, FUN = var)),
         Gene = str_remove(rownames(ccle_counts), pattern = " \\(.*"),
         Entrez = as.character(str_extract(rownames(ccle_counts), pattern = "(?<=\\()\\d+(?=\\))")),
         Ensembl = str_extract(rownames(ccle_counts), 
                               pattern = "ENSG[0-9]+")) %>%
  ## Get rid of remaining spikes
  filter(!grepl(Gene, pattern = "ERCC-", fixed = TRUE))

## If coding_genes = TRUE, just keep protein coding genes
if (coding_genes) {
  hgnc_coding_genes <- data.table::fread(protein_coding_genes) %>%
    rename(Gene = symbol, Ensembl = ensembl_gene_id, Entrez = entrez_id) %>%
    mutate(Entrez = as.character(Entrez)) %>%
    select(Gene, Ensembl, Entrez) %>%
    unique()
  ccle_counts <- bind_rows(
    ccle_counts %>%
      filter(!is.na(Entrez)) %>%
      select(-c(Gene, Ensembl)) %>%
      merge(hgnc_coding_genes, by = "Entrez"),
    
    ccle_counts %>%
      filter(is.na(Entrez)) %>%
      select(-c(Gene, Entrez)) %>%
      merge(hgnc_coding_genes, by = "Ensembl")
)
}

## Keep the most variable gene when two ENSEMBL ids point to the same HGNC
ccle_counts <- ccle_counts %>%
  group_by(Gene) %>%
  filter(gene_variance == max(gene_variance)) %>%
  column_to_rownames("Gene") %>%
  select(-gene_variance, -Entrez, -Ensembl)

## Round counts
ccle_counts <- as.matrix(ccle_counts)

## Save object
saveRDS(ccle_counts, file = raw_gene_counts)
