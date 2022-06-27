log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

suppressMessages(library("limma"))
suppressMessages(library("openxlsx"))
suppressMessages(library("tidyverse"))
suppressMessages(library("GSEABase"))

## SNAKEMAKE I/O ##
eBayes_model <- snakemake@input[["fitted_bayes"]]
drug_info <- snakemake@input[["treatment_info"]]
hgnc <- snakemake@input[["hgnc"]]

geneset_directory <- snakemake@output[["bidirectional_geneset"]]

## SNAKEMAKE PARAMS ##
signature_type <- tolower(as.character(snakemake@params[['signature_type']]))

## FUNCTION ##
create.gmt <- function(x, mode) {
  gset <- GeneSet(x, geneIdType = SymbolIdentifier())
  setName(gset) <- paste(sig_name, mode, sep = "_")
  collapsed_name <- str_replace_all(sig_name, pattern = " ", repl = "")
  toGmt(gset, con = paste0(geneset_directory, '/', collapsed_name, "_", mode, 
                           ".gmt"))
}

## CODE ##
# Create output directory
dir.create(geneset_directory, showWarnings = FALSE)

# Get Bayes result
bayes_results <- readRDS(eBayes_model)

# Get treatment annotation
drug_info <- read.xlsx(drug_info)
drug_info <- drug_info %>% dplyr::select(DRUG_ID, DRUG_NAME) %>% unique

# Get signature name
sigID <- gsub("_eBayes.rds", "", basename(eBayes_model))
common_name <- drug_info$DRUG_NAME[match(sigID, drug_info$DRUG_ID, 
                                         nomatch = sigID)]
sig_name <- paste(common_name, "GDSC", sigID, sep = "_")

# Get probe ID and HUGO symbol correspondence
hgnc <- read.table(hgnc, sep = "\t", header = TRUE)

# Variables depending on signature_type
if (signature_type == "classic") {
  p.val = 1
  magnitude <- "t"
} else if (signature_type == "fold") {
  p.val = 0.05
  magnitude <- "logFC"
}

# Get top genes table
top_genes <- topTable(bayes_results, coef = "AUC", number = Inf, 
                      adjust.method = "fdr", p.value = p.val)

# Continue if the table is not empty
if (nrow(top_genes) != 0) {
  # Make the rownames a new column
  top_genes$probe <- rownames(top_genes)
  
  # Match with HUGO symbols
  top_genes <- na.omit(merge(top_genes, hgnc[, c("probe", "hgnc")], all.x = TRUE))
  
  # If there are several values for the same gene, keep the one with the highest 
  # absolute value of "magnitude"
  top_genes <- top_genes %>% group_by(hgnc) %>% 
    dplyr::slice(which.max(get(magnitude)))
  
  # Create bidirectional signature
  if (signature_type == "classic") {
    sensitivity <- top_genes %>% filter(t < 0)
    resistance <- top_genes %>% filter(t > 0)
  } else if (signature_type == "fold") {
    sensitivity <- top_genes %>% filter(logFC < 0)
    resistance <- top_genes %>% filter(logFC > 0)
  }
  
  # Keep the 500 top/bottom genes
  sensitivity <- head(sensitivity, n = 500) %>% arrange(desc(get(magnitude))) %>% 
    dplyr::select(hgnc) %>% pull
  resistance <- tail(resistance, n = 500) %>% arrange(get(magnitude)) %>%
    dplyr::select(hgnc) %>% pull
  
  # Ugly fix to keep fold signatures at 250. We really should refactor these
  if (signature_type == "fold"){
    sensitivity <- head(sensitivity, n = 250)
    resistance <- head(resistance, n = 250)
  }

  # Create a gmt if the number of genes > 15
  if (length(sensitivity) >= 15) {
    create.gmt(sensitivity, "UP")
  }
  if (length(resistance) >= 15) {
    create.gmt(resistance, "DOWN")
  }
}
