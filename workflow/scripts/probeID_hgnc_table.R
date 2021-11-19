log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

suppressMessages(library("hgu219.db"))
suppressMessages(library("biomaRt"))
suppressMessages(library("dplyr"))

## CODE ##
httr::set_config(httr::config(ssl_verifypeer = FALSE))

# Get ENSEMBL IDs for each probe
probeIDs <- unlist(as.list(hgu219ENSEMBL))
probeIDs <- na.omit(data.frame(probe = names(probeIDs), ensembl = probeIDs, 
                               row.names = NULL))

# BiomaRt
mart <- useMart("ENSEMBL_MART_ENSEMBL")
human <- useDataset("hsapiens_gene_ensembl", mart)

# Get HUGO genes
hugoIDs <- getBM(c("hgnc_symbol", "ensembl_gene_id"), mart = human,
                 filters = "ensembl_gene_id", values = unique(probeIDs$ensembl))
colnames(hugoIDs) <- c("hgnc", "ensembl")
hugoIDs <- hugoIDs %>% mutate(hgnc = na_if(hgnc, "")) %>% na.omit

# Merge tables
merged <- merge(probeIDs, hugoIDs, by = "ensembl")
merged <- unique(merged[c("probe", "hgnc", "ensembl")])

# Save
write.table(merged, file = snakemake@output[["tsv"]], sep = "\t", 
            row.names = FALSE, quote = FALSE)