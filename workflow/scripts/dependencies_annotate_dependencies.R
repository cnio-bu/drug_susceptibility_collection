library("tidyverse")

###  SNAKEMAKE I/O ###
raw_crispr_data = snakemake@input[["crispr_gene_dependency_chronos"]]
lines_info = snakemake@input[["sample_info"]]
ccle_raw_reads = snakemake@input[["expression_matrix"]]
where_to_save = snakemake@output[["model_candidates"]]

ccle_counts <- readRDS(ccle_raw_reads)
crispr <- data.table::fread(raw_crispr_data)
annot <-  data.table::fread(lines_info)

## keep lines whose primary disease is cancer AND with available expr.
annot_cancer <- annot %>%
    filter(
        !primary_disease %in% c("Unknown", "Non-Cancerous") 
    )

## drop lineages with too few CCLs
all_lineages <- table(annot_cancer$lineage)
lineages_to_keep <- names(all_lineages[all_lineages >= 5])

## get the names of the CCLs whose lineages are OK AND with RNA-seq
lines_to_keep <- annot_cancer %>%
    filter(
        lineage %in% lineages_to_keep
    ) %>%
    pull(DepMap_ID)

lines_to_keep <- intersect(lines_to_keep, colnames(ccle_counts))

## Filter crispr data
crispr_filtered <- crispr %>%
    filter(
        ModelID %in% lines_to_keep
    )

crispr_lines <- crispr_filtered$ModelID
crispr_filtered$ModelID <- NULL

# CCLE nomenclature is HUGO (ENTREZ). Keep HUGO ONLY
colnames(crispr_filtered) <- stringr::str_remove_all(string = colnames(crispr_filtered),
                                                     pattern = " \\(.*$"
                                                     )

## Resolve duplicates keeping most variable gene
crispr_filtered <- t(crispr_filtered)
colnames(crispr_filtered) <- crispr_lines

vars <- apply(crispr_filtered, 1, var)
crispr_filtered <- cbind(crispr_filtered, vars)
crispr_filtered <- crispr_filtered[order(crispr_filtered[, "vars"],decreasing=TRUE),]

crispr_filtered_unique <- crispr_filtered[!duplicated(rownames(crispr_filtered)), ]
crispr_filtered_unique <- crispr_filtered_unique[, colnames(crispr_filtered_unique) != "vars"]

## Melt the dataset
crispr_probabilities <- crispr_filtered_unique %>%
    as_tibble(rownames = "Gene") %>%
    pivot_longer(
        cols = starts_with("ACH-"),
        names_to = "cell_line",
        values_to = "probability"
        )


# Filter out genes with less than 5 cell lines dependant on them
crispr_enough_power <- crispr_probabilities %>%
    group_by(Gene) %>%
    mutate(dependent_lines = sum(probability >= 0.5)) %>%
    filter(
        dependent_lines >= 5
    )

# annotate the lineage and stripped name
crispr_enough_power_annotated <- crispr_enough_power %>%
    left_join(
        y = annot_cancer[, c("DepMap_ID",
                             "lineage",
                             "stripped_cell_line_name",
                             "original_lineage"
                             )
                         ],
        by = c("cell_line" = "DepMap_ID")
    )


if(!dir.exists(file.path(where_to_save))){
    dir.create(where_to_save)}

by(crispr_enough_power_annotated,
   crispr_enough_power_annotated$Gene,
   FUN=function(i) write.csv(i, paste0(where_to_save, "/", i$Gene[1], ".csv"), , row.names=FALSE))
