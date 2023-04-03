library("limma")
library("fgsea")

### SNAKEMAKE I/O ###
ebayes_model    <- snakemake@input[["ebayes_model"]]
geneset_collection <- snakemake@input[["drugs"]]
enrichment_table  <- snakemake@output[["enrichments"]]

## Load and set data types manually just in case
ebayes_model <- readRDS(ebayes_model)
geneset_collection <- fgsea::gmtPathways(geneset_collection)

## Subset the collection to keep UP drug genesets only
is_up_geneset <- grepl(pattern = "\\.*_UP", x = names(geneset_collection))
geneset_collection <- geneset_collection[is_up_geneset]

## generate a named vector of gene-level stats
## An actual ranking is not needed
t_rank <- topTable(
    fit = ebayes_model,
    coef = "probability",
    number = Inf,
    sort.by = "t"
    )

named_gene_rank <- t_rank$t
names(named_gene_rank) <- rownames(t_rank)

fgsea_res <- fgsea(
    pathways = geneset_collection,
    stats = named_gene_rank,
    minSize = 15,
    maxSize = 500,
    gseaParam = 1,
    scoreType = "pos"
)

fgsea_res$leadingEdge <- sapply(
    fgsea_res$leadingEdge,
    FUN = paste,
    collapse = " "
    )

new_directory <- dirname(enrichment_table)

if (!dir.exists(file.path(new_directory))) {
    dir.create(new_directory)
    }

write.table(
    x = fgsea_res,
    file = enrichment_table,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE
)
