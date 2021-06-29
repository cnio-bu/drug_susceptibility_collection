library('oligo')
library('oligoClasses')

### SNAKEMAKE I/O ###
raw_cel_path      <- snakemake@input[['cel_files']]
normalized_arrays <- snakemake@output[['normalized_arrays']]

# Setup log file
log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

sdrf_location <- file.path(raw_cel_path, "E-MTAB-3610.sdrf.txt")
SDRF          <- read.delim(sdrf_location)


rownames(SDRF) <- SDRF$Array.Data.File
SDRF           <- AnnotatedDataFrame(SDRF)

raw_exprset <- oligo::read.celfiles(filenames = file.path(raw_cel_path, 
                                    SDRF$Array.Data.File),
                                    verbose = FALSE, phenoData = SDRF)


norm_eset <- oligo::rma(raw_exprset, background = TRUE, normalize = TRUE)

## Save it as a matrix directly, no need to keep the whole object.
## Make sure to change colnames 
pdata <- phenoData(norm_eset) 
pdata <- as.data.frame(pdata@data)

eset            <- exprs(norm_eset)
colnames(eset)  <- pdata[colnames(eset), 'Characteristics.cell.line.']

saveRDS(eset, normalized_arrays)