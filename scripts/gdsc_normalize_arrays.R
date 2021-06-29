library('oligo')
library('oligoClasses')

### SNAKEMAKE I/O ###
raw_cel_path      <- snakemake@input[['cel_files']]
normalized_arrays <- snakemake@output[['normalized_arrays']]


sdrf_location <- file.path(raw_cel_path, "E-MTAB-3610.sdrf.txt")
SDRF          <- read.delim(sdrf_location)


rownames(SDRF) <- SDRF$Array.Data.File
SDRF           <- AnnotatedDataFrame(SDRF)

raw_exprset <- oligo::read.celfiles(filenames = file.path(raw_cel_path, 
                                    SDRF$Array.Data.File),
                                    verbose = FALSE, phenoData = SDRF)


## TADO: download pd.hg.u219_3.12.0
norm_eset <- oligo::rma(raw_exprset, background = TRUE, normalize = TRUE)
saveRDS(norm_eset, normalized_arrays)