library('ArrayExpress')
library('oligo')
library('preprocessCore')

### SNAKEMAKE I/O ###
raw_cel_path    <- snakemake@output[['cel_files']]
norm_array_path <- snakemake@output[['normalized_arrays']]

raw_data_dir <- '/raid/sagarcia/drug_signatures'

# E-MTAB-3610 is the ArrayExpress ID of the GDSC-1000 profiling array
array_set = ArrayExpress::getAE('E-MTAB-3610', path=raw_data_dir, type='raw')


sdrf_location <- file.path(raw_data_dir, "E-MTAB-3610.sdrf.txt")
SDRF          <- read.delim(sdrf_location)


rownames(SDRF) <- SDRF$Array.Data.File
SDRF           <- AnnotatedDataFrame(SDRF)

raw_exprset <- oligo::read.celfiles(filenames = file.path(raw_data_dir, 
                                    SDRF$Array.Data.File),
                                    verbose = FALSE, phenoData = SDRF)


## TADO: download pd.hg.u219_3.12.0

## test outliers
norm_eset <- oligo::rma(raw_exprset, normalize = TRUE)
