library('ArrayExpress')

### SNAKEMAKE I/O ###
raw_cel_path      <- snakemake@output[['raw_cel_files']]

if(!dir.exists(file.path(raw_cel_path))){
    dir.create(raw_cel_path)}

# E-MTAB-3610 is the ArrayExpress ID of the GDSC-1000 profiling array
array_set = ArrayExpress::getAE('E-MTAB-3610', path=raw_cel_path, type='raw')