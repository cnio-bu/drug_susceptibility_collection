rule download_cel_files:
    input:
    output:
        raw_cel_files=directory(f'{results}/gdsc/array_data/raw'),
    threads: 1
    resources:
        mem_mb=1024
    conda:
        '../envs/gdsc_arrayexpress.yaml'
    script:
        '../scripts/gdsc_download_raw_files.R'


rule normalize_gdsc_arrays:
    input:
        cel_files= rules.download_cel_files.output.raw_cel_files
    output:
        normalized_arrays=f'{results}/gdsc/array_data/normalized_arrays.rds'
    threads: 1
    resources:
        mem_mb=20480
    conda:
        '../envs/gdsc_array_normalization.yaml'
    script:
        '../scripts/gdsc_normalize_arrays.R'
