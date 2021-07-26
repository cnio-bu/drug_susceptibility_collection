rule gdsc_download_cel_files:
    input:
    output:
        raw_cel_files=directory(f'{results}/gdsc/array_data/raw'),
        array_metadata=f'{results}/gdsc/array_data/raw/E-MTAB-3610.sdrf.txt'
    threads: 1
    resources:
        mem_mb=1024
    conda:
        '../envs/gdsc_arrayexpress.yaml'
    script:
        '../scripts/gdsc_download_raw_files.R'


rule gdsc_normalize_arrays:
    input:
        cel_files= rules.gdsc_download_cel_files.output.raw_cel_files
    output:
        normalized_arrays=f'{results}/gdsc/array_data/normalized_arrays.rds'
    log:
        'logs/gdsc_normalize_arrays.log'
    threads: 1
    resources:
        mem_mb=20480
    conda:
        '../envs/gdsc_array_normalization.yaml'
    script:
        '../scripts/gdsc_normalize_arrays.R'


#TODO: output routes can be generalized
checkpoint gdsc_generate_compound_curves:
    input:
        dose_response_curves=datasets.loc['gdsc_response_curves', 'directory'],
        array_metadata=rules.gdsc_download_cel_files.output.array_metadata,
        cell_lines_annotation=rules.annotate_cell_lines.output.cell_lines_annotation,
    output:
        auc_models_candidates=directory(f'{results}/gdsc/auc_models_candidates'),
        compounds_lines_profiled=f'{results}/gdsc/compounds_lines_profiled.csv'
    log:
        'logs/gdsc_generate_compound_curves.log'
    threads: 1
    resources:
        mem_mb=16000
    conda:
        '../envs/common_file_manipulation.yaml'
    script:
        '../scripts/gdsc_generate_curves.R'


rule gdsc_compounds_diffexp:
    input:
        normalized_arrays  = rules.gdsc_normalize_arrays.output.normalized_arrays,
        compound_to_test = f'{results}/gdsc/auc_models_candidates/{{drug_id}}.csv'
    output:
        ebayes= f'{results}/gdsc/ebayes/{{drug_id}}_eBayes.rds'
    conda:
        '../envs/prism_limma.yaml'
    script:
        '../scripts/gdsc_generate_ebayes_model.R'