rule prism_get_rnaseq_counts:
    input:
        raw_expected_counts=datasets.loc['raw_ccle_reads', 'directory']
    output:
        raw_gene_counts=f'{results}/prism/raw_ccle_counts.rds'
    conda:
        '../envs/common_file_manipulation.yaml'
    script:
        '../scripts/prism_raw_counts_from_expected_counts.R'


## TODO: This rule can be split in 2
checkpoint prism_annotate_models:
    input:
        response_curves       =  datasets.loc['prism_response_curves', 'directory'],
        cell_lines_annotation =  rules.annotate_cell_lines.output.cell_lines_annotation,
        count_matrix=rules.prism_get_rnaseq_counts.output.raw_gene_counts
    output:
        auc_models_candidates=directory(f'{results}/prism/auc_models_candidates'),
        compounds_lines_profiled=f'{results}/prism/compounds_lines_profiled.csv'
    conda:
        '../envs/common_file_manipulation.yaml'
    script:
        '../scripts/prism_generate_annotation.R'


rule prism_compounds_diffexpr:
    input:
        raw_gene_counts  = rules.prism_get_rnaseq_counts.output.raw_gene_counts,
        compound_to_test = f'{results}/prism/auc_models_candidates/{{broad_id}}.csv'
    output:
        ebayes= f'{results}/prism/ebayes/{{broad_id}}_eBayes.rds'
    conda:
        '../envs/prism_limma.yaml'
    script:
        '../scripts/prism_generate_ebayes_model.R'


##TODO: These two rules could benefit from rule inheritance
rule prism_geneset_from_ebayes_classic:
    input:
        fitted_bayes= rules.prism_compounds_diffexpr.output.ebayes,
        treatment_info = datasets.loc['prism_treatment_info', 'directory']
    output:
        bidirectional_geneset=directory(f'{results}/prism/genesets/classic/{{broad_id}}')
    params:
        signature_type='classic'
    conda:
        '../envs/generate_genesets.yaml'
    script:
        '../scripts/prism_signature_from_ebayes.R'

rule prism_geneset_from_ebayes_fold:
    input:
        fitted_bayes= rules.prism_compounds_diffexpr.output.ebayes,
        treatment_info = datasets.loc['prism_treatment_info', 'directory']
    output:
        bidirectional_geneset=directory(f'{results}/prism/genesets/fold/{{broad_id}}')
    params:
        signature_type='fold'
    conda:
        '../envs/generate_genesets.yaml'
    script:
        '../scripts/prism_signature_from_ebayes.R'