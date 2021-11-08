checkpoint ctrp_annotate_models:
    input:
        curves_data       =  datasets.loc['ctrp_response_curves', 'directory'],
        compount_meta     =  datasets.loc['ctrp_compound_meta', 'directory'],
        cell_meta         =  datasets.loc['ctrp_cell_meta', 'directory'],
        experiment_meta   =  datasets.loc['ctrp_experiment_meta', 'directory'],
        cell_lines_annotation =  rules.annotate_cell_lines.output.cell_lines_annotation,
        rna_count_matrix      =  rules.get_rnaseq_counts.output.raw_gene_counts
    output:
        auc_models_candidates=directory(f'{results}/ctrp/auc_models_candidates'),
        compounds_lines_profiled=f'{results}/ctrp/compounds_lines_profiled.csv'
    log:
        'logs/ctrp_annotate_models.log'
    threads: 1
    resources:
        mem_mb=4096
    conda:
        '../envs/common_file_manipulation.yaml'
    script:
        '../scripts/ctrp_generate_curves.R'


rule ctrp_generate_ebayes:
    input:
        raw_gene_counts  = rules.get_rnaseq_counts.output.raw_gene_counts,
        compound_to_test = f'{results}/ctrp/auc_models_candidates/{{broad_id}}.csv'
    output:
        ebayes= f'{results}/ctrp/ebayes/{{broad_id}}_eBayes.rds'
    conda:
        '../envs/prism_limma.yaml'
    script:
        '../scripts/ctrp_generate_ebayes_model.R'


##TODO: These two rules could benefit from rule inheritance
rule ctrp_geneset_from_ebayes_classic:
    input:
        fitted_bayes= rules.ctrp_generate_ebayes.output.ebayes,
        treatment_info = datasets.loc['ctrp_compound_meta', 'directory']
    output:
        bidirectional_geneset=directory(f'{results}/ctrp/genesets/classic/{{broad_id}}')
    params:
        signature_type='classic'
    conda:
        '../envs/generate_genesets.yaml'
    script:
        '../scripts/ctrp_signature_from_ebayes.R'


rule ctrp_geneset_from_ebayes_fold:
    input:
        fitted_bayes= rules.ctrp_generate_ebayes.output.ebayes,
        treatment_info = datasets.loc['ctrp_compound_meta', 'directory']
    output:
        bidirectional_geneset=directory(f'{results}/ctrp/genesets/fold/{{broad_id}}')
    params:
        signature_type='fold'
    conda:
        '../envs/generate_genesets.yaml'
    script:
        '../scripts/ctrp_signature_from_ebayes.R'