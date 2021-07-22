checkpoint dependencies_annotate_crispr_data:
    input:
        crispr_gene_dependency_chronos=datasets.loc['crispr_gene_dependency_chronos', 'directory'],
        sample_info=rules.annotate_cell_lines.output.cell_lines_annotation,
        raw_expected_counts=datasets.loc['raw_ccle_reads', 'directory']
    output:
        model_candidates=directory(f'{results}/dependencies/model_candidates')
    threads: 1
    resources:
        mem=4096,
        walltime=10
    conda:
        '../envs/common_file_manipulation.yaml'
    script:
        '../scripts/dependencies_annotate_dependencies.py'


rule dependencies_generate_ebayes:
    input:
        raw_gene_counts    = rules.prism_get_rnaseq_counts.output.raw_gene_counts,
        dependency_to_test = f'{results}/dependencies/model_candidates/{{gene}}.csv'
    output:
        ebayes= f'{results}/dependencies/ebayes/{{gene}}_eBayes.rds'
    threads: 1
    resources:
        mem=8192,
        walltime=10
    conda:
        '../envs/prism_limma.yaml'
    script:
        '../scripts/dependencies_generate_ebayes_model.R'


rule dependencies_geneset_from_ebayes:
    input:
        fitted_bayes=rules.dependencies_generate_ebayes.output.ebayes
    output:
        bidirectional_geneset=directory(f'{results}/dependencies/genesets/{{gene}}')
    conda:
        '../envs/generate_genesets.yaml'
    script:
        '../scripts/dependency_signature_from_ebayes.R'