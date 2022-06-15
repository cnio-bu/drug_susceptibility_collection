import glob

checkpoint mets_generate_models:
    input:
        metastasis_data=datasets.loc["met_experiment_data", "directory"],
        cell_lines_annotation=rules.annotate_cell_lines.output.cell_lines_annotation,
        count_matrix=rules.get_rnaseq_counts.output.raw_gene_counts,
    output:
        met_model_candidates=directory(f"{results}/mets/met_model_candidates"),
        met_lines_profiled=f"{results}/mets/met_lines_profiled.csv",
    threads: get_resource("ctrp_annotate_models", "threads"),
    resources:
        mem_mb=get_resource("ctrp_annotate_models", "mem_mb"),
        walltime=get_resource("ctrp_annotate_models", "walltime"),
    conda:
        "../envs/common_file_manipulation.yaml"
    script:
        "../scripts/metastasis_generate_models.R"


rule mets_generate_ebayes:
    input:
        raw_gene_counts=rules.get_rnaseq_counts.output.raw_gene_counts,
        mets_to_test=f"{results}/mets/met_model_candidates/{{met_type}}.csv",
    output:
        ebayes=f"{results}/mets/ebayes/{{met_type}}_eBayes.rds",
    log:
        f"{LOGDIR}/mets_compounds_diffexpr/{{met_type}}.log",
    threads: get_resource("gdsc_compounds_diffexp", "threads"),
    resources:
        mem_mb=get_resource("gdsc_compounds_diffexp", "mem_mb"),
        walltime=get_resource("gdsc_compounds_diffexp", "walltime"),
    conda:
        "../envs/prism_limma.yaml"
    script:
        "../scripts/metastasis_generate_ebayes_model.R"


rule mets_geneset_from_ebayes:
    input:
        fitted_bayes=rules.mets_generate_ebayes.output.ebayes,
    output:
        bidirectional_geneset=directory(f"{results}/mets/genesets/{{met_type}}"),
    log:
        f"{LOGDIR}/mets_genesets/{{met_type}}.log",
    threads: get_resource("ctrp_generate_geneset", "threads"),
    resources:
        mem_mb=get_resource("ctrp_generate_geneset", "mem_mb"),
        walltime=get_resource("ctrp_generate_geneset", "walltime"),
    conda:
        "../envs/generate_genesets.yaml"
    script:
        "../scripts/metastasis_geneset_from_ebayes.R.R"
