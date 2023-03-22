import glob

# TODO: output routes can be generalized
checkpoint gdsc_rna_generate_compound_curves:
    input:
        dose_response_curves=datasets.loc["gdsc_response_curves", "directory"],
        count_matrix=rules.get_rnaseq_counts.output.raw_gene_counts,
        cell_lines_annotation=rules.annotate_cell_lines.output.cell_lines_annotation,
    output:
        auc_models_candidates=directory(f"{results}/gdsc_rna/auc_models_candidates"),
        compounds_lines_profiled=f"{results}/gdsc_rna/compounds_lines_profiled.csv",
    log:
        f"{LOGDIR}/gdsc_rna_generate_compound_curves/log.log",
    threads: get_resource("gdsc_generate_compound_curves", "threads"),
    resources:
        mem_mb=get_resource("gdsc_generate_compound_curves", "mem_mb"),
        walltime=get_resource("gdsc_generate_compound_curves", "walltime"),
    conda:
        "../envs/common_file_manipulation.yaml"
    script:
        "../scripts/gdsc_rna_generate_curves.R"


rule gdsc_rna_compounds_diffexpr:
    input:
        raw_gene_counts=rules.get_rnaseq_counts.output.raw_gene_counts,
        compound_to_test=f"{results}/gdsc_rna/auc_models_candidates/{{drug_id}}.csv",
    output:
        ebayes=f"{results}/gdsc_rna/ebayes/{{drug_id}}_eBayes.rds",
    log:
        f"{LOGDIR}/gdsc_rna_compounds_diffexpr/{{drug_id}}.log",
    threads: get_resource("gdsc_compounds_diffexp", "threads"),
    resources:
        mem_mb=get_resource("gdsc_compounds_diffexp", "mem_mb"),
        walltime=get_resource("gdsc_compounds_diffexp", "walltime"),
    conda:
        "../envs/prism_limma.yaml"
    script:
        "../scripts/gdsc_rna_generate_ebayes_model.R"


##TODO: These two rules could benefit from rule inheritance
rule gdsc_rna_geneset_from_ebayes_classic:
    input:
        fitted_bayes=rules.gdsc_compounds_diffexp.output.ebayes,
        treatment_info=datasets.loc["gdsc_response_curves", "directory"],
    output:
        bidirectional_geneset=directory(f"{results}/gdsc_rna/genesets/classic/{{drug_id}}"),
    params:
        signature_type="classic",
    threads: get_resource("gdsc_generate_geneset", "threads"),
    resources:
        mem_mb=get_resource("gdsc_generate_geneset", "mem_mb"),
        walltime=get_resource("gdsc_generate_geneset", "walltime"),
    log:
        f"{LOGDIR}/gdsc_rna/genesets/classic/{{drug_id}}.log",
    conda:
        "../envs/generate_genesets.yaml"
    script:
        "../scripts/gdsc_rna_signature_from_ebayes.R"


rule gdsc_rna_geneset_from_ebayes_fold:
    input:
        fitted_bayes=rules.gdsc_compounds_diffexp.output.ebayes,
        treatment_info=datasets.loc["gdsc_response_curves", "directory"],
    output:
        bidirectional_geneset=directory(f"{results}/gdsc_rna/genesets/fold/{{drug_id}}"),
    params:
        signature_type="fold",
    threads: get_resource("gdsc_generate_geneset", "threads"),
    resources:
        mem_mb=get_resource("gdsc_generate_geneset", "mem_mb"),
        walltime=get_resource("gdsc_generate_geneset", "walltime"),
    log:
        f"{LOGDIR}/gdsc_rna/genesets/fold/{{drug_id}}.log",
    conda:
        "../envs/generate_genesets.yaml"
    script:
        "../scripts/gdsc_rna_signature_from_ebayes.R"