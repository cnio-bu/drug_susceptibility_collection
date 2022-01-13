import glob

rule gdsc_download_cel_files:
    input:
        HTTP.remote(
            "https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-3610/E-MTAB-3610.sdrf.txt",
            keep_local=True,
        ),
    output:
        raw_cel_files=directory(f"{results}/gdsc/array_data/raw"),
        array_metadata=f"{results}/gdsc/array_data/raw/E-MTAB-3610.sdrf.txt",
    threads: 1
    resources:
        mem_mb=1024,
    conda:
        "../envs/gdsc_arrayexpress.yaml"
    script:
        "../scripts/gdsc_download_raw_files.R"


rule gdsc_normalize_arrays:
    input:
        cel_files=rules.gdsc_download_cel_files.output.raw_cel_files,
    output:
        normalized_arrays=f"{results}/gdsc/array_data/normalized_arrays.rds",
    log:
        f"{LOGDIR}/gdsc_normalize_arrays/log.log",
    threads: get_resource("gdsc_normalize_arrays", "threads"),
    resources:
        mem=get_resource("gdsc_normalize_arrays", "mem"),
        walltime=get_resource("gdsc_normalize_arrays", "walltime"),
    conda:
        "../envs/gdsc_array_normalization.yaml"
    script:
        "../scripts/gdsc_normalize_arrays.R"


# TODO: output routes can be generalized
checkpoint gdsc_generate_compound_curves:
    input:
        dose_response_curves=datasets.loc["gdsc_response_curves", "directory"],
        array_metadata=rules.gdsc_download_cel_files.output.array_metadata,
        cell_lines_annotation=rules.annotate_cell_lines.output.cell_lines_annotation,
    output:
        auc_models_candidates=directory(f"{results}/gdsc/auc_models_candidates"),
        compounds_lines_profiled=f"{results}/gdsc/compounds_lines_profiled.csv",
    log:
        f"{LOGDIR}/gdsc_generate_compound_curves/log.log",
    threads: get_resource("gdsc_generate_compound_curves", "threads"),
    resources:
        mem_mb=get_resource("gdsc_generate_compound_curves", "mem"),
        walltime=get_resource("gdsc_generate_compound_curves", "walltime"),
    conda:
        "../envs/common_file_manipulation.yaml"
    script:
        "../scripts/gdsc_generate_curves.R"


rule gdsc_compounds_diffexp:
    input:
        normalized_arrays=rules.gdsc_normalize_arrays.output.normalized_arrays,
        compound_to_test=f"{results}/gdsc/auc_models_candidates/{{drug_id}}.csv",
    output:
        ebayes=f"{results}/gdsc/ebayes/{{drug_id}}_eBayes.rds",
    threads: get_resource("gdsc_compounds_diffexp", "threads"),
    resources:
        mem=get_resource("gdsc_compounds_diffexp", "mem"),
        walltime=get_resource("gdsc_compounds_diffexp", "walltime"),
    conda:
        "../envs/prism_limma.yaml"
    script:
        "../scripts/gdsc_generate_ebayes_model.R"


rule gdsc_probeID_to_hgnc:
    input:
        cel_files=ancient(rules.gdsc_download_cel_files.output.raw_cel_files),
    output:
        tsv=f"resources/probeID_to_hgnc.tsv",
    log:
        f"{LOGDIR}/gdsc/probeID_to_hgnc.log",
    threads: get_resource("default", "threads"),
    resources:
        mem=get_resource("default", "mem"),
        walltime=get_resource("default", "walltime"),
    conda:
        "../envs/gdsc_probeID_to_hgnc.yaml"
    script:
        "../scripts/probeID_hgnc_table.R"


##TODO: These two rules could benefit from rule inheritance
rule gdsc_geneset_from_ebayes_classic:
    input:
        fitted_bayes=rules.gdsc_compounds_diffexp.output.ebayes,
        treatment_info=datasets.loc["gdsc_response_curves", "directory"],
        hgnc=rules.gdsc_probeID_to_hgnc.output.tsv,
    output:
        bidirectional_geneset=directory(f"{results}/gdsc/genesets/classic/{{drug_id}}"),
    params:
        signature_type="classic",
    threads: get_resource("gdsc_generate_geneset", "threads"),
    resources:
        mem=get_resource("gdsc_generate_geneset", "mem"),
        walltime=get_resource("gdsc_generate_geneset", "walltime"),
    log:
        f"{LOGDIR}/gdsc/genesets/classic/{{drug_id}}.log",
    conda:
        "../envs/generate_genesets.yaml"
    script:
        "../scripts/gdsc_signature_from_ebayes.R"


rule gdsc_geneset_from_ebayes_fold:
    input:
        fitted_bayes=rules.gdsc_compounds_diffexp.output.ebayes,
        treatment_info=datasets.loc["gdsc_response_curves", "directory"],
        hgnc=rules.gdsc_probeID_to_hgnc.output.tsv,
    output:
        bidirectional_geneset=directory(f"{results}/gdsc/genesets/fold/{{drug_id}}"),
    params:
        signature_type="fold",
    threads: get_resource("gdsc_generate_geneset", "threads"),
    resources:
        mem=get_resource("gdsc_generate_geneset", "mem"),
        walltime=get_resource("gdsc_generate_geneset", "walltime"),
    log:
        f"{LOGDIR}/gdsc/genesets/fold/{{drug_id}}.log",
    conda:
        "../envs/generate_genesets.yaml"
    script:
        "../scripts/gdsc_signature_from_ebayes.R"
