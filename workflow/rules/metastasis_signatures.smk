import glob

checkpoint met_generate_models:
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