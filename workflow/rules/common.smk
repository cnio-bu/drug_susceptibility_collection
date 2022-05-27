rule annotate_cell_lines:
    input:
        sample_info=datasets.loc["ccle_sample_info", "directory"],
        celligner_lines_data=datasets.loc["celligner_data", "directory"],
    output:
        cell_lines_annotation=f"{results}/annotation/cell_line_annotation.csv",
    threads: get_resource("annotate_cell_lines", "threads"),
    resources:
        mem_mb=get_resource("annotate_cell_lines", "mem_mb"),
        walltime=get_resource("annotate_cell_lines", "walltime"),
    conda:
        "../envs/common_file_manipulation.yaml"
    script:
        "../scripts/common_annotate_lines.py"


rule get_rnaseq_counts:
    input:
        raw_expected_counts=datasets.loc["raw_ccle_reads", "directory"],
    output:
        raw_gene_counts=f"{results}/prism/raw_ccle_counts.rds",
    threads: get_resource("default", "threads"),
    resources:
        mem_mb=get_resource("default", "mem_mb"),
        walltime=get_resource("default", "walltime"),
    conda:
        "../envs/common_file_manipulation.yaml"
    script:
        "../scripts/prism_raw_counts_from_expected_counts.R"