rule ctrp_build_db:
    input:
        compound_data=glob.glob(f"{results}/ctrp/auc_models_candidates/*.csv"),
        lines_compounds=rules.ctrp_annotate_models.output.compounds_lines_profiled,
        compound_meta=datasets.loc["ctrp_compound_meta", "directory"],
    output:
        csv_db=f"{results}/ctrp/drug_data.csv",
        rdata_db=f"{results}/ctrp/drug_data.rdata",
    log:
        f"{LOGDIR}/ctrp_build_db/log.txt",
    threads: get_resource("annotate_cell_lines", "threads"),
    resources:
        mem_mb=get_resource("annotate_cell_lines", "mem_mb"),
        walltime=get_resource("annotate_cell_lines", "walltime"),
    conda:
        "../envs/common_file_manipulation.yaml"
    script:
        "../scripts/ctrp_generate_drug_db.R"


rule gdsc_build_db:
    input:
        signatures_data=glob.glob(f"{results}/gdsc/auc_models_candidates/*.csv"),
        lines_compounds=rules.gdsc_generate_compound_curves.output.compounds_lines_profiled,
        cell_line_annotation = rules.annotate_cell_lines.output.cell_lines_annotation,
        compound_meta=datasets.loc["gdsc_compound_meta", "directory"],
    output:
        csv_db=f"{results}/gdsc/drug_data.csv",
        rdata_db=f"{results}/gdsc/drug_data.rdata",
    log:
        f"{LOGDIR}/gdsc_build_db/log.txt",
    threads: get_resource("annotate_cell_lines", "threads"),
    resources:
        mem_mb=get_resource("annotate_cell_lines", "mem_mb"),
        walltime=get_resource("annotate_cell_lines", "walltime"),
    conda:
        "../envs/common_file_manipulation.yaml"
    script:
        "../scripts/gdsc_generate_drug_db.R"


rule gdsc_rna_build_db:
    input:
        signatures_data=glob.glob(f"{results}/gdsc_rna/auc_models_candidates/*.csv"),
        lines_compounds=rules.gdsc_generate_compound_curves.output.compounds_lines_profiled,
        cell_line_annotation = rules.annotate_cell_lines.output.cell_lines_annotation,
        compound_meta=datasets.loc["gdsc_compound_meta", "directory"],
    output:
        csv_db=f"{results}/gdsc_rna/drug_data.csv",
        rdata_db=f"{results}/gdsc_rna/drug_data.rdata",
    log:
        f"{LOGDIR}/gdsc_build_db/log.txt",
    threads: get_resource("annotate_cell_lines", "threads"),
    resources:
        mem_mb=get_resource("annotate_cell_lines", "mem_mb"),
        walltime=get_resource("annotate_cell_lines", "walltime"),
    conda:
        "../envs/common_file_manipulation.yaml"
    script:
        "../scripts/gdsc_generate_drug_db.R"


rule prism_build_db:
    input:
        compound_data=glob.glob(f"{results}/prism/auc_models_candidates/*.csv"),
        lines_compounds=rules.prism_annotate_models.output.compounds_lines_profiled,
    output:
        csv_db=f"{results}/prism/drug_data.csv",
        rdata_db=f"{results}/prism/drug_data.rdata",
    log:
        f"{LOGDIR}/prism_build_db/log.txt",
    threads: get_resource("annotate_cell_lines", "threads"),
    resources:
        mem_mb=get_resource("annotate_cell_lines", "mem_mb"),
        walltime=get_resource("annotate_cell_lines", "walltime"),
    conda:
        "../envs/common_file_manipulation.yaml"
    script:
        "../scripts/prism_generate_drug_db.R"


rule build_drug_database:
    input:
        ctrp_db = rules.ctrp_build_db.output.csv_db,
        gdsc_db = rules.gdsc_build_db.output.csv_db,
        prism_db = rules.prism_build_db.output.csv_db,
        cell_line_annotation = rules.annotate_cell_lines.output.cell_lines_annotation,
    output:
        full_database=f"{results}/drug_database.csv",
        full_database_rd=f"{results}/drug_database.rdata",
        drug_summary=f"{results}/drug_summary.csv",
        drug_summary_rd=f"{results}/drug_summary.rdata",
    threads: get_resource("default", "threads"),
    resources:
        mem_mb=get_resource("default", "mem_mb"),
        walltime=get_resource("default", "walltime"),
    conda:
        "../envs/common_file_manipulation.yaml"
    script:
        "../scripts/common_build_drug_db.R"