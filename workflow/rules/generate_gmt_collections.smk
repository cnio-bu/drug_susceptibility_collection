rule get_gmt_collection:
    input:
        get_ctrp_genesets_classic,
        get_gdsc_genesets_classic,
        get_prism_genesets_classic,
        f"{results}/drug_database.rdata",
    output:
        f"{results}/drug_signatures_classic.gmt",
    threads: get_resource("get_gmt_collections", "threads"),
    resources:
        mem_mb=get_resource("get_gmt_collection", "mem_mb"),
        walltime=get_resource("get_gmt_collection", "walltime"),
    shell:
        "cat {results}/ctrp/genesets/classic/*/*.gmt {results}/gdsc/genesets/classic/*/*.gmt {results}/prism/genesets/classic/*/*.gmt > {output}"