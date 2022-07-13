rule drugbank_synonyms:
    input: 
        f"{DRUGBANKDIR}/drugBank.RData",
    output:  
        f"{results}/drugbank_synonyms.tsv",
    log:
        f"{LOGDIR}/drugbank_synonyms.log",
    params:
        blacklist = config["parameters"]["drugbank"]["blacklist"],
    threads: get_resource("drugbank_synonyms", "threads"),
    resources:
        mem_mb=get_resource("drugbank_synonyms", "mem_mb"),
        walltime=get_resource("drugbank_synonyms", "walltime"),
    conda:
        "../envs/beyondcell_signatures.yaml"
    script:
        "../scripts/drugbank.R"