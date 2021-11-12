import pandas as pd


def main():

    ###  SNAKEMAKE I/O ###
    raw_celligner = snakemake.input["celligner_lines_data"]
    lines_info = snakemake.input["sample_info"]
    where_to_save = snakemake.output["cell_lines_annotation"]

    raw_celligner = pd.read_csv(
        raw_celligner,
        sep=",",
        usecols=["sampleID", "sampleID_CCLE_Name", "undifferentiated_cluster"],
    )

    ## Get rid of non ccle data
    is_line = raw_celligner["sampleID"].str.startswith("ACH-")

    undifferentiated_lines = raw_celligner.loc[
        (is_line & raw_celligner["undifferentiated_cluster"]), "sampleID"
    ]

    del raw_celligner

    ## Annotate lines data with data from the undifferentiated clusters
    lines_info = pd.read_csv(lines_info, sep=",")
    lines_info["is_undifferentiated"] = lines_info["DepMap_ID"].isin(
        undifferentiated_lines
    )

    # rename the lineages from the undif. lines to "undifferentiated"
    lines_info.loc[lines_info["is_undifferentiated"], "lineage"] = "undifferentiated"

    lines_info.to_csv(where_to_save, sep=",", index=False)


if __name__ == "__main__":
    main()
