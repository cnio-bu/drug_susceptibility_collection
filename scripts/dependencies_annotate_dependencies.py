import pandas as pd
import os

def main():

    ###  SNAKEMAKE I/O ###
    raw_crispr_data    = snakemake.input['crispr_gene_dependency_chronos']
    lines_info         = snakemake.input['sample_info']
    ccle_raw_reads     = snakemake.input['raw_expected_counts']
    where_to_save      = snakemake.output['model_candidates']
    

    # Prepare a dict. to convert DepMap IDs to stripped cell line name
    lines  = pd.read_csv(lines_info,
                         sep=',',
                         usecols=['DepMap_ID', 'stripped_cell_line_name', 'lineage'],
                         index_col=0
                         )
    
    human_readable = lines['stripped_cell_line_name'].to_dict()
    lines_lineage  = lines['lineage'].to_dict()

    raw_crispr_data = pd.read_csv(raw_crispr_data, sep=',', index_col=0)

    # Drop CRISPR profiled lines not present in RNA-Seq assay
    rna_seq_lines         = pd.read_csv(ccle_raw_reads, usecols=[0])
    rna_seq_lines.columns = ['profiled_lines']

    raw_crispr_data = raw_crispr_data[raw_crispr_data.index.isin(rna_seq_lines['profiled_lines'])]
    
    del rna_seq_lines

    # Transpose for better performance. Genes are rows now
    raw_crispr_data = raw_crispr_data.T

    # Get rid of genes with NaN values, usually due to defective processing/contamination
    filtered_crispr_data = raw_crispr_data.dropna(axis=0)

      # Drop engineered cell lines, If present
    is_engineered = lines.loc[lines['lineage'].str.contains('engineered'), 'stripped_cell_line_name']
    filtered_crispr_data = filtered_crispr_data.loc[:,~filtered_crispr_data.columns.isin(is_engineered)]

    # Get index as columns to avoid R calling it 'X'
    filtered_crispr_data.reset_index(drop=False, inplace=True)
    filtered_crispr_data.rename({'index':'Gene'}, axis=1, inplace=True)

    # CCLE nomenclature is HUGO (ENTREZ). Keep HUGO ONLY
    filtered_crispr_data['variance'] = filtered_crispr_data.var(axis=1)
    filtered_crispr_data['Gene']     = [x.split(' ')[0] for x in filtered_crispr_data['Gene'].values.tolist()]

    ## Resolve duplicates keeping most variable gene
    filtered_crispr_data_uniques = filtered_crispr_data.sort_values('variance', ascending=False).drop_duplicates('Gene', keep='first').sort_index()	
    filtered_crispr_data_uniques.drop('variance', axis=1, inplace=True)

    # melt the dataset 
    melted_crispr = filtered_crispr_data_uniques.melt(id_vars=['Gene'], var_name='cell_line', value_name='probability')

    # Filter out genes with less than 5 cell lines dependant on them
    melted_crispr            = melted_crispr.groupby('Gene').filter(lambda x: len(x[x['probability'] >= 0.95]) >= 5)

    # annotate the lineage and stripped name
    melted_crispr['lineage']            = melted_crispr['cell_line'].map(lines_lineage)
    melted_crispr['stripped_line_name'] = melted_crispr['cell_line'].map(human_readable)

    # TODO: we could add lineages before saving the models in the future
    # for now, It's pancancer
    melted_by_gene = melted_crispr.groupby('Gene')

    for gene_model in melted_by_gene.groups.items():

        model_data = melted_crispr.loc[gene_model[1]]
        gene_name  = model_data['Gene'].unique()[0]
        
        if not os.path.exists(where_to_save):
            os.mkdir(where_to_save)

        model_data.to_csv(f"{where_to_save}/{gene_name}.csv",sep=',', index=False, header=True)

    

if __name__ == "__main__":
    main()