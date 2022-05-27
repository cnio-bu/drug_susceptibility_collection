library('GSEABase')
library('limma')
library('dplyr')


## SNAKEMAKE I/O  ##

eBayes_model       <- snakemake@input[['fitted_bayes']]
drug_info          <- snakemake@input[['treatment_info']]

geneset_directory <- snakemake@output[['bidirectional_geneset']]

compound_id       <- snakemake@wildcards[['broad_id']]

## SNAKEMAKE PARAMS ##
signature_type    <- tolower(as.character(snakemake@params[['signature_type']]))

## If signature_type == 'Classic',
## then we'll take the top/bottom 250 genes sorted by T.statistics
## If signature_type == 'fold', we'll perform FDR + log fold selection

generate_bidirectional_signature <- function(sig_name, deg_genes){
    
    ## Split the table between S_genes (negative logfold) 
    ## and R genes (positive)
    if(signature_type=='classic'){
        sensitivity_genes <- deg_genes[deg_genes$t <0, 'ID']
        resistance_genes  <- deg_genes[deg_genes$t >0, 'ID']
    }else{
        sensitivity_genes <- deg_genes[deg_genes$logFC <= -1, 'ID']
        resistance_genes  <- deg_genes[deg_genes$logFC >= 1, 'ID']
    }
    if(length(sensitivity_genes) >= 15){
       candidate_genes <- extract_top_genes(deg_genes, 'sensitivity', signature_type=signature_type)
       sensitivity_gset <- GSEABase::GeneSet(candidate_genes, geneIdType=SymbolIdentifier())
       setName(sensitivity_gset) <- paste(sig_name, 'UP', sep='_')
       
       GSEABase::toGmt(sensitivity_gset, con = paste(geneset_directory, '/', sig_name, '_UP', '.gmt', sep=''))
    }
    if(length(resistance_genes) >= 15){
        candidate_genes <- extract_top_genes(deg_genes, 'resistance', signature_type=signature_type)
        resistance_set <- GSEABase::GeneSet(candidate_genes, geneIdType=SymbolIdentifier())
        setName(resistance_set) <- paste(sig_name, 'DOWN', sep='_')
        
        GSEABase::toGmt(resistance_set, con = paste(geneset_directory, '/', sig_name, '_DOWN', '.gmt', sep='')) 
    }
}

#TODO: This function here is extremely ugly. Refactor it
extract_top_genes <- function(deg_genes, mode='sensitivity', signature_type='classic'){
    
    ## Extract and sort the top hits for the selected direction/mode
    if(mode=='sensitivity'){

        if(signature_type=='classic'){
            deg_genes <- head(deg_genes[order(deg_genes$t), 'ID'], n=500)
        }else{
            deg_genes <- deg_genes[deg_genes$logFC <= -1,]
            deg_genes <- head(deg_genes[order(deg_genes$logFC), 'ID'], n=250)
        }
    }else{
        if(signature_type=='classic'){
            deg_genes <- head(deg_genes[order(-deg_genes$t), 'ID'], n=500)
        }else{
        deg_genes <- deg_genes[deg_genes$logFC >= 1,]
        deg_genes <- head(deg_genes[order(-deg_genes$logFC), 'ID'], n=250)
        }
    }
    return(deg_genes)
}


drug_info <- data.table::fread(drug_info) %>%
    as.data.frame()

drug_info <- drug_info[!duplicated(drug_info$broad_cpd_id),]

bayes_results <- readRDS(eBayes_model)

if(signature_type == 'classic'){
    all_genes <- topTable(bayes_results, coef='area_under_curve', number = Inf, adjust.method='fdr')
}else{ 
    ## Set a threshold of the adjusted p-value if signature_type is not classic
    all_genes <- topTable(bayes_results, coef = 'area_under_curve', 
                     number = Inf, adjust.method = 'fdr',
                     p.value=0.05)
}
# Set the genes as a column
all_genes$ID <- rownames(all_genes)

if(!dir.exists(file.path(geneset_directory))){
    dir.create(geneset_directory)}

## Check if we have at least 15 genes left after FDR-cutoff. Does not matter
## If both directions are < 15 individually, we'll account for that later.
## Just make sure we are not evaluating empty results 

if(nrow(all_genes) >= 15){

    common_name <- drug_info[drug_info$broad_cpd_id == compound_id, 'cpd_name']
        
    if(length(common_name) == 0){
        common_name <- compound_id}

    brd_split   <- strsplit(x=compound_id, split='-', fixed=TRUE)[[1]][2]
    
    # This code deals with combinations of drugs, to keep names short
    if(grepl(pattern = "mol/mol", x = common_name, fixed = TRUE)){
     common_name <- strsplit(x = common_name, split = " ")[[1]][1]}

    sig_name <- paste(sep='_', common_name,'CTRP', brd_split)

    generate_bidirectional_signature(sig_name, all_genes)

}
