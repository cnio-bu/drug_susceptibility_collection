library('GSEABase')
library('limma')
library('dplyr')


## SNAKEMAKE I/O  ##

eBayes_model       <- snakemake@input[['fitted_bayes']]
drug_info          <- snakemake@input[['treatment_info']]

geneset_directory <- snakemake@output[['bidirectional_geneset']]

compound_id       <- snakemake@wildcards[['broad_id']]

generate_bidirectional_signature <- function(sig_name, deg_genes){
    
    ## Split the table between S_genes (negative logfold) 
    ## and R genes (positive)
    sensitivity_genes <- deg_genes[deg_genes$logFC <= -1, 'ID']
    resistance_genes  <- deg_genes[deg_genes$logFC >= 1,  'ID']
    
    if(length(sensitivity_genes) >= 15){
        
       candidate_genes <- extract_top_genes(deg_genes, 'sensitivity')
       sensitivity_gset <- GSEABase::GeneSet(candidate_genes, geneIdType=SymbolIdentifier())
       setName(sensitivity_gset) <- paste(sig_name, 'UP', sep='_')
       
        GSEABase::toGmt(sensitivity_gset, con = paste(geneset_directory, '/', sig_name, '_UP', '.gmt', sep=''))
    }
    if(length(resistance_genes) >= 15){
        candidate_genes <- extract_top_genes(deg_genes, 'resistance')
        resistance_set <- GSEABase::GeneSet(candidate_genes, geneIdType=SymbolIdentifier())
        setName(resistance_set) <- paste(sig_name, 'DN', sep='_')
        
        GSEABase::toGmt(resistance_set, con = paste(geneset_directory, '/', sig_name, '_DN', '.gmt', sep='')) 
    }
}

extract_top_genes <- function(deg_genes, mode='sensitivity'){
    
    if(mode=='sensitivity'){
        deg_genes <- deg_genes[deg_genes$logFC <= -1,]
        deg_genes <- head(deg_genes[order(deg_genes$logFC), 'ID'], n=500)
    }else{
        deg_genes <- deg_genes[deg_genes$logFC >= 1,]
        deg_genes <- head(deg_genes[order(-deg_genes$logFC), 'ID'], n=500)
    }
    
    return(deg_genes)
}


drug_info <- data.table::fread(drug_info) %>%
    as.data.frame()

drug_info <- drug_info[!duplicated(drug_info$broad_id),]
drug_info <- drug_info[drug_info$screen_id == 'HTS002',]

bayes_results <- readRDS(eBayes_model)

all_genes <- topTable(bayes_results, coef = 'auc', 
                     number = Inf, adjust.method = 'fdr',
                     p.value=0.05)

if(!dir.exists(file.path(geneset_directory))){
    dir.create(geneset_directory)}


if(nrow(all_genes) >= 15){

    common_name <- drug_info[drug_info$broad_id == compound_id, 'name']
        
    if(length(common_name) == 0){
        common_name <- compound_id}

    brd_split   <- strsplit(x=compound_id, split='-', fixed=TRUE)[[1]][2]
    
    sig_name <- paste(sep='_', common_name,'PRISM',brd_split)

    generate_bidirectional_signature(sig_name, all_genes)

}

