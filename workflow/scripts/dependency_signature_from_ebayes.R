suppressMessages(library('GSEABase'))
suppressMessages(library('limma'))
suppressMessages(library('dplyr'))


## SNAKEMAKE I/O  ##

eBayes_model       <- snakemake@input[['fitted_bayes']]
geneset_directory  <- snakemake@output[['bidirectional_geneset']]

gene_name          <- snakemake@wildcards[['gene']]

## Logging
log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")


 
generate_bidirectional_signature <- function(sig_name, deg_genes){
    
    ## Split the table between dep_associated_genes (positive logfold) 
    ## and non.dep associated genes (negative
    dependency_upregulated         <- deg_genes[deg_genes$logFC >= 0.5, 'ID']
    dependency_downregulated       <- deg_genes[deg_genes$logFC <= -0.5,'ID']
    
    if(length(dependency_upregulated) >= 15){
        
       candidate_genes  <- extract_top_genes(deg_genes, 'dependent_up')
       sensitivity_gset <- GSEABase::GeneSet(candidate_genes, geneIdType=SymbolIdentifier())
       setName(sensitivity_gset) <- paste(sig_name, 'UP', sep='_')
       
        GSEABase::toGmt(sensitivity_gset, con = paste(geneset_directory, '/', sig_name, '_dependency_UP', '.gmt', sep=''))
    }
    if(length(dependency_downregulated) >= 15){
        candidate_genes <- extract_top_genes(deg_genes, 'dependent_down')
        resistance_set  <- GSEABase::GeneSet(candidate_genes, geneIdType=SymbolIdentifier())
        setName(resistance_set) <- paste(sig_name, 'DOWN', sep='_')
        
        GSEABase::toGmt(resistance_set, con = paste(geneset_directory, '/', sig_name, '_dependency_DOWN', '.gmt', sep='')) 
    }
}

extract_top_genes <- function(deg_genes, mode='dependent_up'){
    
    if(mode=='dependent_up'){
        deg_genes <- deg_genes[deg_genes$logFC >= 0.5,]
        deg_genes <- head(deg_genes[order(deg_genes$logFC), 'ID'], n=500)
    }else{
        deg_genes <- deg_genes[deg_genes$logFC <= -0.5,]
        deg_genes <- head(deg_genes[order(-deg_genes$logFC), 'ID'], n=500)
    }
    
    return(deg_genes)
}


bayes_results <- readRDS(eBayes_model)

all_genes <- topTreat(bayes_results, coef = 'probability', 
                     number = Inf, adjust.method = 'fdr')

all_genes <- all_genes[all_genes$adj.P.Val <= 0.05, ]

all_genes$ID <- rownames(all_genes)

if(!dir.exists(file.path(geneset_directory))){
    dir.create(geneset_directory)}


if(nrow(all_genes) >= 15){

    sig_name <- paste(sep='_', gene_name, 'DepMap', '22Q2')
    generate_bidirectional_signature(sig_name, all_genes)
}

