suppressMessages(library("GSEABase"))
suppressMessages(library("limma"))
suppressMessages(library("dplyr"))


## SNAKEMAKE I/O  ##

eBayes_model       <- snakemake@input[["fitted_bayes"]]
geneset_directory  <- snakemake@output[["bidirectional_geneset"]]

met_name          <- snakemake@wildcards[["met_type"]]
model_type        <- snakemake@params[["model_type"]]

## Logging
log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")


 
generate_bidirectional_signature <- function(sig_name, deg_genes){
    
    ## Split the table between dep_associated_genes (positive logfold) 
    ## and non.dep associated genes (negative
    dependency_upregulated         <- deg_genes[deg_genes$logFC > 0, "ID"]
    dependency_downregulated       <- deg_genes[deg_genes$logFC < 0, "ID"]
    
    if(length(dependency_upregulated) >= 5){
        
       candidate_genes  <- extract_top_genes(deg_genes, "met_UP")
       sensitivity_gset <- GSEABase::GeneSet(candidate_genes, geneIdType=SymbolIdentifier())
       setName(sensitivity_gset) <- paste(sig_name, "UP", sep="_")
       
        GSEABase::toGmt(sensitivity_gset, con = paste(geneset_directory, "/", sig_name, "_", model_type, "_met_UP", ".gmt", sep=""))
    }
    if(length(dependency_downregulated) >= 5){
        candidate_genes <- extract_top_genes(deg_genes, "met_DOWN")
        resistance_set  <- GSEABase::GeneSet(candidate_genes, geneIdType=SymbolIdentifier())
        setName(resistance_set) <- paste(sig_name, "DOWN", sep="_")
        
        GSEABase::toGmt(resistance_set, con = paste(geneset_directory, "/", sig_name, "_", model_type, "_met_DOWN", ".gmt", sep="")) 
    }
}

extract_top_genes <- function(deg_genes, mode="met_UP"){
    
    if(mode=="met_UP"){
        deg_genes <- head(deg_genes[order(deg_genes$logFC, decreasing = TRUE), "ID"], n=250)
    }else{
        deg_genes <- head(deg_genes[order(deg_genes$logFC, decreasing = FALSE), "ID"], n=250)
    }
    
    return(deg_genes)
}


bayes_results <- readRDS(eBayes_model)

all_genes <- topTreat(bayes_results, coef = "mean", 
                     number = Inf, adjust.method = "fdr")

all_genes <- all_genes[all_genes$adj.P.Val <= 0.05, ]

all_genes$ID <- rownames(all_genes)

if(!dir.exists(file.path(geneset_directory))){
    dir.create(geneset_directory)}


if(nrow(all_genes) >= 5){

    sig_name <- paste(sep="_", met_name, model_type, "MetMap", "2022")
    generate_bidirectional_signature(sig_name, all_genes)
}

