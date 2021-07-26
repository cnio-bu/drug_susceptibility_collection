library('biomart')
library('GSEABase')
library('limma')
library('tidyverse')

## SNAKEMAKE I/O  ##
eBayes_model       <- snakemake@input[['fitted_bayes']]
drug_info          <- snakemake@input[['dose_response_curves']]

geneset_directory <- snakemake@output[['bidirectional_geneset']]

compound_id       <- snakemake@wildcards[['drug_id']]

## Load annot.
ensembl_grch37 <- useMart(biomart="ENSEMBL_MART_ENSEMBL",
                          host="grch37.ensembl.org",
                          dataset="hsapiens_gene_ensembl")

human_grch38  <- useMart("ensembl", dataset="hsapiens_gene_ensembl")

## Function definitions ##

map_grch37_to_38 <- function(x){
    
    genesV2 = getLDS(attributes = c("ENSEMBL"),
                     filters = "ENSEMBL", values = x ,
                     mart = ensembl_grch37,
                     attributesL = c("ENSEMBL"),
                     martL = human_grch38, uniqueRows=FALSE)
    
    humanx <- unique(genesV2[])
    
    return(humanx)
}

annotateResults <- function(dataset){
    
    mouse_ensembl <- rownames(dataset)
    dataset       <- as.data.frame(dataset)
    
    dataset$mouse_ensembl <- mouse_ensembl
    
    dataset$mouse_symbol <- mapIds(org.Mm.eg.db,
                                   keys=rownames(dataset),
                                   column="SYMBOL",
                                   keytype="ENSEMBL",
                                   multiVals="first")
    
    mouse_symbol <- unique(dataset$mouse_symbol)
    human_genes <- convertMouseGeneList(mouse_symbol)
    
    dataset <- merge(x= dataset, y = human_genes, by.x = 'mouse_symbol', by.y = 'MGI.symbol', all.x = TRUE)
    names(dataset)[names(dataset) == "HGNC.symbol"] <- "human_symbol"
    
    return(dataset)
}

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

drug_info <-  readxl::read_xlsx(drug_info) %>%
                    as.data.frame()

drug_info <- drug_info[!duplicated(drug_info$DRUG_ID),]

bayes_results <- readRDS(eBayes_model)

all_genes <- topTable(bayes_results, coef = 'auc', 
                     number = Inf, adjust.method = 'fdr',
                     p.value=0.05)

if(!dir.exists(file.path(geneset_directory))){
    dir.create(geneset_directory)}


if(nrow(all_genes) >= 15){

    #TADO: annotate probes here

    common_name <- drug_info[drug_info$DRUG_ID == compound_id, 'DRUG_NAME']
        
    if(length(common_name) == 0){
        common_name <- compound_id}

    brd_split   <- strsplit(x=compound_id, split='-', fixed=TRUE)[[1]][2]
    
    sig_name <- paste(sep='_', common_name,'GDSC',brd_split)

    generate_bidirectional_signature(sig_name, all_genes)
}