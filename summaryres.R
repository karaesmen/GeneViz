summary.res <- function(my.gene, my.disease, my.genome, my.cohort, my.outcome){
        header <- shiny::tags$h1(paste0("VEGAS2 gene-based results for ",  my.disease))#, " (", my.cohort, "/", my.outcome, ")"))
        gene.symbol <- paste0("Gene Symbol: ", my.gene)
        full.name <- paste0("Gene Name: ", mapIds(org.Hs.eg.db, toupper(my.gene),"GENENAME", "SYMBOL"))
        #gene.based <- summaryData(my.gene, my.disease, my.genome, my.cohort, my.outcome, "gene.based")
        #topsnp <- summaryData(my.gene, my.disease, my.genome, my.cohort, my.outcome, "topsnp")
        #topsnp.pval <- summaryData(my.gene, my.disease, my.genome, my.cohort, my.outcome, "toppval")
        #HTML(paste(header, gene.symbol, full.name, gene.based, topsnp, topsnp.pval, sep='<br/>'))
        HTML(paste(header, gene.symbol, full.name, sep='<br/>'))
}
