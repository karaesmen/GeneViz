summaryData <- function(data, my.gene, my.disease, my.genome, my.cohort, my.outcome, result=c("gene.based", "topsnp", "toppval")){
        p <- data[gene==my.gene & disease==my.disease & genome==my.genome & cohort==my.cohort & outcome==my.outcome,]
        p <- unique(data.frame(p[,list(gene, geneBasedPvalue, topSNP, topSNPpVal)]))
        if(result=="gene.based"){
                paste0("Gene based p-value: ", round(p$geneBasedPvalue, 10))
        }
        else if(result=="topsnp"){
                paste0("Top SNP: ", p$topSNP)
                
        } else{
                paste0("Top SNP p-value: ", round(p$topSNPpVal, 10))
                
        }
}
