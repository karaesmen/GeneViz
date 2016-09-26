#BiocInstaller::biocLite("TxDb.Hsapiens.UCSC.hg19.knownGene")
#BiocInstaller::biocLite("org.Hs.eg.db")
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)


summaryData <- function(my.gene, my.disease, my.genome, my.cohort, my.outcome, result=c("gene.based", "topsnp", "toppval")){
        p <- x[gene==my.gene & disease==my.disease & genome==my.genome & cohort==my.cohort & outcome==my.outcome,]
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



