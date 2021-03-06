vegas2Summary <- function(x, my.gene, my.disease, my.genome){
        cols <- c("gene", "cohort", "outcome", "geneBasedPvalue", "topSNP", "topSNPpVal")
        p <- x[gene==my.gene & disease==my.disease & genome==my.genome, cols, with=FALSE]
        p <- data.frame(p)
        colnames(p) <- c("Gene Symbol", "Cohort", "Survival Outcome", "Gene Based P-value", "Top SNP", "Top SNP P-value")
        unique(p)
}