genedata <- function(data,genes, my.disease, my.genome, my.cohort, my.outcome){
        data[gene==genes & disease==my.disease & genome==my.genome & cohort==my.cohort & outcome==my.outcome,-15,with=F]
}

