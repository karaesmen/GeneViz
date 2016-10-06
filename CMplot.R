# circular manhattan plot

#install.packages("CMplot")
library(CMplot)
library(shiny)
library(data.table)
library(Gviz)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
setwd("~/Desktop/awesome/")

x <- fread("./data/FULL_table.txt",  header="auto", sep="auto")
#setkey(x, gene)
#genes <- data.frame(unique(x)[,"gene",with=F])
x[, c("Pvalue.log10") := -log10(x[,"Pvalue", with = F])]
xx <- dcast(x, gene + rsID + chr + BP + impute + disease + genome + outcome + geneBasedPvalue + topSNP + topSNPpVal + BPfromtopSNP 
            ~ cohort, value.var="Pvalue")
setorder(xx, gene, BP)


## making a circular manhattan plot using CMplot
# Pmap = df with 4 columns: 1. rsid, 2. chr, 3. chr location, 4. pvalues of each trait
library(plyr)
my.disease <- "AMLonly"
my.cohort <- "c2"
my.genome <- "D"
cm.plot <- function(my.disease, my.cohort, my.genome){
        cols <- c("rsID", "chr", "BP", "Pvalue", "outcome", "impute", "genome")
        df <- x[disease == my.disease & cohort==my.cohort & genome == my.genome,cols,with=F]
        df <- ddply(df, .(outcome), function(x) {x$row <- 1:nrow(x); x})
        dt <- dcast(df, rsID + chr + BP + impute + genome  ~ outcome, value.var="Pvalue")
        dt$impute <- NULL
        dt$row <- NULL
        dt$genome <- NULL
        chr <- strsplit(dt$chr, "[^[:digit:]]")
        dt$chr <- na.omit(as.numeric(unlist(chr)))
        colors <- c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99')
        CMplot(dt, col = colors, plot.type="c", multracks=TRUE, threshold=5e-08*max(dim(dt)),
               threshold.col=c("red", "red", "red", "red" ), cex=0.2)
}



?CMplot

