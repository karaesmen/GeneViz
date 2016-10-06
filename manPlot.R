## Interactive Manhattan Plot
#devtools::install_github("sahirbhatnagar/manhattanly", build_vignettes = TRUE)
library(manhattanly)
library(data.table)

manPlot <- function(CG.data, my.disease, my.genome, my.cohort, my.outcome){
        CG.data.raw = CG.data[order(gene, BP),]
        CG.data.df <- as(CG.data.raw, "data.frame")
        
        gene_loc <- read.table("./data/gene_locations.txt", stringsAsFactors = F, header = T)
        genes2excl <- gene_loc[gene_loc$seqnames == "chr1", "gene_symbol"]
        
        CG.data.df = CG.data.df[!c(CG.data.df$gene %in% genes2excl & !c(CG.data$chr == "chr1")), ]
        
        # change Pvalue to "P", chr to "CHR" couldn't find the argument to specify!
        # might be required for some other columns as well...
        colnames(CG.data.df) <- c("gene","rsID" ,"CHR", "BP", "impute", "P" , "cohort" , "disease",
                                  "genome" , "outcome" , "geneBasedPvalue", "topSNP", "topSNPpVal",     
                                  "BPfromtopSNP" ,"Pvalue.log10")
        # CHR column need to be numeric!
        CG.data.df$CHR <- as.numeric(substr(CG.data.df$CHR, 4, 100))
        
        chr.col =  c('#b3cde3','#8c96c6','#8856a7','#810f7c')
        suggestiveline_col = "#08519c"
        genomewideline_col = "#a50f15"
        
        #my.ylim.upper = 8
        CG.data.df.viz <- CG.data.df[c( #CG.data.df$gene == my.gene & 
                CG.data.df$disease == my.disease & 
                        CG.data.df$genome == my.genome & CG.data.df$cohort == my.cohort &
                        CG.data.df$outcome == my.outcome), c("gene", "rsID", 'CHR','BP', 'geneBasedPvalue', 'topSNP', 'P')]
        
        head(CG.data.df.viz[order(CG.data.df.viz$P, decreasing = T),])
        
        manhattanly(CG.data.df.viz, snp = "rsID", gene = "gene", col = chr.col,
                    suggestiveline_color = suggestiveline_col, suggestiveline_width = 2,
                    genomewideline_color = genomewideline_col, genomewideline_width = 2,
                    title= "Interactive Manhattan Plot", showgrid = T,
                    annotation1 = "geneBasedPvalue", 
                    annotation2 = "topSNP"
 
        )
        
}

