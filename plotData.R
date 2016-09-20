plotData <- function(gene, disease, genome, cohort, outcome){
        
        gr.viz <- gr[c(gr$gene == gene & gr$disease == disease & gr$genome == genome), ]
        
        gene_loc <- read.table("./data/gene_locations.txt", stringsAsFactors = F, header = T)
        my.gene.range = gene_loc[gene_loc$gene_symbol == gene, c("start_position", "end_position")]
        my.chr <- unique(as.character(seqnames(gr.viz)))
        
        getOption("Gviz.scheme")
        scheme <- getScheme()
        scheme$DataTrack$type <- "p"
        scheme$DataTrack$cex <- 0.5
        scheme$DataTrack$ylim = c(0,9)
        scheme$DataTrack$cex.legend = 1
        addScheme(scheme, "GeneViz-Shiny")
        options(Gviz.scheme="GeneViz-Shiny")
        
        bgr.col <-  c('#a6cee3','#1f78b4','#b2df8a','#33a02c')
        
        dTrack.DD <- DataTrack( gr.viz[gr.viz$outcome == "DD", cohort], 
                                name="DD (-log10)",  background.title=bgr.col[1], col.frame=bgr.col[1]) 
        
        dTrack.OS <- DataTrack( gr.viz[gr.viz$outcome == "OS", cohort], 
                                name="OS (-log10)", background.title=bgr.col[2]) 
        
        dTrack.PFS <- DataTrack( gr.viz[gr.viz$outcome == "PFS", cohort], 
                                 name="PFS (-log10)", background.title=bgr.col[3]) 
        
        dTrack.TRM <- DataTrack( gr.viz[gr.viz$outcome == "TRM", cohort], 
                                 name="TRM (-log10)", background.title=bgr.col[4]) 
        
        
        gtrack <- GenomeAxisTrack(background.title="DarkGray")
        itrack <- IdeogramTrack(genome= "hg19", chromosome =  my.chr, background.title="DarkGray",
                                rotation.title=0,  showBandId=T
                                #,name="Chr location", showTitle =F
        )
        
        
        biomTrack <- BiomartGeneRegionTrack(genome = "hg19",
                                            name = "ENSEMBL\n hg19", #symbol = my.gene, 
                                            chromosome=my.chr, start=my.gene.range$start_position, end=my.gene.range$end_position,
                                            background.panel = "#FFFEDB",
                                            background.title="DarkGray", min.height=10, arrowHeadMaxWidth=20, lwd=2, fill="salmon")
        
        
        myTracks <- list(itrack, gtrack, biomTrack, dTrack.DD, dTrack.OS, dTrack.PFS, dTrack.TRM)
        names(myTracks) <- c("itrack", "gtrack", "biomTrack", "DD", "OS", "PFS", "TRM")
        
        cohort.cols <- c('#66c2a5','#fc8d62','#8da0cb')
        
        
        plotTracks(myTracks[c("itrack", "gtrack", "biomTrack", outcome)], 
                   type=c( "p"), legend=T, cex=0.9, col.line = "gray", col= NULL,
                   groups = cohort, col.symbol = cohort.cols,
                   transcriptAnnotation="symbol", stackHeight = 0.7, 
                   fontsize=16, collapseTranscripts = "meta",
                   baseline=c(0:8, -log10(5e-5), -log10(5e-8)),
                   col.baseline=c(rep('gray', 9), '#4daf4a','#e41a1c'), lwd.baseline = c(rep(1, 9), 3, 3),
                   from= my.gene.range$start_position - 20 , to=my.gene.range$end_position + 20, 
                   legend.cex=0.3)
        
}