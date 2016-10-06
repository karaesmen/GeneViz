library(data.table)
library(Gviz)


plotData <- function(data, my.gene, my.disease, my.genome, my.cohort, my.outcome){
        my.cohort.col = c('#66c2a5','#fc8d62','#8da0cb')
        my.ylim.upper = 8
        cg.df <- data.frame(data)
        
        #############################
        ### exclude the wrong chrs ###
        # genes to exclude for now 
        gene_loc <- read.table("./data/gene_locations.txt", stringsAsFactors = F, header = T)
        genes2excl <- gene_loc[gene_loc$seqnames == "chr1", "gene_symbol"]
        cg.df <- cg.df[!c(cg.df$gene %in% genes2excl & !c(cg.df$chr == "chr1")), ]
        
        #############################
        # need to have an end to our SNP location
        cg.df[["end"]] <- cg.df[["BP"]]+1
        # make chromosomes into factors
        cg.df[["chr"]] <- as.factor(cg.df$chr)
        # create GRanges object
        gr <- makeGRangesFromDataFrame(cg.df, ignore.strand = T, seqnames.field = "chr",
                                       start.field = "BP", end.field = "end", starts.in.df.are.0based = T, keep.extra.columns = T) # is it 0based??
        names(mcols(gr)) <- c(names(mcols(gr))[1:6], "Meta", "Cohort 1", "Cohort 2")
        
        ## inputs
        gr.viz <- gr[c(gr$gene == my.gene & gr$disease == my.disease & gr$genome == my.genome), ]
        my.gene.range = gene_loc[gene_loc$gene_symbol == my.gene, c("start_position", "end_position")]
        my.chr <- unique(as.character(seqnames(gr.viz)))
        
        ## Gviz Options
        getOption("Gviz.scheme")
        scheme <- getScheme()
        scheme$DataTrack$type <- "p"
        scheme$DataTrack$cex <- 0.5
        scheme$DataTrack$ylim = c(0,9)
        scheme$DataTrack$cex.legend = 1
        addScheme(scheme, "GeneViz-Shiny")
        options(Gviz.scheme="GeneViz-Shiny")
        
        bgr.col <-  c('#a6cee3','#1f78b4','#b2df8a','#33a02c')
        
        dTrack.DD <- DataTrack( gr.viz[gr.viz$outcome == "DD", my.cohort], 
                                name="DD (-log10)",  background.title=bgr.col[1], 
                                col.frame=bgr.col[1]) 
        
        dTrack.OS <- DataTrack( gr.viz[gr.viz$outcome == "OS", my.cohort], 
                                name="OS (-log10)", background.title=bgr.col[2],
                                col.frame=bgr.col[2])  
        
        dTrack.PFS <- DataTrack( gr.viz[gr.viz$outcome == "PFS", my.cohort], 
                                 name="PFS (-log10)", background.title=bgr.col[3],
                                 col.frame=bgr.col[3]) 
        
        dTrack.TRM <- DataTrack( gr.viz[gr.viz$outcome == "TRM", my.cohort], 
                                 name="TRM (-log10)", background.title=bgr.col[4],
                                 col.frame=bgr.col[4]) 
        
        gtrack <- GenomeAxisTrack(background.panel="transparent")
        itrack <- IdeogramTrack(genome= "hg19", chromosome =  my.chr, showBandId=T,
                                background.panel="transparent", size=1.8, fontsize=12)
        #rotation.title=0,  
        #,name="Chr location", showTitle =F
        biomTrack <- BiomartGeneRegionTrack(genome = "hg19",
                                            name = "ENSEMBL\n hg19", #symbol = my.gene, 
                                            chromosome=my.chr, start=my.gene.range$start_position, end=	my.gene.range$end_position,
                                            background.panel = "#FFFEDB",
                                            background.title="salmon", min.height=10, arrowHeadMaxWidth=20, lwd=2, fill="salmon",
                                            cex=1, size=2, fontsize=16)
        
        myTracks <- list(itrack, gtrack, biomTrack, dTrack.DD, dTrack.OS, dTrack.PFS, dTrack.TRM)
        names(myTracks) <- c("itrack", "gtrack", "biomTrack", "DD", "OS", "PFS", "TRM")
        
        plotTracks(myTracks[c("itrack", "gtrack", "biomTrack", my.outcome)], 
                   type=c( "p"), legend=T, cex=1.2, col.line = "gray", col= NULL,
                   groups=my.cohort, 
                   col.symbol = my.cohort.col,
                   transcriptAnnotation="symbol", stackHeight = 0.8, 
                   collapseTranscripts = "meta",
                   baseline=c(0:8, -log10(5e-5), -log10(5e-8)),
                   col.baseline=c(rep('gray', 9), '#4daf4a','#e41a1c'),
                   lwd.baseline=c(rep(1, 9), 3, 3),
                   from=my.gene.range$start_position - 20,
                   to=my.gene.range$end_position + 20, 
                   cex.title=1.2,
                   cex.axis=0.8,
                   frame=T)
        
}
