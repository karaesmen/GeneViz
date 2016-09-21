library(shiny)
library(data.table)
library(Gviz)
setwd("~/Desktop/OSU_PHD/Candidate_gene/gene_viz/")

x <- fread("./data/FULL_table.txt",  header="auto", sep="auto")


setkey(x, gene)
genes <- data.frame(unique(x)[,"gene",with=F])
x[, c("Pvalue.log10") := -log10(x[,"Pvalue", with = F])]

xx <- dcast(x, gene + rsID + chr + BP + impute + disease + genome + outcome + geneBasedPvalue + topSNP + topSNPpVal + BPfromtopSNP  
            ~ cohort, value.var="Pvalue")

setorder(xx, gene, BP)

## Generate the GRanges object
test <- DataFrame(xx)
test[["end"]] <- test[["BP"]] +1


gr <- makeGRangesFromDataFrame(test, ignore.strand = T, seqnames.field = "chr", 
                               start.field = "BP", end.field = "end", starts.in.df.are.0based = T, keep.extra.columns = T) # is it 0based??


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


ui <- fluidPage(
        titlePanel("DISCoVERY-BMT Replication of Candidate Gene Studies"),
        sidebarPanel(
                selectInput(inputId = "gene", label = "Select gene:", choices=genes),
                selectInput(inputId = "disease", label = "Select disease:", choices=c("ALLonly", "AMLonly", "mixed", "noALL")),
                selectInput(inputId = "genome", label = "Select genome:", choices=c("D", "R", "S")),
                checkboxGroupInput(inputId="outcome", label="Select survival outcomes:", choices=c("DD", "PFS", "OS", "TRM")),
                checkboxGroupInput(inputId="cohort", label="Select cohort:", choices=c("c1", "c2", "M")),
                actionButton("plotbutton", "Show plot")
        ),
        mainPanel(
                tabsetPanel(
                        tabPanel("Plot", plotOutput("gviz", width="100%", height="1000px")),
                        tabPanel("Summary", htmlOutput("gene_name"))
                )
                
        )
)       

server <- function(input, output){
        output$gene_name <- renderUI({
                header <- shiny::tags$h1(paste0(input$gene, " gene-based results for ",  input$disease, " (", input$cohort, "/", input$outcome, ")")) 
                gene.symbol <- paste0("Gene Symbol: ", input$gene)
                full.name <- paste0("Gene Name: ", mapIds(org.Hs.eg.db, toupper(gene),"GENENAME", "SYMBOL"))
                gene.based <- summaryData(input$gene, input$disease, input$genome, input$cohort, input$outcome, "gene.based")
                topsnp <- summaryData(input$gene, input$disease, input$genome, input$cohort, input$outcome, "topsnp")
                topsnp.pval <- summaryData(input$gene, input$disease, input$genome, input$cohort, input$outcome, "toppval")
                HTML(paste(header, gene.symbol, full.name, gene.based, topsnp, topsnp.pval, sep='<br/>'))
        })
        
        output$gviz <- renderPlot({
                input$plotbutton
                isolate(plotData(input$gene, input$disease, input$genome, input$cohort, input$outcome))
        })
        
}


shinyApp(ui = ui, server = server)
