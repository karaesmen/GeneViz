library(shiny)
library(data.table)
library(Gviz)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
library(plotly)
library(manhattanly)


# read in full data table
x <- readRDS("./data/FULL_table.rds")


# append -log10(pvalue) column
x[, c("Pvalue.log10") := -log10(x[,"Pvalue", with = F])]

# move pvalues from rows to columns
x.dcast <- dcast(x, gene + rsID + chr + BP + impute + disease + genome + outcome ~ cohort, value.var="Pvalue.log10")

# get unique list of genes
setkey(x, gene)
genes <- data.frame(unique(x[,"gene",with=F]))

source("summaryres.R")
source("plotData.R")
source("vegas2Summary.R")
source("manPlot.R")
source("genedata.R")

ui <- shinyUI(pageWithSidebar(
        headerPanel("DISCOVeRY-BMT Replication of Candidate Gene Studies"),
        sidebarPanel(
                conditionalPanel(condition="input.conditionedPanels==1",
                                 selectInput(inputId = "gene",
                                             label = "Select gene:",
                                             choices=genes),
                                 selectInput(inputId = "disease",
                                             label = "Select disease:",
                                             choices=c("ALLonly", "AMLonly", "mixed", "noALL")),
                                 selectInput(inputId = "genome",
                                             label = "Select genome:",
                                             choices=c("D", "R", "S")),
                                 checkboxGroupInput(inputId="outcome",
                                                    label="Select survival outcomes:",
                                                    choices=c("DD", "PFS", "OS", "TRM")),
                                 checkboxGroupInput(inputId="cohort",
                                                    label="Select cohort:",
                                                    choices=c("Cohort 1", "Cohort 2", "Meta")),
                                 actionButton("plotbutton", "Show Plot")
                ),
                conditionalPanel(condition="input.conditionedPanels==2",
                                 selectInput(inputId = "gene2",
                                             label = "Select gene:",
                                             choices=genes),
                                 selectInput(inputId = "disease2",
                                             label = "Select disease:",
                                             choices=c("ALLonly", "AMLonly", "mixed", "noALL")),
                                 selectInput(inputId = "genome2",
                                             label = "Select genome:",
                                             choices=c("D", "R", "S")),
                                 # checkboxGroupInput(inputId="outcome2",
                                 #                    label="Select survival outcomes:",
                                 #                    choices=c("DD", "PFS", "OS", "TRM")),
                                 # checkboxGroupInput(inputId="cohort2",
                                 #                    label="Select cohort:",
                                 #                    choices=c("c1", "c2", "M")),
                                 actionButton("summarybutton", "Generate Summary")
                ),
                conditionalPanel(condition="input.conditionedPanels==3",
                                 selectInput(inputId = "disease3",
                                             label = "Select disease:",
                                             choices=c("ALLonly", "AMLonly", "mixed", "noALL")),
                                 selectInput(inputId = "genome3",
                                             label = "Select genome:",
                                             choices=c("D", "R", "S")),
                                 selectInput(inputId="outcome3",
                                                     label="Select survival outcomes:",
                                                     choices=c("DD", "PFS", "OS", "TRM")),
                                 selectInput(inputId="cohort3",
                                                     label="Select cohort:",
                                                     choices=c("c1", "c2", "M")),
                                 actionButton("manbutton", "Generate Manhattan Plot")
                ),
                conditionalPanel(condition="input.conditionedPanels==4",
                                 selectInput(inputId = "gene4",
                                             label = "Select gene:",
                                             choices=genes),
                                 selectInput(inputId = "disease4",
                                             label = "Select disease:",
                                             choices=c("ALLonly", "AMLonly", "mixed", "noALL")),
                                 selectInput(inputId = "genome4",
                                             label = "Select genome:",
                                             choices=c("D", "R", "S")),
                                 checkboxGroupInput(inputId="outcome4",
                                                    label="Select survival outcomes:",
                                                    choices=c("DD", "PFS", "OS", "TRM")),
                                 checkboxGroupInput(inputId="cohort4",
                                                    label="Select cohort:",
                                                    choices=c("c1", "c2", "M")),
                                 actionButton("genebutton", "Generate Raw Summary"))
                ),
        mainPanel(
                tabsetPanel(
                        id = "conditionedPanels",
                        tabPanel("Gene Viewer",
                                 uiOutput("ui_plot"),
                                 value=1),
                        tabPanel("VEGAS2 Summary",
                                 htmlOutput("gene_info"),
                                 dataTableOutput("table"),
                                 downloadButton('downloadvegas', 'Download'),
                                 value=2),
                        tabPanel("Manhattan Plot",
                                 plotlyOutput("manhattanly"),
                                 uiOutput("ui_man"),
                                 value=3),
                        tabPanel("Gene Info",
                                 dataTableOutput("gene_table"),
                                 downloadButton("genedownload", 'Download'),
                                 value=4)
                )
        )
))

server <- function(input, output){
        plot_height <- function(){
                ceiling((length(input$outcome)/2.5)*1000)
        }

        output$gviz <- renderPlot({
                input$plotbutton
                isolate(plotData(x.dcast, input$gene, input$disease, input$genome, input$cohort, input$outcome))
        })
        
        output$ui_plot <- renderUI({
                input$plotbutton
                isolate(plotOutput("gviz", height=plot_height(), width="100%"))
        })
        
        output$gene_info <- renderUI({
                input$summarybutton
                isolate(summary.res(input$gene2, input$disease2, input$genome2, input$cohort2, input$outcome2))
        })
        
        output$table <- renderDataTable({
                input$summarybutton
                isolate(vegas2Summary(x, input$gene2, input$disease2, input$genome2, c("c1", "c2", "M"), c("DD", "PFS", "OS", "TRM")))
        },
        options = list(bLengthChange=FALSE, sDom  = '<"top">lrt<"bottom">ip') #remove the entries dropdown box
        )
        
        output$downloadvegas <- downloadHandler(
                filename = function () { c(paste(input$gene2, input$disease2, input$genome2, sep ="-"), paste(".txt", sep=""))},
                content = function(file){write.table(vegas2Summary(input$gene2, input$disease2, input$genome2, c("c1", "c2", "M"), c("DD", "PFS", "OS", "TRM")), file, quote=F, sep="\t")}
        )

        output$manhattanly <- renderPlotly({
                input$manbutton
                isolate(manPlot(x, input$disease3, input$genome3, input$cohort3, input$outcome3))
        })

        output$man_info <- renderUI({
                input$manbutton
                shiny::tags$h1(paste0("Interactive Manhattan Plot for ", input$disease3, " (", input$genome3, "/", input$cohort3, ")"))
        })

        # output$ui_man <- renderUI({
        #         input$manbutton
        #         isolate(plotlyOutput("manhattanly"))
        # })
        
        output$gene_table <- renderDataTable({
                input$genebutton
                isolate(genedata(x, input$gene4, input$disease4, input$genome4, input$cohort4, input$outcome4))
        })
        
        output$genedownload <- downloadHandler(
                filename = function () { c(paste(input$gene4, input$disease4, input$genome4, "FULL", "GENE", sep ="-"), paste(".txt", sep=""))},
                content = function(file){write.table(genedata(x, input$gene4, input$disease4, input$genome4, input$cohort4, input$outcome4), file, quote=F, sep="\t")}
        )
        
}

shinyApp(ui = ui, server = server)
