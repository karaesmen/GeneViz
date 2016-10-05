OSU_PC <- "C:/Users/karaesmen.1/Box Sync/Sucheston-Campbell Lab/Candidate Gene/awesome/"
OSU_VB <- "/media/sf_Sucheston-Campbell_Lab/Candidate Gene/awesome/"
ezgi_MB <- "~/Box Sync/Sucheston-Campbell Lab/Candidate Gene/awesome/"
setwd(ezgi_MB)

library(data.table)
library(Gviz)
library(manhattanly)

gene_loc <- read.table("./gene_locations.txt", stringsAsFactors = F, header = T)


## data.table of FULL.table -- Candidate Gene (CG) data
CG.data.raw <- fread("./FULL_table.txt",  header="auto", sep="auto")
CG.data.raw[, c("Pvalue.log10") := -log10(CG.data.raw[,"Pvalue", with = F])]
save(CG.data.raw, file="Full_table.RData")


##### Explore P-values ##### 
sig_SNPs = CG.data.raw[Pvalue.log10 > 7,] 
head(CG.data.raw[gene == "ABCB1",])
nrow(CG.data.raw)
anyNA(CG.data.raw)

p[rsID == "rs9516551" & impute == "typed" & disease == "noALL" & genome == "R"  & outcome == "TRM" ,]

test = CG.data.raw[gene == "ABCB1" & rsID == "rs28364280",]
anyNA(test)
test2 = dcast(test, gene + rsID + chr + BP + impute + disease + genome + outcome
              ~ cohort, value.var="Pvalue.log10")
anyNA(test2)
nas = which(is.na(test2$M))
test2[nas,]
###########################

CG.data <- dcast(CG.data.raw, gene + rsID + chr + BP + impute + disease + genome + outcome
            ~ cohort, value.var="Pvalue.log10")


# anyNA(CG.data)
# max(CG.data$c1)
# nas=which(is.na(CG.data$c1))
# CG.data[nas[1:10], ]


#setorder(CG.data, gene, BP)
CG.data = CG.data[order(gene, BP),]

## Generate the GRanges object
CG.data <- as(CG.data, "data.frame")

#############################
### exclude the wrong chrs ###

# genes to exclude for now 
genes2excl <- gene_loc[gene_loc$seqnames == "chr1", "gene_symbol"]

CG.data = CG.data[!c(CG.data$gene %in% genes2excl & !c(CG.data$chr == "chr1")), ]
#############################

#### explore CG.data####
head(CG.data)
tail(CG.data)

cutoff = -log10(5e-5)
CG.data[M >= cutoff ,]#| CG.data$c1 >= cutoff | CG.data$c2 >= cutoff), ]
head(CG.data[!is.na(CG.data),])

table(CG.data$gene)

## Prep CG.data for GRanges trans
CG.data$end <- CG.data$BP+1
CG.data$chr <- as.factor(CG.data$chr)

gr <- makeGRangesFromDataFrame(CG.data, ignore.strand = T, start.field = "BP", end.field = "end", 
                               starts.in.df.are.0based = T, keep.extra.columns = T) # is it 0based??

names(mcols(gr)) <- c(names(mcols(gr))[1:6], "Meta", "Cohort 1", "Cohort 2")

###################
###### Input ######

my.gene = "CTH"
my.disease = "ALLonly"
my.genome = "S"
my.cohort = c("Meta", "Cohort 1", "Cohort 2")
my.cohort.col = c('#66c2a5','#fc8d62','#8da0cb')
my.cohort.label = c("Cohort 1", "Cohort 2", "Meta")
my.outcome = c("DD", "OS", "PFS", "TRM")
my.ylim.upper = 8
gr.viz <- gr[c(gr$gene == my.gene & gr$disease == my.disease & gr$genome == my.genome), ]

my.gene.range = gene_loc[gene_loc$gene_symbol == my.gene, c("start_position", "end_position")]
my.chr <- unique(as.character(seqnames(gr.viz)))
### Gviz

getOption("Gviz.scheme")
scheme <- getScheme()
scheme$DataTrack$type <- "p"
scheme$DataTrack$ylim = c(0,my.ylim.upper)
addScheme(scheme, "GeneViz-Shiny")
options(Gviz.scheme="GeneViz-Shiny")

#availableDisplayPars(dTrack.DD)

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
           groups = my.cohort.label, 
           col.symbol = my.cohort.col,
           transcriptAnnotation="symbol", stackHeight = 0.8, 
           collapseTranscripts = "meta",
           baseline=c(0:8, -log10(5e-5), -log10(5e-8)),
           col.baseline=c(rep('gray', 9), '#4daf4a','#e41a1c'), lwd.baseline = c(rep(1, 9), 3, 3),
           from= my.gene.range$start_position - 20 , to=my.gene.range$end_position + 20, 
           cex.title=1.2, cex.axis=0.8, frame=T)




# # setkey(xx, genome)
# # xx.D <- xx["D",]
# # xx.D[,'threshold' := 8]
# # # xx.D <- xx[,c("gene", "rsID", "genome") := NULL]
# # xx.D <- xx.D[order(gene, BP)]
# 
# ## tidyR
# x <- read_tsv("./FULL_table.txt")
# head(x)
# class(x)
# 
# x_sub <- filter(x, gene == my.gene, disease == my.disease, genome==my.genome)
# x_sub_sp <- spread(x_sub, key=cohort, value=Pvalue, fill=NA)
# head(x_sub_sp)
# 
# xx_sub$cohort == "c2"]
# table(x_sub$cohort)
# 
# library(dplyr)
# stocks <- data.frame(
#   time = as.Date('2009-01-01') + 0:9,
#   X = rnorm(10, 0, 1),
#   Y = rnorm(10, 0, 2),
#   Z = rnorm(10, 0, 4)
# )
# stocksm <- stocks %>% gather( price, stock, -time)
# stockss <- stockss %>% spread(stock, price)
# stocksm %>% spread(time, price)
# 
# head(stocks)
# head(stocksm)

