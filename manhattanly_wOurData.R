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

## Don't apply dcast!

#setorder(CG.data, gene, BP)
CG.data.raw = CG.data.raw[order(gene, BP),]

## Generate the GRanges object
CG.data.df <- as(CG.data.raw, "data.frame")

#############################
### exclude the wrong chrs ###

# genes to exclude for now 
genes2excl <- gene_loc[gene_loc$seqnames == "chr1", "gene_symbol"]

CG.data.df = CG.data.df[!c(CG.data.df$gene %in% genes2excl & !c(CG.data$chr == "chr1")), ]
#############################

head(CG.data.df)

# change Pvalue to "P", chr to "CHR" couldn't find the argument to specify!
# might be required for some other columns as well...
colnames(CG.data.df) <- c("gene","rsID" ,"CHR", "BP", "impute", "P" , "cohort" , "disease",
                          "genome" , "outcome" , "geneBasedPvalue", "topSNP", "topSNPpVal",     
                          "BPfromtopSNP" ,"Pvalue.log10")
# CHR column need to be numeric!
CG.data.df$CHR <- as.numeric(substr(CG.data.df$CHR, 4, 100))

head(CG.data.df)

############################
### Exploring the data ####

CG.data.df[CG.data.df$Pvalue.log10 > 7.3,]
sig.SNPs <- CG.data.df[CG.data.df$Pvalue.log10 > 3, ]
sig.SNPs <- sig.SNPs[order(sig.SNPs$rsID, sig.SNPs$cohort, sig.SNPs$disease, sig.SNPs$outcome),]
sig.SNPs[sig.SNPs$rsID == "rs113939899",]

### ABCC4 from the talk
ABCC4 = CG.data.df[CG.data.df$gene == "ABCC4" & CG.data.df$Pvalue < 0.05 , ]
ABCC4 = ABCC4[order(ABCC4$Pvalue, ABCC4$disease, ABCC4$cohort, ABCC4$disease),]
head(ABCC4, 10)

###################
###### Input ######

my.gene = "CTH"
my.disease = "AMLonly"
my.genome = "S"
# my.cohort = c("Meta", "Cohort 1", "Cohort 2")
my.cohort = "c2"
# my.cohort.col = c('#66c2a5','#fc8d62','#8da0cb')
# my.cohort.label = c("Cohort 1", "Cohort 2", "Meta")
# my.outcome = c("DD", "OS", "PFS", "TRM")
my.outcome = "TRM"
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
            annotation2 = "topSNP",
            #annotation3 = "", 
            #annotation4 = "topSNPpVal"
            )

