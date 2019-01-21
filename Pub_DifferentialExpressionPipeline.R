#------------------------------------------#
# R codes for RNA-Seq QC and Normalization on CDTA dataset #
#------------------------------------------#
#Adapted from the RNAseq workshop led by Dr. Jenny Drnevich at the University
#of Illinois at Urbana-Champaign's High Performance Computing Biocluster in Fall 2016

#Assigning workspace on computer
getwd()
setwd("~/OneDrive - University of Florida/Grad School/Pubs/MastersThesisPub/Data/")

#Setting up workspace
options(stringsAsFactors = FALSE)

#Install necessary packages
source("https://bioconductor.org/biocLite.R")
biocLite('sva')
biocLite('affycoretools')
biocLite('rgl')
biocLite('limma')
biocLite('edgeR')
biocLite('rtracklayer')
biocLite("impute")
biocLite("statmod")

#Relevant citation information, though incomplete
citation("edgeR")
citation("sva")
citation("limma")

#Load package libraries
library(sva)
library(WGCNA) 
library(affycoretools)
library(limma)
library(rgl)
library(edgeR)
library(rgl)
library(rtracklayer)

#Clear workspace
ls()
rm(list = ls())


#Set file input directory for quantified read files ending in 'abundance.tsv'
FileName <- dir(path = "~/OneDrive - University of Florida/Grad School/Pubs/MastersThesisPub/Data/",pattern = "abundance.tsv")
FileName

#Keeping track of filenames with our count data
write.csv(FileName, file = "CountFileNames.csv")

#Manually created .txt file with details on samples, sampling scheme, filenames, read totals, and mapped reads
targets <- readTargets("TargetsCDTA.txt")
targets

#Setup categorical variables to explore data
targets$group <- paste(targets$Treatment, targets$Rep, sep = ".")
targets

targets$GpF <- factor(targets$group, levels = unique(targets$group))
targets

targets$col <- as.numeric(targets$GpF)
targets

#Calculate read fates across samples
read.fate <- data.frame(Unaligned = targets$Total - targets$Mapped,
                        Mapped = targets$Mapped,
                        row.names = targets$Label) 
read.fate

read.fate <- read.fate / targets$Total
read.fate

rowSums(read.fate)

apply(read.fate, 2, summary)

#Construct DGE object containing all necessary read information across samples and explore the data
d <- readDGE(targets, labels = targets$Label, header = TRUE,
             columns = c(1,5))
d

dim(d)

#Plot your read densities by sample
x11()
plotDensities(log2(d$counts + 0.1), group = d$samples$GpF, col = 1:8)
?plotDensities

x11()
plot(colSums(d$counts), colSums(d$counts == 0), pch = 16, col = d$samples$col)

#Clustering algorithm to show expression similarity
hc.raw <- hclust(dist(t(log2(d$counts + 0.1))), method = "average")

#Now plot it:
x11(width = 12 , height = 6 )
plot(hc.raw, hang = -1, main = "CDTA Transcriptome, raw values", sub = "", xlab = "", cex = 0.9)

#To change the labels to anything else you can do:
x11(width = 12 , height = 6 )
plot(hc.raw, hang = -1, main = "CDTA Transcriptome, raw values", sub = "", xlab = "", cex = 0.9,
     labels = d$samples$group)

# First, check the screeplot to see how many of the PCs explain much variation.
# By definition, PC1 > PC2 > PC3, etc.
x11()  
plotPCA(log2(d$counts + 0.1), screeplot = T)
?plotPCA

# Second, plot PC1 vs. PC2, coloring the samples by group and adding
# on the sample names:
x11()
plotPCA(log2(d$counts + 0.1), groups = d$samples$col, groupnames = levels(d$samples$GpF),
        addtext = d$samples$Sample, main = "PCA on raw counts, CDTA" )

#3D plot
plotPCA(log2(d$counts + 0.1), pch = 16, col = d$samples$col, groupnames = levels(
  d$samples$GpF), plot3d = T, pcs = 1:3)

#Histogram of counts per transcript
x11()
hist(log2(rowSums(d$counts + 0.1)), 100, main = "total counts per transcript, log2 scale")

#Now we normalize by library size using TMM method
d <- calcNormFactors(d)
d$samples

#Plot normalization
x11()
barplot(d$samples$norm.factors, col = d$samples$col, main = "TMM Normalization Factors by Sample", xlab ="Sample",
        ylab = "Normalization Value",names.arg = d$samples$Label, las = 2)

#Calculating counts per million reads
cpm.values <- cpm(d$counts)
head(cpm.values)

above1cpm <- rowSums(cpm.values  >= 1)
length(above1cpm)
table(above1cpm)

x11()
hist(above1cpm, xlab = "Number of samples with > 1 CPM, CDTA", ylab = "Number of genes")

sum(above1cpm >= 1)

sum(above1cpm >= 1) / length(above1cpm)
temp <- log2(rowSums(d$counts + 0.1))
x11()
layout(matrix( 1:3 , 3, 1 ))
hist(temp, 1000, main = "all genes", xlim = c(0, 20))
hist(temp[above1cpm >= 1], 1000, main = "kept genes", xlim = c(0, 20))
hist(temp[above1cpm < 1], 1000, main = "filtered genes", xlim = c(0, 20))

#Keeping only those reads with >= 1 cpm
d.filt <- d[above1cpm >= 1, , keep.lib.sizes = FALSE]
filteredTags <- d.filt$counts[,0]
filteredTags
write.csv(filteredTags, "filteredTags.txt")

d.filt$samples$lib.size / d$samples$lib.size
nrow(d.filt) / nrow(d)

d.filt <- calcNormFactors(d.filt)
d.filt$samples

d.filt$samples$norm.factors / d$samples$norm.factors
x11()
par(mar = c(8,4.1,4.1,2.1))
barplot( d.filt$samples$lib.size / 1e6, main = "Actual Library Sizes, CDTA", 
         las = 2, names.arg = d$samples$Label, col = d$samples$col, ylab = "millions of reads", 
         cex.axis = 0.8 )

x11()
par(mar = c(8,4.1,4.1,2.1))
barplot( d.filt$samples$lib.size * d.filt$samples$norm.factors / 1e6, main = "Effective Library
         sizes, CDTA", las = 2, names.arg = d$samples$Label, col = d$samples$col, ylab = "millions of reads", 
         cex.axis = 0.8 )

x11()
plotMDS(d.filt, main = "MDS plot, CDTA", col = d.filt$samples$col, top = 1000 )


#### The edgeR vignette suggests using the cpm() function to get modified CPM values
#### like before, except on the log2 scale for clustering or heatmaps. Since you can't
#### take a log of 0, a small value needs to be added. The default is to add an average 
#### of 0.25 counts to each gene, proportional to library size:

logCPM <- cpm(d.filt, log = T)
class(logCPM)
dim(logCPM)
head(logCPM)
logCPM
d.filt

#### Examine the logCPM values; first, the distribution of values:

x11()
plotDensities( logCPM, group = d$samples$GpF, col = 1:4 )


#### Second, basic hierarchical clustering:

hc.edgeR <- hclust(dist(t(logCPM)), method = "average")

x11( 12, 6 )
plot(hc.edgeR, hang = -1, main = "CDTA, edgeR logCPM values", sub = "", xlab = 
       "", cex = 0.9)


#### Third, can do both PCA plots. First, check the screeplot to see how much variation 
#### each PC explains:

x11()   
plotPCA(logCPM,screeplot=T)


# Second, plot PC1 vs. PC2, coloring the samples by group and adding
# on the sample names, and also do the movable 3D plot:

x11()  
plotPCA(logCPM, pch = 16, col = rep(1:2,each=4), groupnames = c("C","E"), 
        addtext = d.filt$samples$Sample, main = "PCA edgeR logCPM, CDTA")

plotPCA(logCPM, pch = 16, col = rep(1:2,each=4), groupnames = c("C","E"), 
        plot3d = TRUE, pcs = 1:3)


design <- cbind(C=c(1,1,1,1,1,1,1,1),E=c(0,0,0,0,1,1,1,1))
design

d.filt$samples$group <- as.character(d.filt$samples$group)
d.filt$samples$group <- factor(c("C","C","C","C","E","E","E","E"))
d.filt$samples$group

design2 <- model.matrix(~d.filt$samples$group) #run after changing d.filt sample$groups to C and E
design2
colnames(design2) <- c("C","EvsC")
design2

rownames(design2) <- d.filt$sample$Label
design2

cont.matrix <- makeContrasts(CvsE = C-EvsC, levels=design2)
cont.matrix

##For voom and removing batch effects
x11()
d.filt$counts
v <- voom(d.filt$counts,design,plot=T)
summary(v)

#######################################
#Taking out batch effect using surrogate variable analysis, sva package

n.sv <- num.sv(v$E,design)
n.sv
head(v$E)
design <- design2[,1:2]
design
svobj <- sva(v$E,design,design[,1],n.sv=n.sv)
svobj
names(svobj)
colnames(svobj$sv) <- c("SV1","SV2","SV3") 
colnames(svobj$sv)

design2sv <- cbind(design,svobj$sv)
design2sv
x11()
design
design2
v2 <- voom(d.filt,design2sv,plot=T)
v2
no.sv2 <- removeBatchEffect(v2, design= design, covariates = svobj$sv,plot=T)
summary(no.sv2)

fit.sv2 <- eBayes(lmFit(no.sv2,design),trend=T)
x11()
plotSA(fit.sv2)
summary(fit.sv2)
plot(fit.sv2$coefficients)

summary(decideTests(fit.sv2,p.value=0.05,adjust.method="BH",lfc=1))
topEBayes <- topTable(fit.sv2,coef=2,lfc=1,number=305,confint=T,p.value=0.05,sort.by="p",adjust.method="BH")
summary(topEBayes)

x11()
hist(topEBayes)
hist(fit.sv2$p.value[,2], 1000)

x11()
plotPCA(no.sv2,groups=rep(1:2,each=4),groupnames = c("C","E"),addtext = 1:8)

topEBayes$FC <- (2^abs(topEBayes$logFC)) * sign(topEBayes$logFC)

x11()
plot(topEBayes$FC)
hist(topEBayes$FC)

write.csv(topEBayes, file = "topEBayesResults.csv")


#To construct heatmap
biocLite("Heatplus")
library(gplots)
summary(codes.voom <- decideTests(fit.sv2,p.value=0.05,lfc=1))
x11()
layout(matrix(1:3,3,1)) #this divides the graph into 3 parts
vennDiagram(codes.voom[,2], include = "both", main = "both UP and DOWN")
vennDiagram(codes.voom[,2], include = "up", main = "UP only")
vennDiagram(codes.voom[,2], include = "down", main = "DOWN only")

#codes2.CvsE
summary(codes.voom)
G.up <- rownames(codes.voom)[codes.voom[,2] != 0 &
                               codes.voom[,2] != -1]
G.up
write.csv(G.up, file = "CDTAUpRegulatedTranscripts.csv")

G.down <- rownames(codes.voom)[codes.voom[,2] != 0 &
                                 codes.voom[,2] != 1]
G.down
write.csv(G.down, file = "CDTADownRegulatedTranscripts.txt")

length(G.up)  #check on how many selected
length(G.down)

G.all <- c(G.up, G.down)
length(G.all)
G.all
heatdata <- logCPM[G.all, ]
heatdata
dim(heatdata)

heatdata.baseFC <- heatdata - rowMeans(heatdata[ , d.filt$samples$group=="C"])
range(heatdata.baseFC)

#### Check distribution of values to determine best color scale:

x11()
hist(heatdata.baseFC, 100)

#### Load in a function I got from the BioC mailing list to make a centered color
#### scale for a particular set of data:
biocLite("MKmisc")
library("MKmisc")
source("heatmapCol.R")

max(heatdata.baseFC)
min(heatdata.baseFC)
color.baseFC <- heatmapCol(data = heatdata.baseFC, lim = 6, 
                           col =colorRampPalette(c("blue","white","red"))(128))

heatdata.meanFC <- heatdata - rowMeans(heatdata)

x11()
hist(heatdata.meanFC, 100)

#smaller range and more symmetric. To be
#consistent, let's also use a color scale range of -3 to 3
min(heatdata.meanFC)
max(heatdata.meanFC)

color.meanFC <- heatmapCol(data = heatdata.meanFC, lim = 7, 
                           col =colorRampPalette(c("blue","white","red"))(128))

heatdata.scaled <- t(scale(t(heatdata)))
dim(heatdata.scaled)
x11()
hist(heatdata.scaled, 100)

i.names <- c("C1","C2","C3","C4","E1","E2","E3","E4")
i.names
colnames(heatdata.scaled)[i.names]
x11()
heatmap.2(heatdata.meanFC[ , i.names], col = color.meanFC, scale = "none", 
          labRow = F, trace = "none", density.info = "none", margins = c(5,2), 
          symbreaks= FALSE, key.xlab = "logFC to mean of all samples")

cluster.all <- hclust(dist(heatdata.baseFC))