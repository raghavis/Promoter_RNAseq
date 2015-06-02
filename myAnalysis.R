#this is where you'll build your analysis script throughout the workshop
library(Rsubread)
library(limma)
library(edgeR)
library(ShortRead)
options(digits=2)
#library(Biostrisource("http://bioconductor.org/biocLite.R")
biocLite()source("http://bioconductor.org/biocLite.R")
biocLite()source("http://bioconductor.org/biocLite.R")
biocLite()ngs)

#read in your study design file
targets <- read.delim("retinaDegen_studyDesign.txt", sep="\t")
targets
myGroups <- factor(targets$description)
#create some more human-readable labels for your samples using the info in this file
sampleLabels <- paste(targets$name, targets$description, targets$rep, sep=".")

#set-up your experimental design
design <- model.matrix(~0+myGroups)
colnames(design) <- levels(myGroups)
design

# ##############################################################################################################################
# #you can check read quality using shortRead package
# #but I usually find it is better to do this on the sequencer or Illumina's BaseSpace website
# ##############################################################################################################################
# myFastq <- targets$fastq
# #collecting statistics over the files
# qaSummary <- qa("B6-WT-untreat-rep2_S2_mergedLanes_read1.fastq", type="fastq")
# #create and view a report
# browseURL(report(qaSummary))

##############################################################################################################################
#build index from your reference genome (expect this to take about 20 min on 8G RAM for mouse genome)
#you must have already downloaded the fasta file for your genome of interest and have it in the working directory
#this only needs to be done once, then index can be reused for future alignments
##############################################################################################################################
buildindex(basename="mouse",reference="Mus_musculus.GRCm38.dna.primary_assembly.fa")

##############################################################################################################################
#align your reads (in the fastq files) to your indexed reference genome that you created above
#expect this to take about 45min for a single fastq file containing 25 million reads
#the output from this is a .bam file for each of your original fastq files
##############################################################################################################################
reads1 <- targets$fastq[9]
reads2 <- targets$fastq[21] 
align(index="mouse", readfile1=reads1, readfile2=reads2, input_format="gzFASTQ",output_format="BAM",
      output_file="alignmentResultsPE_sample9.BAM", tieBreakHamming=TRUE,unique=TRUE,indels=5, nthreads=8)

##############################################################################################################################
#use the 'featureCounts' function to summarize read counts to genomic features (exons, genes, etc)
#will take about 1-2min per .bam file.
#for total transcriptome data summarized to mouse/human .gtf, expect about 50-60% of reads summarize to genes (rest is non-coding)
##############################################################################################################################
#read in text file with bam file names
MyBAM <- read.delim("BAMfiles_names.txt", header=T)
MyBAM <- as.character(MyBAM[,1])
#summarize aligned reads to genomic features (i.e. exons)
fc <- featureCounts(files=MyBAM, annot.ext="Mus_musculus.GRCm38.79.gtf", isGTFAnnotationFile=TRUE, GTF.featureType = "exon",
                    GTF.attrType="gene_id", useMetaFeatures=TRUE, isPairedEnd=TRUE, requireBothEndsMapped=TRUE, strandSpecific=2, nthreads=8)
#use the 'DGEList' function from EdgeR to make a 'digital gene expression list' object
DGEList <- DGEList(counts=fc$counts, genes=fc$annotation)
load("DGEList")
#retrieve all your gene/transcript identifiers from this DGEList object
myEnsemblIDs <- DGEList$genes$GeneID
#dim(fc$counts)
#tail(fc$counts)

##############################################################################################################################
#Normalize unfiltered data using 'voom' function in Limma package
#This will normalize based on the mean-variance relationship
#will also generate the log2 of counts per million based on the size of each library (also a form of normalization)
##############################################################################################################################
normData.unfiltered <- voom(DGEList, design, plot=TRUE)
exprs.unfiltered <- normData.unfiltered$E
#note that because you're now working with Log2 CPM, much of your data will be negative number (log2 of number smaller than 1 is negative) 
head(exprs.unfiltered)
dim(exprs.unfiltered)

#if you need RPKM for your unfiltered, they can generated as follows
#Although RPKM are commonly used, not really necessary since you don't care to compare two different genes within a sample
rpkm.unfiltered <- rpkm(DGEList, DGEList$genes$Length)
tail(rpkm.unfiltered)

##############################################################################################################################
#Filtering your dataset
#Only keep in the analysis those genes which had >10 reads per million mapped reads in at least two libraries.
##############################################################################################################################
cpm.matrix.filtered <- rowSums(cpm(DGEList) > 10) >= 2
DGEList.filtered <- DGEList[cpm.matrix.filtered,]
dim(DGEList.filtered)
#Use Voom again to normalize this filtered dataset
normData <- voom(DGEList.filtered, design, plot=TRUE)
exprs.filtered <- normData$E
head(exprs.filtered)

##############################################################################################################################
#annotate your normalized data using the organism-specific database package
##############################################################################################################################
library(org.Cf.eg.db)
library(AnnotationDbi)
#If we want to know what kinds of data are retriveable via the 'select' command, look at the columns of the annotation database
columns(org.Cf.eg.db)
#If we want to know what kinds of fields we could potentially use as keys to query the database, use the 'keytypes' command
keytypes(org.Cf.eg.db)
#transform you identifiers to entrezIDs
myAnnot.unfiltered <- select(org.Cf.eg.db, keys=rownames(exprs.unfiltered), keytype="ENSEMBL", columns=c("ENTREZID", "SYMBOL", "GENENAME"))
myAnnot.filtered <- select(org.Cf.eg.db, keys=rownames(exprs.filtered), keytype="ENSEMBL", columns=c("ENTREZID", "SYMBOL", "GENENAME"))
resultTable.unfiltered <- merge(myAnnot.unfiltered, exprs.unfiltered, by.x="ENSEMBL", by.y=0)
resultTable.filtered <- merge(myAnnot.filtered, exprs.filtered, by.x="ENSEMBL", by.y=0)
head(resultTable.unfiltered)
colnames(resultTable.unfiltered)
#add more appropriate sample names as column headers
colnames(resultTable.unfiltered) <- c("ensembl", "entrezID", "symbol", "name", sampleLabels)
colnames(resultTable.filtered) <- c("ensembl", "entrezID", "symbol", "name", sampleLabels)
#now write these annotated datasets out
write.table(resultTable.unfiltered, "normalizedUnfiltered.txt", sep="\t", quote=FALSE)
write.table(resultTable.filtered, "normalizedFiltered.txt", sep="\t", quote=FALSE)
head(resultTable.unfiltered)

#goal of this script is to using multivariate statisical approaches to explore the structure of your data
#begin by taking a look at the text expression matrix you created at the end of the last class
head(exprs.filtered)


###############################################################################################
#carry out hierarchical clustering on filtered data
###############################################################################################
#make some sample labels 
sampleLabels <- paste(targets$name, targets$description, targets$rep, sep=".")
distance <- dist(t(exprs.filtered),method="maximum")
clusters <- hclust(distance, method = "complete") 
plot(clusters, label = sampleLabels, hang = -1)

###############################################################################################
#Principal component analysis of the filtered data matrix
###############################################################################################
pca.res <- prcomp(t(exprs.filtered), scale.=F, retx=T)
ls(pca.res)
summary(pca.res) # Prints variance summary for all principal components.
head(pca.res$rotation) #$rotation shows you how much each GENE influenced each PC (callled 'eigenvalues', or loadings)
head(pca.res$x) #$x shows you how much each SAMPLE influenced each PC (called 'scores')
plot(pca.res, las=1)
pc.var<-pca.res$sdev^2 #sdev^2 gives you the eigenvalues
pc.per<-round(pc.var/sum(pc.var)*100, 1)
pc.per

#make some graphs to visualize your PCA result
library(ggplot2)
#lets first plot any two PCs against each other
#turn your scores for each gene into a data frame
data.frame <- as.data.frame(pca.res$x)
ggplot(data.frame, aes(x=PC1, y=PC2, colour=factor(myGroups))) +
  geom_point(size=5) +
  theme(legend.position="right")

#create a 'small multiples' chart to look at impact of each variable on each pricipal component
library(reshape2)
melted <- cbind(myGroups, melt(pca.res$x[,1:9]))
#look at your 'melted' data
ggplot(melted) +
  geom_bar(aes(x=Var1, y=value, fill=myGroups), stat="identity") +
  facet_wrap(~Var2)

#this script walks thorough some basic data wrangling for organizing expression data spreadsheets and ends with
#how to create publication quality graphics from transcriptomic data generated (regardless of platform used)
#to start this script you need a file with all your expression data and some non-redundant identifiers as row names (usually gene symbols)
#you also need a study design file

#for creating graphics, we'll use the ggplot2 and ggvis packages which employ a 'grammar of graphics' approach
#load the packages
library(ggplot2)
library(dplyr)
library(ggvis)

#read in your data from a text file that contains genes (or probe IDs) as rows, and samples as columns. 
myData <- read.delim("normalizedFiltered.txt", header=TRUE, row.names=1)
myData
#If you prefer, this data could be converted to a 'local' dataframe using the dplyr "tbl_df" function
#this simply makes reading the table a bit easier
#myData.local <- tbl_df(myData)
#myData.local

#column headers a bit cumbersome, so we'll change these to something more human-readable
geneSymbols <- myData[,1]
colnames(myData)

#use the dplyr 'mutate' command to get averages and fold changes for all your replicates
myData <- mutate(myData,
                 normal.AVG = (CEACID.normal.1 + CEACIV.normal.2 + CEACMI.normal.3)/3,
                 RCD1.AVG = (X2149.RCD1.1 + X2150.RCD1.2 + X2151.RCD1.3)/3,
                 XLPRA2.AVG = (Z478.XLPRA2.1 + Z479.XLPRA2.2 + Z480.XLPRA2.3)/3,
                 my.LogFC.RCD1 = (RCD1.AVG - normal.AVG),
                 my.LogFC.XLPRA2 = (XLPRA2.AVG - normal.AVG),
                 geneSymbols = geneSymbols)

#use dplyr "arrange" and "select" functions to sort by LogFC column of interest (arrange) 
#and then display only the columns of interest (select) to see the most differentially expressed genes
myData.sort <- myData %>%
  arrange(desc(my.LogFC.RCD1)) %>%
  select(geneSymbols, my.LogFC.RCD1, symbol, name)
head(myData.sort)

#use dplyr "filter" and "select" functions to pick out genes of interest (filter)
#and again display only columns of interest (select)
myData.filter <- myData %>%
  filter(geneSymbols=="DES") %>%
  select(geneSymbols, my.LogFC.RCD1, my.LogFC.XLPRA2)
head(myData.filter)

#make simple line or bar graph
#first reorder and clean up the data you want to graph
library(dplyr)
row.names(myData.filter) <- myData.filter[,1]
myData.filter <- select(myData.filter, -geneSymbols)
myData.filter.transpose <- as.data.frame(t(myData.filter))
treatments <- row.names(myData.filter.transpose)
myData.filter.transpose <- mutate(myData.filter.transpose, treatments, 
                                  genotype=c("WT", "het","KO"))
myData.filter.transpose

source("http://bioconductor.org/biocLite.R")
biocLite("dplyr")
source("http://bioconductor.org/biocLite.R")
biocLite("dplyr")
source("http://bioconductor.org/biocLite.R")
biocLite("ggvis")


#the goal of this script is to identify differentially expressed genes (DEG)
#you should already know what pairwise comparisons are most important to you

#I prefer to use a table of data with a non-redundant set of gene identifiers as input for this script
#this keeps my heatmaps and lists of DEGS as simple of as possible
#if you are coming to this script with array data, you've already reduced your list in this way
#RNAseq data is a bit more complicated (many transcripts, but the same gene), 
#but you achieve the same level of reduction by filtering based on annotation data
#################################################################
#additional filtering based on annotation
#if you have array data, you can skip this section
#################################################################
#our initial filtering of RNAseq data was based on solely cpm, which only removes relatively lowly expressed genes
#now that you have annotation info, you can further filter to reduce to single line of data per gene symbol
#first, check out how many unqiue genes are represented in your data
dupFiltered <- unique(resultTable.filtered$symbol)
#use collapseRows function from WGCNA package to collapse your dataset
library(WGCNA)
#pull your rownames and unique identifiers
myIDs <- rownames(resultTable.filtered)
#retrieve your gene symbols from your data
colnames(resultTable.filtered)
mySymbols <- resultTable.filtered[,3]
#remove all annotation columns so you're left with only numeric data
resultTable <- resultTable.filtered[,-1:-4]
myCollapsed <- collapseRows(resultTable, mySymbols, myIDs, method = "MaxMean")
myCollapsed <- myCollapsed$datETcollapsed
#now that the matrix is collapsed to give non-redudant list of genes,
#you could set the symbols to be the row names and move on

#################################################################################################################
#if you have no biological replicates, you will not be able to leverage statistical tools to identify DE genes
#Instead, you will ONLY rely on fold changes
#if you DO have replicates, skip this section and proceed to the next part 
#################################################################################################################
#use the dplyr 'filter' command to capture all the genes that are up/down regulated x-fold in n conditions
#in this case, 'myData' is a dataframe that you generated with Log2 expression and annotation
myData.filter <- myData %>%
  filter((abs(Ecdysone.vs.PBS_18hr_gut) >= 1) | (abs(Ecdysone.vs.PBS_5hr_gut) >= 1)) %>%
  select(geneID, Ecdysone.vs.PBS_5hr_carcass, Ecdysone.vs.PBS_18hr_carcass, Ecdysone.vs.PBS_5hr_gut, Ecdysone.vs.PBS_18hr_gut)
head(myData.filter)


###############################################################################################
# use Limma to find differentially expressed genes between two or more conditions
###############################################################################################
# fit the linear model to your filtered expression data
library(limma)
fit <- lmFit(myCollapsed, design)
#add annotation into your linear model fit
#don't really need to do this if you have RNAseq data
library(annotate)
fit$genes$Symbol <- getSYMBOL(probeList, "lumiMouseAll.db")
fit$genes$Entrez <- getEG(probeList, "lumiMouseAll.db")

# set up a contrast matrix based on the pairwise comparisons of interest
contrast.matrix <- makeContrasts(RCD1 = RCD1 - normal, XLPRA2 = XLPRA2 - normal, levels=design)

# check each contrast matrix
contrast.matrix

# extract the linear model fit for the contrast matrix that you just defined above
fits <- contrasts.fit(fit, contrast.matrix)
#get bayesian stats for your linear model fit
ebFit <- eBayes(fits)


###############################################################################################
# use the topTable and decideTests functions to see the differentially expressed genes
###############################################################################################

# use topTable function to take a look at the hits
myTopHits <- topTable(ebFit, adjust ="BH", coef=3, number=20, sort.by="logFC")
myTopHits

# use the 'decideTests' function to show Venn diagram for all diffexp genes for up to three comparisons
results <- decideTests(ebFit, method="global", adjust.method="BH", p.value=0.05, lfc=0.59)
#stats <- write.fit(ebFit)
vennDiagram(results, include="up") #all pairwise comparisons on a B6 background


# take a look at what the results of decideTests looks like
results

# now pull out probeIDs from selected regions of the Venn diagram.  In this case, I want all genes in the venn.
diffProbes <- which(results[,1] !=0 | results[,2] !=0 | results[,3] !=0)
diffSymbols <- fit$genes$Symbol[results[,1] !=0 | results[,2] !=0 | results[,3] !=0]
diffEntrez <- fit$genes$Entrez[results[,1] !=0 | results[,2] !=0 | results[,3] !=0]

#before pulling out expression data for differentially expressed genes, convert matrix to eset with annotation
library(Biobase)
myEset.ALL <- new("ExpressionSet", exprs = filtered.matrix)
annotation(myEset.ALL) <- "lumiMouseAll.db"

# retrieve expression data for the probes from above
diffData <- filtered.eset[results[,1] !=0 | results[,2] !=0 | results[,3] !=0]


#pull the expression data back out of the eset object
diffData <- exprs(diffData)

#combine probeIDs, gene symbols and expression data for differentially expressed genes into one file
write.table(cbind(diffProbes, diffSymbols, diffEntrez, diffData),"DiffGenes.xls", sep="\t", quote=FALSE)

# take a look at each expression matrix
dim(diffData)


