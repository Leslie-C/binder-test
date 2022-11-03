#     Example DESeq2 workflow developed for Launch, Nov 4, 2022
#     by Leslie Coonrod, Assoc. Dir, Bioinformatics and Genomics Master's Program, 
#                        Knight Campus Graduate Internship Program, University of Oregon
#              https://internship.uoregon.edu/bioinformatics
#
#     See https://bioconductor.org/packages/release/bioc/html/DESeq2.html
#     or https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
#     for more information.
#
#                                            UOUOUOUOUOUOU                                 
#                               UOUOUOUOUOUOUOUOUOUOUOUOUOUOUOUOUOUOUOU                    
#                        UOUOUOUOUOUOUOUOUOUOUOUOUOUOUOUOUOUOUOUOUOUOUOUOUOUO              
#                    UOUOUOUOUOUOUOUOUOU                      UOUOUOUOUOUOUOUOUOU          
#                 UOUOUOUOUOUOUOUOUOU                             UOUOUOUOUOUOUOUOUO       
#               UOUOUOUOUOUOUOUOUO                                  UOUOUOUOUOUOUOUOUOU    
#             UOUOUOUOUOUOUOUOUOU                                     UOUOUOUOUOUOUOUOUO   
#            UOUOUOUOUOUOUOUOUOU                                      UOUOUOUOUOUOUOUOUOU  
#            UOUOUOUOUOUOUOUOUOU                                       UOUOUOUOUOUOUOUOUOU 
#           UOUOUOUOUOUOUOUOUOUO                                       UOUOUOUOUOUOUOUOUOUO
#           UOUOUOUOUOUOUOUOUOU                                        UOUOUOUOUOUOUOUOUOUO
#           UOUOUOUOUOUOUOUOUOU                                        UOUOUOUOUOUOUOUOUOUO
#          UOUOUOUOUOUOUOUOUOUO                                        UOUOUOUOUOUOUOUOUOUO
#          UOUOUOUOUOUOUOUOUOUO                                        UOUOUOUOUOUOUOUOUOUO
#          UOUOUOUOUOUOUOUOUOUO                                        UOUOUOUOUOUOUOUOUOUO
#           UOUOUOUOUOUOUOUOUOU                                        UOUOUOUOUOUOUOUOUOUO
#           UOUOUOUOUOUOUOUOUOU                                        UOUOUOUOUOUOUOUOUOUO
#           UOUOUOUOUOUOUOUOUOUO                                       UOUOUOUOUOUOUOUOUOU 
#            UOUOUOUOUOUOUOUOUOU                                       UOUOUOUOUOUOUOUOUOU 
#            UOUOUOUOUOUOUOUOUOU                                      UOUOUOUOUOUOUOUOUOU  
#              UOUOUOUOUOUOUOUOUO                                     UOUOUOUOUOUOUOUOUO   
#               UOUOUOUOUOUOUOUOUOU                                 UOUOUOUOUOUOUOUOUO     
#                  UOUOUOUOUOUOUOUOUO                             UOUOUOUOUOUOUOUOUO       
#                     UOUOUOUOUOUOUOUOUOU                     UOUOUOUOUOUOUOUOUOU          
#                         UOUOUOUOUOUOUOUOUOUOUOUOUOUOUOUOUOUOUOUOUOUOUOUOUOU              
#                               UOUOUOUOUOUOUOUOUOUOUOUOUOUOUOUOUOUOUO                     


#install necessary packages -> Not required if using Binder image
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("DESeq2")
# BiocManager::install("apeglm")
# install.packages("ashr")
# install.packages("ggplot2")
# install.packages("pheatmap")

library(DESeq2)

#Create table with sample metadata
directory <- paste(getwd(),"/data/", sep="")
sampleFiles <- grep("mouse",list.files(directory),value=TRUE)
sampleCondition1 <- sub("mouse","\\1",sampleFiles)
sampleCondition <- sub("\\d.*","\\1",sampleCondition1)
type <- sub("\\w+\\d.","\\1",sampleCondition1)
type <- sub("*.genecount","\\1",type)
sampleTable <- data.frame(sampleName=sub("*.\\w+_\\w+.genecount","\\1",sampleCondition1),
                          fileName=sampleFiles,
                          condition=sampleCondition)
print(sampleTable)

#Import data from HTSeqCount output
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable=sampleTable,
                                       directory=directory,
                                       design=~ condition)

#set factor levels
ddsHTSeq$condition <- factor(ddsHTSeq$condition,
                             levels = c("Saline","Hept","Furam"))
#Run DESeq
dds <- DESeq(ddsHTSeq)

#remove rows with 0 or 1 reads
dds <- dds[rowSums(counts(dds)) > 1, ]

#contrast saline v. treatments
resHept <- results(dds, contrast=c("condition","Hept","Saline"))
resFur <- results(dds, contrast=c("condition","Furam","Saline"))
resHeptOrdered <- resHept[order(resHept$padj),]
resFurOrdered <- resFur[order(resFur$padj),]
summary(resHeptOrdered)
summary(resFurOrdered)

resHept_shrink <- lfcShrink(dds, contrast=c("condition","heptamidine","saline"),res=resHept, type = "ashr")
resFur_shrink <- lfcShrink(dds, contrast=c("condition","furamidine","saline"),res=resFur, type = "ashr")

par(mfrow=c(1,2))
plotMA(resHept_shrink, alpha=0.1, ylim=c(-4.5,4.5), xlim=c(0.1,10000000), main="Heptamidine vs Saline")
plotMA(resFur_shrink, alpha=0.1, ylim=c(-4.5,4.5), xlim=c(0.1,10000000), main="Furamidine vs Saline")


#identify individual genes by clicking on plot
#  idx <- identify(resHept$baseMean, resHept$log2FoldChange) #click on plot, esc to exit
#  rownames(resHept)[idx]

#comparing various normalization methods, more details here: https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#moreshrink
resultsNames(dds)
par(mfrow=c(1,3))
H_norm <- lfcShrink(dds, coef=3, type = "normal")
H_apeglm<- lfcShrink(dds, coef=3, type = "apeglm")
H_ashr <- lfcShrink(dds, coef=3, type = "ashr")
plotMA(H_norm, alpha=0.1, ylim=c(-4.5,4.5), xlim=c(0.1,10000000), main="normal")
plotMA(H_apeglm, alpha=0.1, ylim=c(-4.5,4.5), xlim=c(0.1,10000000), main="apeglm")
plotMA(H_ashr, alpha=0.1, ylim=c(-4.5,4.5), xlim=c(0.1,10000000), main="ashr")

#Plotting normalized counts of specific genes
par(mfrow=c(2,2))
plotCounts(dds, gene = which.min(resHept$padj), intgroup = "condition")
plotCounts(dds, gene = which.min(resFur$padj), intgroup = "condition")
plotCounts(dds, gene = "ENSMUSG00000005506", intgroup = "condition") #CELF1
plotCounts(dds, gene = "ENSMUSG00000027763", intgroup = "condition") # MBNL

library(ggplot2)
data_count <- plotCounts(dds, gene="ENSMUSG00000074264", intgroup=c("condition"), 
                         returnData=TRUE)
ggplot(data_count, aes(x=condition, y=count)) + scale_y_log10() + 
  geom_point(position=position_jitter(width = .1, height = 0),size=3)

#transform values
rld <- rlog(dds, blind=FALSE)

#PCA
library("ggplot2")
par(mfrow=c(1,1))
plotPCA(rld, intgroup=c("condition"))

#PCA: with label names
plotPCA(rld, intgroup=c("condition")) + geom_text(aes(label=name))

#PCA: ggplot
data <- plotPCA(rld, intgroup=c("condition"), returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))
ggplot(data, aes(PC1, PC2, color=condition)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  geom_text(aes(label=substring(name,nchar(name))), nudge_x = 1, nudge_y = 1)

#To remove row 6 (hept3) from sample table, paste the next line onto the console without #
# sampleTable <- sampleTable[-6,]
#Now re-run starting from line 62

#Comparisons
sig_fur <- subset(resFurOrdered, resFurOrdered$padj <= 0.1)
sig_hept <- subset(resHeptOrdered, resHeptOrdered$padj <= 0.1)

fur_hept <- merge(as.data.frame(sig_fur), as.data.frame(sig_hept), by="row.names", suffixes = c(".fur",".hept"))

plot(fur_hept$log2FoldChange.fur,fur_hept$log2FoldChange.hept,
     col = "grey10", type = "p", cex = .4, main = "fur v hept", ylim=c(-5,5), xlim=c(-5,5))
abline(a=0,b=1,col="grey")
abline(h=0,col="red")
abline(v=0, col="red")

