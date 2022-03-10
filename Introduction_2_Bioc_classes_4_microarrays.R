## ----include=FALSE--------------------------------------------------------------------------------------------------------
require(knitr)
opts_chunk$set(
concordance=FALSE, echo=TRUE,  warning=FALSE, error=FALSE, message=FALSE)


## ----readData, print=FALSE, echo=TRUE-------------------------------------------------------------------------------------
# setwd(" ") # Put  here your working directory
datadir <- "."
info <-readLines(file.path(datadir,"GSE58435_series_matrix.txt"), n=70)
rows2read <- 54743 -66 -2
x <-read.table(file.path(datadir,"GSE58435_series_matrix.txt"), skip=66, header=TRUE, sep="\t",row.names=1, nrows = rows2read)


## ----relabelX-------------------------------------------------------------------------------------------------------------
dim(x)
colnames(x) <- c(paste("Turner",1:5, sep="_"), paste("Control",1:5, sep="_"))
colnames(x)
head(x)


## ----summarize, print=FALSE,echo=TRUE-------------------------------------------------------------------------------------
round(apply(x,2, summary),3)  # Column-wise summary statistics,3)


## ----boxplot1,fig.align='center', fig.cap='',echo=F-----------------------------------------------------------------------
boxplot(x, col=c(rep("red", 5) , rep("green", 5)),main="Expression values for\n 5 Turner  and 5 Control samples",
    xlab="Slides",
    ylab="Expression", las=2, cex.axis=0.7, cex.main=0.7)
abline(0,0, col="black")


## ----pca------------------------------------------------------------------------------------------------------------------
pcX<-prcomp(t(x), scale=TRUE) 
loads<- round(pcX$sdev^2/sum(pcX$sdev^2)*100,1)


## ----plotPCA, fig=TRUE----------------------------------------------------------------------------------------------------
xlab<-c(paste("PC1",loads[1],"%"))
ylab<-c(paste("PC2",loads[2],"%"))
plot(pcX$x[,1:2],xlab=xlab,ylab=ylab, xlim=c(-150, 150))
title("Principal components (PCA)")
text(pcX$x[,1],pcX$x[,2],colnames(x), pos=4)


## ----codedendrogramcomputeHC----------------------------------------------------------------------------------------------
clust.euclid.average <- hclust(dist(t(x)),method="ward.D2")


## ----plotdendrograms, fig=T-----------------------------------------------------------------------------------------------
plot(clust.euclid.average, hang=-1)


## ----loadPackage----------------------------------------------------------------------------------------------------------
require(Biobase)


## ----ExpressionSet, fig.cap="Structure of the <tt>ExpressionSet</tt> class, showing its slots and their meaning. Reproduced from Klaus, B., & Reisenauer, S. (2018)", echo=FALSE----
knitr::include_graphics("images/Structure-of-Bioconductors-ExpressionSet-class.png")


## ----simulateData---------------------------------------------------------------------------------------------------------
expressionValues <- matrix (rnorm (300), nrow=30)
colnames(expressionValues) <- paste0("sample",1:10)
head(expressionValues)


## ----simulateCovariates---------------------------------------------------------------------------------------------------
targets <- data.frame(sampleNames = paste0("sample",1:10),
                      group=c(paste0("CTL",1:5),paste0("TR",1:5)),
                      age = rpois(10, 30), 
                      sex=as.factor(sample(c("Male", "Female"),10,replace=TRUE)),
                      row.names=1)
head(targets, n=10)


## ----simulateGeneInfo-----------------------------------------------------------------------------------------------------
myGenes <-  paste0("gene",1:30)


## ----simulateInfo---------------------------------------------------------------------------------------------------------
myInfo=list(myName="Alex Sanchez", myLab="Bioinformatics Lab",
          myContact="alex@somemail.com", myTitle="Practical Exercise on ExpressionSets")
show(myInfo)


## -------------------------------------------------------------------------------------------------------------------------
pcs <- prcomp(expressionValues)
names(pcs)
barplot(pcs$sdev)
plot(pcs$rotation[,1], pcs$rotation[,2], main="Representation of first two principal components")
text(pcs$rotation[,1], pcs$rotation[,2], targets$group, cex=0.8, pos=3)


## -------------------------------------------------------------------------------------------------------------------------
variab <- apply(expressionValues, 1, sd)
orderedGenes <- myGenes[order(variab, decreasing=TRUE)]
head(variab[order(variab, decreasing=TRUE)])
head(orderedGenes)


## ----subsetExpressions----------------------------------------------------------------------------------------------------
newExpress<- expressionValues[,-9]
newTargets <- targets[-9,]
wrongNewTargets <- targets [-10,]


## ----creaExpressionSet1---------------------------------------------------------------------------------------------------
myEset <- ExpressionSet(expressionValues)
class(myEset)
show(myEset)


## ----AnnotatedDataFrame2--------------------------------------------------------------------------------------------------
columnDesc <-  data.frame(labelDescription= c("Treatment/Control", 
                                                "Age at disease onset", 
                                                "Sex of patient (Male/Female"))
myAnnotDF <- new("AnnotatedDataFrame", data=targets, varMetadata= columnDesc)
show(myAnnotDF)


## -------------------------------------------------------------------------------------------------------------------------
phenoData(myEset) <- myAnnotDF


## ----creaEset2------------------------------------------------------------------------------------------------------------
myEset <- ExpressionSet(assayData=expressionValues, phenoData=myAnnotDF)
show(myEset)


## -------------------------------------------------------------------------------------------------------------------------
myEset <- ExpressionSet(assayData=expressionValues,
                        phenoData=myAnnotDF,
                        featureNames =myGenes)
# show(myEset)


## ----label=MIAME----------------------------------------------------------------------------------------------------------
myDesc <- new("MIAME", name= myInfo[["myName"]],
            lab= myInfo[["myLab"]],
            contact= myInfo[["myContact"]] ,
            title=myInfo[["myTitle"]])
print(myDesc)


## -------------------------------------------------------------------------------------------------------------------------
myEset <- ExpressionSet(assayData=expressionValues,
                        phenoData=myAnnotDF,
                        fetureNames =myGenes,
                        experimentData = myDesc)
# show(myEset)


## ----usingExpressionSets--------------------------------------------------------------------------------------------------
dim(exprs(myEset))
class(phenoData(myEset))
class(pData(phenoData(myEset)))
head(pData(phenoData(myEset)))
head(pData(myEset))


## -------------------------------------------------------------------------------------------------------------------------
smallEset <- myEset[1:15,c(1:3,6:8)]
dim(exprs(smallEset))
dim(pData(smallEset))
head(pData(smallEset))
all(colnames(exprs(smallEset))==rownames(pData(smallEset)))


## -------------------------------------------------------------------------------------------------------------------------
youngEset <- myEset[,pData(myEset)$age<30]
dim(exprs(youngEset))
head(pData(youngEset))


## -------------------------------------------------------------------------------------------------------------------------
if (!require(GEOquery)) {
  BiocManager::install("GEOquery")
}
require(GEOquery)
gse <- getGEO("GSE58435")
class(gse)
names(gse)
gse[[1]]
esetFromGEO <- gse[[1]]

