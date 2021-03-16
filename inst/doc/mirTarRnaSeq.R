## ----style, echo=FALSE, results='asis'----------------------------------------
BiocStyle::markdown()

## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  dpi=300,
  warning = FALSE,
  collapse = TRUE,
  error = FALSE,
  comment = "#>",
  echo=TRUE
)
library(mirTarRnaSeq)
library(knitr)
library(rmarkdown)
library(mirTarRnaSeq)
library(dplyr)
library(ggplot2)
library(reshape2)
library(tidyr)
library(viridis)

doctype <- opts_knit$get("rmarkdown.pandoc.to")

## ----eval=TRUE, echo=TRUE-----------------------------------------------------
DiffExp<-read.table("test/EBV_mRNA.txt", as.is=TRUE, header=T, row.names=1)
miRNAExp<-read.table("test/EBV_miRNA.txt", as.is=TRUE, header=T, row.names=1)

## ----eval=TRUE, echo=TRUE-----------------------------------------------------
miRanda <- getInputSpecies("Epstein_Barr", threshold = 140)
DiffExpmRNASub <- miRanComp(DiffExp, miRanda)

## ----eval=TRUE, echo=TRUE-----------------------------------------------------
miRNA_select<-c("ebv-mir-bart9-5p")

## ----eval=TRUE, echo=TRUE-----------------------------------------------------
Combine <- combiner(DiffExp, miRNAExp, miRNA_select)
geneVariant <- geneVari(Combine, miRNA_select)

## ----eval=TRUE, echo=TRUE-----------------------------------------------------
j <- runModel(`LMP-1` ~ `ebv-mir-bart9-5p`, 
              Combine, model = glm_poisson(), 
              scale = 100)
# Association between between mRNA and miRNA
print(modelTermPvalues(j)) 

## ----gauss, eval=TRUE, echo=TRUE,out.width="60%", fig.align="center"----------
blaGaus <- runModels(Combine,
                     geneVariant, miRNA_select,
                     family = glm_gaussian(), 
                     scale = 100)
par(mfrow=c(2,2))
plot(blaGaus[["all_models"]][["BHLF1"]])
plot(modelData(blaGaus[["all_models"]][["BHLF1"]]))
#To test AIC model performance
G <- do.call(rbind.data.frame, blaGaus[["AICvalues"]]) 
names(G) <- c("AIC_G")
#Low values seems like a reasonable model
plot(density(G$AIC_G)) 
GM <- melt(G)

## ----poisson, eval=TRUE, echo=TRUE, out.width="60%", fig.align="center"-------
blaPois <- runModels(Combine, 
                     geneVariant, miRNA_select, 
                     family = glm_poisson(), 
                     scale = 100)
par(mfrow=c(2,3))
plot(blaPois[["all_models"]][["LMP-2A"]])
plot(modelData(blaPois[["all_models"]][["LMP-2A"]]))
P <- do.call(rbind.data.frame, blaPois[["AICvalues"]])
names(P) <- c("AIC_Po")
PM <- melt(P)

## ----negbin, eval=TRUE, echo=TRUE, out.width="60%", fig.align="center"--------
blaNB <- runModels(Combine, 
                   geneVariant, miRNA_select, 
                   family = glm_nb(), scale = 100)
plot(modelData(blaNB[["all_models"]][["BALF1"]]))
B <- do.call(rbind.data.frame, blaNB[["AICvalues"]])
names(B) <- c("AIC_NB")
BM <- melt(B)

## ----zeroinfnegbin, eval=TRUE, echo=TRUE, out.width="60%", fig.align="center"----
blazeroinflNB <- runModels(Combine, geneVariant,
                            miRNA_select,
                            family = glm_zeroinfl(dist = "negbin"),
                            scale = 100)
# To test AIC model performance
ZNB <- do.call(rbind.data.frame, blazeroinflNB[["AICvalues"]])
names(ZNB) <- c("AIC_ZNB")
plot(density(ZNB$AIC_ZNB))
ZNBM<-melt(ZNB)

## ----zeroinfpoisson, message=FALSE, echo=TRUE, cache=FALSE, results=TRUE, out.width="60%", fig.align="center"----
blazeroinfl <- runModels(Combine, geneVariant,
                         miRNA_select, 
                         family = glm_zeroinfl(), 
                         scale = 100)
# To test AIC model performance
Zp <- do.call(rbind.data.frame, blazeroinfl[["AICvalues"]])
names(Zp) <- c("AIC_Zp")
plot(density(Zp$AIC_Zp))
ZpM <- melt(Zp)

## ----plots, eval=TRUE, echo=TRUE, fig.width = 5, fig.height=5, out.width="80%", dpi=300, fig.align="center"----
bindM <- rbind(PM, BM, GM, ZpM, ZNBM)
p2 <- ggplot(data = bindM, aes(x = value, group = variable, 
                               fill = variable)) +
  geom_density(adjust = 1.5, alpha = .3) +
  xlim(-400, 2000)+
  ggtitle("Plot of of AIC for ebv-mir-bart9-5p regressed all mRNAs ")+
  ylab("Density")+ xlab ("AIC Value")
p2

## ----message=FALSE,echo=TRUE, cache=FALSE, results=TRUE, warning=FALSE, comment=TRUE, warning=FALSE----
miRNA_select<-c("ebv-mir-bart9-5p","ebv-mir-bart6-3p")
Combine <- combiner(DiffExp, miRNAExp, miRNA_select)
geneVariant <- geneVari(Combine, miRNA_select)
MultiModel <- runModels(Combine, geneVariant,
                        miRNA_select, family = glm_multi(), 
                        mode="multi", scale = 10)
print(table(unlist(lapply(MultiModel$all_models, modelModelName)))) #print the name of the models used for the analysis

## ----message=FALSE,echo=TRUE, cache=FALSE, results=TRUE, warning=FALSE, comment=TRUE, warning=FALSE----
miRNA_select<-c("ebv-mir-bart9-5p","ebv-mir-bart6-3p")
Combine <- combiner(DiffExp, miRNAExp, miRNA_select)
geneVariant <- geneVari(Combine, miRNA_select)
InterModel <- runModels(Combine,
                        geneVariant, miRNA_select,
                      family = glm_multi(models=list(glm_gaussian, glm_poisson())),mode="inter", scale = 10)
print(table(unlist(lapply(InterModel$all_models, modelModelName)))) #print the name of the models used for the analysis

## ----message=FALSE,echo=TRUE, cache=FALSE, results=TRUE, warning=FALSE, comment=TRUE----
vMiRNA<-rownames(miRNAExp)
All_miRNAs_run<-runAllMirnaModels(mirnas =vMiRNA[1:5] ,
                                  DiffExpmRNA = DiffExpmRNASub, 
                                  DiffExpmiRNA = miRNAExp,
                                  miranda_data = miRanda,prob=0.75,cutoff=0.05,fdr_cutoff = 0.1, method = "fdr",
                                  family = glm_multi(), scale = 2, mode="multi")

hasgenes <- lapply(All_miRNAs_run, function(x) nrow(x$SigFDRGenes)) > 0
All_miRNAs_run <- All_miRNAs_run[hasgenes]
print(table(unlist(lapply(All_miRNAs_run$all_models, modelModelName))))
print(table(unlist(lapply
                   (All_miRNAs_run[[1]][["FDRModel"]][["all_models"]],
                     modelModelName))))
print(table(unlist(lapply
                   (All_miRNAs_run[["ebv-mir-bart1-5p and ebv-mir-bart11-3p"]][["FDRModel"]][["all_models"]],modelModelName))))

## ----eval=TRUE, echo=TRUE-----------------------------------------------------
files <- local({
  filenames <- list.files(path="test", pattern="^.*\\.txt$", full.names=T)
  files <- lapply(filenames, read.table, as.is=T, header=T, sep="\t")
  names(files) <- gsub("^.*/(.*)\\.txt$", "\\1", filenames)
  return(files)
})

## ----eval=TRUE, echo=TRUE-----------------------------------------------------
mrna_files <- files[grep("^mRNA", names(files))]

## ----eval=TRUE, echo=TRUE-----------------------------------------------------
mrna_files <- files[grep("^mRNA", names(files))]
mrna <- one2OneRnaMiRNA(mrna_files, pthreshold = 0.05)$foldchanges

## ----eval=TRUE, echo=TRUE-----------------------------------------------------
mirna_files <- files[grep("^miRNA", names(files))]
mirna <- one2OneRnaMiRNA(mirna_files)$foldchanges

## ----eval=TRUE, echo=TRUE-----------------------------------------------------
corr_0 <- corMirnaRna(mrna, mirna,method="pearson")

## ----eval=TRUE, echo=TRUE, results=FALSE--------------------------------------
outs <- sampCorRnaMirna(mrna, mirna,method="pearson",
                        Shrounds = 100, Srounds = 1000)

## ----eval=TRUE, echo=TRUE, fig.width = 3, fig.height=3, out.width="80%", dpi=300, fig.align="center"----
#Draw density plot
mirRnaDensityCor(corr_0, outs)

## ----eval=TRUE, echo=TRUE-----------------------------------------------------
#Identify significant correlation
sig_corrs <- threshSig(corr_0, outs,pvalue = 0.05)

## ----eval=TRUE, echo=TRUE-----------------------------------------------------
#Import concordant miRanda file
miRanda <- getInputSpecies("Mouse", threshold = 150)

## ----eval=TRUE, echo=TRUE, fig.width = 8, fig.height=5, out.width="90%", dpi=300, fig.align="center"----
#Extract your target correlations based on miRanda and correlation threshold.
newcorr <- corMirnaRnaMiranda(mrna, mirna, -0.7, miRanda)
mirRnaHeatmap(newcorr,upper_bound = -0.6)

## ----eval=TRUE, echo=TRUE, fig.width = 5, fig.height=4, out.width="80%", dpi=300, fig.align="center"----
#Make final results file for signficant correlations intersecting with miRanda file
results <- miRandaIntersect(sig_corrs, outs, mrna, mirna, miRanda)
#Draw correlation heatmap
p<- mirRnaHeatmap(results$corr,upper_bound =-0.99)
p

## ----eval=TRUE, echo=TRUE-----------------------------------------------------
files <- local({
  filenames <- list.files(path="test", pattern="^.*\\.txt$", full.names=T)
  files <- lapply(filenames, read.table, as.is=T, header=T, sep="\t")
  names(files) <- gsub("^.*/(.*)\\.txt$", "\\1", filenames)
  return(files)
})

## ----eval=TRUE, echo=TRUE-----------------------------------------------------
mirna_files <- files[grep("^miRNA0_5", names(files))]
mrna_files <- files[grep("^mRNA0_5", names(files))]

## ----eval=TRUE, echo=TRUE-----------------------------------------------------
# Parse Fold Change Files for P value and Fold Change.
mrna <- one2OneRnaMiRNA(mrna_files, pthreshold = 0.05)$foldchanges
mirna <- one2OneRnaMiRNA(mirna_files)$foldchanges

## ----eval=TRUE, echo=TRUE-----------------------------------------------------
# Estimate the miRNA mRNA FC differences for your dataset
inter0 <- twoTimePoint(mrna, mirna)

## ----eval=TRUE, echo=TRUE, message=FALSE, results=FALSE-----------------------
#Make a background distribution for your miRNA mRNA FC differences
outs <- twoTimePointSamp(mrna, mirna,Shrounds = 10 )

## ----eval=TRUE, echo=TRUE-----------------------------------------------------
#Import concordant miRanda file
miRanda <- getInputSpecies("Mouse", threshold = 140)

## ----eval=TRUE, echo=TRUE-----------------------------------------------------
#Identify miRNA mRNA relationships bellow a P value threshold, default is 0.05
sig_InterR <- threshSigInter(inter0, outs)

## ----eval=TRUE, echo=TRUE-----------------------------------------------------
#Intersect the mirRanda file with your output results
results <- mirandaIntersectInter(sig_InterR, outs, mrna, mirna, miRanda)

## ----eval=TRUE, echo=TRUE, fig.width = 4, fig.height=4, out.width="90%", dpi=300, fig.align="center"----
#Create a results file for heatmap
final_results <- finInterResult(results)
#Draw plots of miRNA mRNA fold changes for your results file
par(mar=c(4,4,2,1))
drawInterPlots(mrna,mirna,final_results)

## ----eval=TRUE, echo=TRUE, fig.width = 3, fig.height=2, out.width="80%", dpi=300, fig.align="center"----
CorRes<-results$corrs
#Draw heatmap for miRNA mRNA significant differences
#Note: you do not have to use the upper_bound function unless you want  
#investigate a particular range for miRNA mRNA differences/relationships
mirRnaHeatmapDiff(CorRes,upper_bound = 9.9)

