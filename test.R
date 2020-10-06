## Written by Mercedeh Movassagh <mercedeh@ds.dfci.harvard.edu>, April 2020

#devtools::document()
#devtools::load_all()
library(dplyr)
library(mirTarRnaSeq)


#Part1

DiffExp<-read.table("test/EBV_mRNA.txt", as.is=TRUE, header=T, row.names=1)
miRNAExp<-read.table("test/EBV_miRNA.txt", as.is=TRUE, header=T, row.names=1)

#tzTransFunction
#DiffExpmiRNA <- tzTrans(miRNAExp)
#DiffExpmRNA <- tzTrans(DiffExp)


# get miranda data
#Add thershold for other variables
miranda <- getInputSpecies("Epstein_Barr", threshold = 140)
DiffExpmRNASub <- miRanComp(DiffExp, miranda)


miRNA_select<-c("ebv-mir-bart9-5p")
#PlotDensityAllmRNA
MatDiffExp<-as.vector(as.matrix(DiffExp))*10
plot(density(MatDiffExp))
#HowManyZeros
table(MatDiffExp==0)

#combiner
Combine<-combiner(DiffExp,miRNAExp,miRNA_select)

#GeneVari
geneVariant<-geneVari(Combine,miRNA_select)

#PlotDataDistribution
library(purrr)
library(tidyr)
Combine[,1:9] %>%
  keep(is.numeric) %>%
  gather() %>%
  ggplot(aes(value)) +
  facet_wrap(~ key, scales = "free") +
  geom_density()


#ActualModelRunForNegativeBinomialModel
MRun<- runModelsNb(Combine,geneVariant,miRNA_select,scale=10)

##ActualModelRunForNegativeBinomialModelZreoInflated
MRun<- runModelsZInf(Combine,geneVariant,miRNA_select, scale=10, dist="poisson")
#Seems to work for EBV mRNA
MRun<- runModelsZInf(Combine,geneVariant,miRNA_select, scale=10, dist="negbin")

#RunCombinedModelForNB
MRun<-runModelsNbCombined(Combine, geneVariant, miRNA_select,  scale = 10)


###PValueDoesNotWork
#ActualModelRan
MRun<- runModels(Combine,geneVariant,miRNA_select,family="poisson")
#Gaussian
MRun<- runModels(Combine,geneVariant,miRNA_select,family = "gaussian")


# NOTE: examples for models_filter()

# filter model based on pvalue (i.e., in MRun$pvalues)
MRun <- modelsFilter(MRun, pvalues < 0.05)
# if there is more than one variable in the model (i.e., length(miRNA_select) > 1)
# in that case there will be multiple pvalues for each model, so we need to use
# sapply in that case.
# at least one significant feature in the model
MRun <- modelsFilter(MRun, sapply(pvalues, function(x) any(x < 0.05)))
# all features in the model significant.
MRun <- modelsFilter(MRun, sapply(pvalues, function(x) all(x < 0.05)))


#FDR sig Genes
FDRModel<-fdrSig(MRun, value=0.1,method="fdr")
table(FDRModel$FDR_significant)
SigFDRGenes<-as.data.frame(cbind(FDRModel$pvalues,FDRModel$FDR_Value))
names(SigFDRGenes)<-c("P Value","FDR")

## make a plots
mymodel <- FDRModel$all_models[[1]]  # pick the first
plotFit(mymodel)
plotResiduals(mymodel)
plotTerms(mymodel)

## make interaction

#miRNA selection
miRNA_select<-as.vector(c("ebv-mir-BART2","ebv-mir-BART9"))

#combiner
Combine<-combiner(DiffExpmRNA,DiffExpmiRNA,miRNA_select)

#GeneVari
geneVariant<-geneVari(Combine,miRNA_select)

MRun<- runModels(Combine,geneVariant,miRNA_select,mode="inter")
mymodel <- MRun$all_models[[6]]  # pick the first
plotFit(mymodel)
plotResiduals(mymodel)
plotTerms(mymodel)

#ForAllmiRNAs
miranda <- getInputSpecies("Epstein_Barr", threshold = 172)
DiffExpmRNASub <- miRanComp(DiffExpmRNA, miranda)
uni_miRanda<-unique(miranda$V1)
uni_miRanda<-uni_miRanda[-28]
mirtrial <- getInputSpecies("Epstein_Barr",threshold = 140)

##Note you are doing only negative but there might be positive relationships
testme <- run_all_mirna_models(uni_miRanda,DiffExpmRNA,DiffExpmiRNA, mirtrial, prob=0.90,
                               method="fdr",all_coeff=TRUE,mode="inter")
everything1 <- testme[sapply(testme, function(x) nrow(x$SigFDRGenes)) > 0]

everything2<-rsquRes(everything1)
DrawCorPlot(everything2)


#Part2
# complicated way of saying: load the files from test/ directory.
files <- local({
  filenames <- list.files(path="test", pattern="^.*\\.txt$", full.names=T)
  files <- lapply(filenames, read.table, as.is=T, header=T, sep="\t")
  names(files) <- gsub("^.*/(.*)\\.txt$", "\\1", filenames)
  return(files)
})

# get mRNAs
mrna_files <- files[grep("^mRNA", names(files))]
#U can also import individually and feed it in like list[(mRNA1,mRNA2,mRNA)]
mrna <- one2OneRnaMiRNA(mrna_files, pthreshold = 0.05)$foldchanges

# get miRNAs
mirna_files <- files[grep("^miRNA", names(files))]
mirna <- one2OneRnaMiRNA(mirna_files)$foldchanges

# get correlations
corr_0 <- corMirnaRna(mrna, mirna)

# make background correlation distr.
outs <- sampCorRnaMirna(mrna, mirna)

# plot density plots for background and corrs in our data
mirRnaDensityCor(corr_0, outs)

# get correlations below threshold
sig_corrs <- threshSig(corr_0, outs)

# get miranda data
miranda <- getInputSpecies("Mouse", threshold = 150)

# make heatmap
newcorr <- corMirnaRnaMiranda(mrna, mirna, -0.7, miranda)
mirRnaHeatmap(newcorr,upper_bound =-0.69)

# now intersect the "significant" correlations with miranda
results <- miRandaIntersect(sig_corrs, outs, mrna, mirna, miranda)

# or...
p<-mirRnaHeatmap(results$corr,upper_bound =0)

####Section3
library(dplyr)
files <- local({
  filenames <- list.files(path="test", pattern="^.*\\.txt$", full.names=T)
  files <- lapply(filenames, read.table, as.is=T, header=T, sep="\t")
  names(files) <- gsub("^.*/(.*)\\.txt$", "\\1", filenames)
  return(files)
})

mirna_files <- files[grep("^miRNA0_2", names(files))]#Only 0-2  time point difference
mrna_files <- files[grep("^mRNA0_2", names(files))]#Only 0-2 time point difference

# get FoldChangesFile
mrna <- one2OneRnaMiRNA(mrna_files, pthreshold = 0.05)$foldchanges
mirna <- one2OneRnaMiRNA(mirna_files)$foldchanges

# get interrelation
inter0 <- twoTimePoint(mrna, mirna)

# make background correlation distr.
outs <- twoTimePointSamp(mrna, mirna,Shrounds = 10 )
# plot density polots for background and corrs in our data
mirRnaDensityInter(inter0, outs)

# get correlations below threshold
sig_InterR <- threshSigInter(inter0, outs)

# get miranda data
miranda <- getInputSpecies("Mouse", threshold = 180)
# now intersect the "significant" correlations with miranda

###Sth Wrong
results <- mirandaIntersectInter(sig_InterR, outs, mrna, mirna, miranda)
#Make a data frame for results
final_results <- finInterResult(results)
#Draw Par Plots for final results dataframe
drawInterPlots(mrna,mirna,final_results)

