## Written by Mercedeh Movassagh <mercedeh@ds.dfci.harvard.edu>, April 2020

#devtools::document()
#devtools::load_all()
library(dplyr)
library(mirTarRnaSeq)


#Part1

DiffExp<-read.table("test/EBV_high_L.txt", as.is=TRUE, header=T, row.names=1)
miRNAExp<-read.table("test/m_Stom_data_auto_miRNA.txt", as.is=TRUE, header=T, row.names=1)

#tzTransFunction
DiffExpmiRNA <- tzTrans(miRNAExp)
DiffExpmRNA <- tzTrans(DiffExp)
DiffExpmRNA <- DiffExpmRNA %>% dplyr::select(-stom_286)

# get miranda data
#Add thershold for other variables
miranda <- getInputSpecies("Epstein_Barr", threshold = 140)
DiffExpmRNASub <- miRanComp(DiffExpmRNA, miranda)

#miRNA selection
miRNA_select<-c("ebv-mir-BART8","ebv-mir-BART9")
#miRNA_select<-c("ebv-mir-BART9")

#combiner
Combine<-combiner(DiffExpmRNA,DiffExpmiRNA,miRNA_select)

#GeneVari
geneVariant<-geneVari(Combine,miRNA_select)

#ActualModelRan
MRun<- runModels(Combine,geneVariant,miRNA_select,mode="multi")

#FFDRsig
FDRModel<-fdrSig(MRun, value=0.1,method="fdr")
table(FDRModel$FDR_significant)

## make a plottings
mymodel <- FDRModel$all_models[[20]]  # pick the first
plotFit(mymodel)
plotResiduals(mymodel)
plotTerms(mymodel)

## make a plottings
MRun<- runModels(Combine,geneVariant,miRNA_select,mode="inter")
mymodel <- MRun$all_models[[6]]  # pick the first
plotFit(mymodel)
plotResiduals(mymodel)
plotTerms(mymodel)


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

# plot density polots for background and corrs in our data
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
outs <- twoTimePointSamp(mrna, mirna,Shrounds =100 )

# plot density polots for background and corrs in our data
mirRnaDensityInter(inter0, outs)

# get correlations below threshold
sig_InterR <- threshSigInter(inter0, outs)

# get miranda data
miranda <- getInputSpecies("Mouse", threshold = 150)
# now intersect the "significant" correlations with miranda
results <- miRandaIntersectInter(sig_InterR, outs, mrna, mirna, miranda)
#Make a data frame for results
final_results <- finInterResult(results)
#Draw Par Plots for final results dataframe
drawInterPlots(mrna,mirna,final_results)

