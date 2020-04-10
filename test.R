## Written by Mercedeh Movassagh <mercedeh@ds.dfci.harvard.edu>, April 2020

#devtools::document()
#devtools::load_all()
library(dplyr)
library(miRTarRNASeq)


#Part1

DiffExp<-read.table("test/EBV_high_L.txt", as.is=TRUE, header=T, row.names=1)
miRNAExp<-read.table("test/m_Stom_data_auto_miRNA.txt", as.is=TRUE, header=T, row.names=1)

#TztransFunction
DiffExpmiRNA <- Tztrans(miRNAExp)
DiffExpmRNA <- Tztrans(DiffExp)
DiffExpmRNA <- DiffExpmRNA %>% dplyr::select(-stom_286)

# get miranda data
#Add thershold for other variables
miranda <- getInputSpecies("Epstein_Barr", threshold = 60)
DiffExpmRNASub <- miRanComp(DiffExpmRNA, miranda)

#miRNA selection
miRNA_select<-c("ebv-mir-BART8","ebv-mir-BART9")
#miRNA_select<-c("ebv-mir-BART9")

#Combiner
Combine<-Combiner(DiffExpmRNA,DiffExpmiRNA,miRNA_select)

#GeneVari
gene_variant<-Gene_vari(Combine,miRNA_select)

#ActualModelRan
MRun<- Run_models(Combine,gene_variant,miRNA_select,mode="multi")

#FFDRsig
FDRModel<-FDRsig(MRun, value=0.1,method="fdr")
table(FDRModel$FDR_significant)

## make a plottings
mymodel <- FDRModel$all_models[[2]]  # pick the first
plotFit(mymodel)
plotResiduals(mymodel)
plotTerms(mymodel)

## make a plottings
MRun<- Run_models(Combine,gene_variant,miRNA_select,mode="inter")
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
mrna <- One2OnemRNAmiRNA(mrna_files, pthreshold = 0.05)$foldchanges

# get miRNAs
mirna_files <- files[grep("^miRNA", names(files))]
mirna <- One2OnemRNAmiRNA(mirna_files)$foldchanges

# get correlations
corr_0 <- CormiRNAmRNA(mrna, mirna)

# make background correlation distr.
outs <- SampCormiRNAmRNA(mrna, mirna)

# plot density polots for background and corrs in our data
miRmRNADensityCor(corr_0, outs)

# get correlations below threshold
sig_corrs <- ThreshSig(corr_0, outs)

# get miranda data
miranda <- getInputSpecies("Mouse", threshold = 60)

# make heatmap
newcorr <- CormiRNAmRNAmiRanda(mrna, mirna, -0.7, miranda)
miRmRNAHeatmap(newcorr,upper_bound =-0.69)

# now intersect the "significant" correlations with miranda
results <- MiRandaIntersect(sig_corrs, outs, mrna, mirna, miranda)

# or...
p<-miRmRNAHeatmap(results$corr,upper_bound =0)

####Section3

files <- local({
  filenames <- list.files(path="test", pattern="^.*\\.txt$", full.names=T)
  files <- lapply(filenames, read.table, as.is=T, header=T, sep="\t")
  names(files) <- gsub("^.*/(.*)\\.txt$", "\\1", filenames)
  return(files)
})

mirna_files <- files[grep("^miRNA0_2", names(files))]#Only 0-2  time point difference
mrna_files <- files[grep("^mRNA0_2", names(files))]#Only 0-2 time point difference

# get FoldChangesFile
mrna <- One2OnemRNAmiRNA(mrna_files, pthreshold = 0.05)$foldchanges
mirna <- One2OnemRNAmiRNA(mirna_files)$foldchanges

# get interrelation
inter0 <- TwoTimePoint(mrna, mirna)

# make background correlation distr.
outs <- TwoTimePointSamp(mrna, mirna,Shrounds =10 )

# plot density polots for background and corrs in our data
miRmRNADensityInter(inter0, outs)

# get correlations below threshold
sig_InterR <- ThreshSigInter(inter0, outs)

# get miranda data
miranda <- getInputSpecies("Mouse", threshold = 60)
# now intersect the "significant" correlations with miranda
results <- MiRandaIntersectInter(sig_InterR, outs, mrna, mirna, miranda)

final_results <- data.frame(
  mRNA=rownames(results$mrna),
  miRNA=rownames(results$mirna),
  FC_mRNA=results$mrna$FC1,
  FC_miRNA=results$mirna$FC1,
  pvalue=results$corrs$pvalue)

# plot density of mRNA FCs
plot(density(mrna$FC1))
abline(v=0, col="grey80", lty=2)

# plot density of miRNA FCs
plot(density(mirna$FC1))
abline(v=0, col="grey80", lty=2)

# plot final results (i.e. p < 0.05) mRNA vs miRNA FCs
with(final_results, plot(FC_mRNA, FC_miRNA))


