setwd("data-raw")

miRNA0_2 <- read.table("test/miRNA0_2.txt", sep="\t", header=TRUE, as.is=TRUE)
miRNA2_5 <- read.table("test/miRNA2_5.txt", sep="\t", header=TRUE, as.is=TRUE)
miRNA0_5 <- read.table("test/miRNA0_5.txt", sep="\t", header=TRUE, as.is=TRUE)
mRNA0_2 <- read.table("test/mRNA0_2.txt", sep="\t", header=TRUE, as.is=TRUE)
mRNA2_5 <- read.table("test/mRNA2_5.txt", sep="\t", header=TRUE, as.is=TRUE)
mRNA0_5 <- read.table("test/mRNA0_5.txt", sep="\t", header=TRUE, as.is=TRUE)
mRNA <- read.table("test/EBV_mRNA.txt", sep="\t", header=TRUE, as.is=TRUE)
miRNA <- read.table("test/EBV_miRNA.txt", sep="\t", header=TRUE, as.is=TRUE)

mRNA_fc <- read.table("internal/mrna_fc.txt", sep="\t", header=TRUE, as.is=TRUE, row.names = 1)
miRNA_fc <- read.table("internal/mirna_fc.txt", sep="\t", header=TRUE, as.is=TRUE, row.names = 1)

mRNA_fc2<-mRNA_fc[,1,drop=FALSE]
miRNA_fc2<-miRNA_fc[,1,drop=FALSE]
  
miRanda<-getInputSpecies("Epstein_Barr", threshold = 140)
miRNA_select <- c("ebv-mir-bart9-5p")
Combine <- combiner(mRNA, miRNA, miRNA_select)
geneVariant <- geneVari(Combine, miRNA_select)


corr_0 <- corMirnaRna(mRNA_fc, miRNA_fc,method="pearson")
outs <- sampCorRnaMirna(mRNA_fc, miRNA_fc,method="pearson",
                        Shrounds = 10, Srounds = 10)
sig_corrs <- threshSig(corr_0, outs,pvalue = 0.05)
miRandaM <- getInputSpecies("Mouse", threshold = 150)

inter0 <- twoTimePoint(mRNA_fc2, miRNA_fc2)
outs2 <- twoTimePointSamp(mRNA_fc2, miRNA_fc2,Shrounds = 5)
sig_InterR <- threshSigInter(inter0, outs2)
results <- mirandaIntersectInter(sig_InterR, outs2, mRNA_fc2, miRNA_fc2, miRandaM)
final_results <- finInterResult(results)
some_model <- runModels(Combine, geneVariant, "ebv-mir-bart9-5p")$all_models[[1]]


usethis::use_data(miRNA0_2, miRNA2_5, miRNA0_5,
                  mRNA0_2, mRNA2_5, mRNA0_5,
                  mRNA, miRNA,Combine,some_model,
                  mRNA_fc, miRNA_fc,miRanda,
                  geneVariant,mRNA_fc2,miRNA_fc2,
                  corr_0,outs,sig_corrs, sig_InterR,
                  miRandaM,inter0,outs2,results,final_results,
                  overwrite = TRUE,
                  internal = FALSE,
                  compress = TRUE)
