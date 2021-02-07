library(dplyr)
library(ggplot2)
library(mirTarRnaSeq)
library(reshape2)
library(tidyr)
library(viridis)

DiffExp <- read.table("test/EBV_mRNA.txt", as.is = TRUE, header = T, row.names = 1)
miRNAExp <- read.table("test/EBV_miRNA.txt", as.is = TRUE, header = T, row.names = 1)
miranda <- getInputSpecies("Epstein_Barr", threshold = 140)
DiffExpmRNASub <- miRanComp(DiffExp, miranda)
miRNA_select <- c("ebv-mir-bart9-5p")
Combine <- combiner(DiffExp, miRNAExp, miRNA_select)
geneVariant <- geneVari(Combine, miRNA_select)

j <- runModel(`LMP-1` ~ `ebv-mir-bart9-5p`, Combine, model = glm_poisson(), scale = 100)
print(modelTermPvalues(j)) # Association between between mRNA and miRNA

# Gaussian Model
blaGaus <- runModels(Combine, geneVariant, miRNA_select, family = glm_gaussian(), scale = 100) # all models, regardless of coeff sign.
#blaGaus <- runModels(Combine, geneVariant, miRNA_select, family = glm_gaussian(), scale = 100, all_coeff = FALSE) # at least one neg. coeff
#blaGaus <- runModels(Combine, geneVariant, miRNA_select, family = glm_gaussian(), scale = 100, all_coeff = TRUE)  # all coeffs neg.

plot(blaGaus[["all_models"]][["BHLF1"]])
plot(modelData(blaGaus[["all_models"]][["BHLF1"]]))

# To test AIC model performance
G <- do.call(rbind.data.frame, blaGaus[["AICvalues"]])
names(G) <- c("AIC_G")
plot(density(G$AIC_G)) # Does not seem like a very good model
GM <- melt(G)

# zeroinfl Model Poisson
blazeroinfl <- runModels(Combine, geneVariant,miRNA_select, family = glm_zeroinfl(), scale = 100)
plot(blazeroinfl[["all_models"]][["BBRF3"]]) ## Fit does not plot
plot(modelData(blazeroinfl[["all_models"]][["BDLF2"]]))
# To test AIC model performance
Zp <- do.call(rbind.data.frame, blazeroinfl[["AICvalues"]])
names(Zp) <- c("AIC_Zp")
plot(density(Zp$AIC_Zp)) # Does not seem like a very good model
ZpM <- melt(Zp)

# zeroinfl Model Negative Binomial
blazeroinflNB <- runModels(Combine, geneVariant, miRNA_select, family = glm_zeroinfl(dist = "negbin"), scale = 100)
plot(blazeroinflNB[["all_models"]][["BTRF1"]]) ## Fit does not plot
plot(modelData(blazeroinflNB[["all_models"]][["BTRF1"]]))
# To test AIC model performance
ZNB <- do.call(rbind.data.frame, blazeroinflNB[["AICvalues"]])
names(ZNB) <- c("AIC_ZNB")
plot(density(ZNB$AIC_ZNB)) # Does not seem like a very good model
ZNBM <- melt(ZNB)

# How about NB and Poisson
# Negative Binomial
blaNB <- runModels(Combine, geneVariant, miRNA_select, family = glm_nb(), scale = 100)
plot(blaNB[["all_models"]][["BALF1"]])
plot(modelData(blaNB[["all_models"]][["BALF1"]]))
B <- do.call(rbind.data.frame, blaNB[["AICvalues"]])
names(B) <- c("AIC_NB")
BM <- melt(B)
# Poission
blaPois <- runModels(Combine, geneVariant, miRNA_select, family = glm_poisson(), scale = 100)
plot(blaPois[["all_models"]][["LMP-2A"]])
plot(modelData(blaPois[["all_models"]][["LMP-2A"]]))
P <- do.call(rbind.data.frame, blaPois[["AICvalues"]])
names(P) <- c("AIC_Po")
PM <- melt(P)

bindM <- rbind(PM, BM, GM, ZpM, ZNBM)

# Overall Models Gaussian Seems to work better in this case however user can choose
# the all model option for picking the best AIC per genes and reporting only P values
# for that (smaller AIC values are generally better in a model)

p2 <- ggplot(data = bindM, aes(x = value, group = variable, fill = variable)) +
  geom_density(adjust = 1.5, alpha = .3) +
  xlim(-400, 2000)+ggtitle("Plot of of AIC for ebv-mir-bart9-5p regressed all mRNAs ")+
  ylab("Density")+ xlab ("AIC Value")
p2

# Run various models compare AIC per gene and use best model per mRNA-miRNA
# interaction model


# test glm_multi() and for each model in "all_models" get the model name.
# use 'table()' to count how many of each.

blaMulti <- runModels(Combine, geneVariant, miRNA_select, family = glm_multi(), mode="multi", scale = 10)
print(table(unlist(lapply(blaMulti$all_models, modelModelName))))

# now just w/ four model types (in all the ways they can be specified... function, the result of function call, name of model as character string)

blaMulti <- runModels(Combine, geneVariant, miRNA_select, family = glm_multi(models = list(glm_gaussian, glm_nb(), glm_zeroinfl("negbin"), "glm_poisson")), scale = 1)
print(table(unlist(lapply(blaMulti$all_models, modelModelName))))

# Run all models runAllMirnaModels
vMiRNA <- rownames(miRNAExp)

# For inter and multi mode this function only supports combinations of 2
# registerDoParallel(4)

All_miRNAs_run<-runAllMirnaModels(mirnas =vMiRNA[1:2] , DiffExpmRNA = DiffExpmRNASub, DiffExpmiRNA = miRNAExp,
                                  miranda_data = miranda,prob=0.75,cutoff=0.05,fdr_cutoff = 0.1, method = "fdr",
                                  family = glm_multi(), scale = 2, mode="multi")

hasgenes <- lapply(All_miRNAs_run, function(x) nrow(x$SigFDRGenes)) > 0
All_miRNAs_run <- All_miRNAs_run[hasgenes]
print(table(unlist(lapply(blaMulti$all_models, modelModelName))))

#Print AICs used for any of the miRNA combinations if selecting item use code bellow.
print(table(unlist(lapply(All_miRNAs_run[[1]][["FDRModel"]][["all_models"]],modelModelName))))

#Print AICs used for any of the miRNA combinations if selecting particular miRNA-combinations use code bellow.
#Note currently we only support combinations of 2 for runAllMirnaModel() if you would like more combinations of miRNAs
#you can select the runModels() function with your select miRNAs.
print(table(unlist(lapply(All_miRNAs_run[["ebv-mir-bart1-5p and ebv-mir-bart11-3p"]][["FDRModel"]][["all_models"]],modelModelName))))






