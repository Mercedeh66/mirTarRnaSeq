devtools::document()
devtools::load_all()

library(miRTarRNASeq)


# test my function
DifExp<-read.table("/Users/mercedeh/Desktop/miRNA_w_Sara/inputmRNAforZscore.txt",
                   header=T,row.names = 1)
z<-tztrans(DifExp)

