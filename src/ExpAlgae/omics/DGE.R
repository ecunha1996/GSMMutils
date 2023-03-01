
library("TCC")



args = commandArgs(trailingOnly=TRUE)

vargs <- strsplit(args, ",")

GeTMMPath = vargs[[1]]

ouputPath = vargs[[2]]

groups = factor(vargs[[3]])
groups


counts = read.delim(GeTMMPath, header = TRUE, row.names = "GeneID")


tcc = new("TCC", counts, groups)

de = estimateDE(tcc, test.method="edger", FDR=0.1)

result <- getResult(de, sort = TRUE)

head(result)


table(de$estimatedDEG)

degs = result[result$estimatedDEG==1,]


write.table(degs, ouputPath, row.names = FALSE, sep="\t", quote = FALSE)