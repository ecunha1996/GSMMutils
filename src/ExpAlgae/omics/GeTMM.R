
library(limma)
library(edgeR)


args = commandArgs(trailingOnly=TRUE)

vargs <- strsplit(args, ",")

featureCountsAbsolutePath = vargs[[1]]

rawDataPath = vargs[[2]]

GeTMM_path = vargs[[3]]

groups = factor(vargs[[4]])



counts = read.delim(featureCountsAbsolutePath, header = TRUE, row.names = "GeneID")


metadata = read.delim(rawDataPath, header = TRUE, row.names = "GeneID")


rpk <- (counts/metadata$Length)

rpk.norm <- DGEList(counts=rpk,group=groups)

rpk.norm <- calcNormFactors(rpk.norm)

norm.counts.rpk_edger <- cpm(rpk.norm)

df = data.frame(rownames(norm.counts.rpk_edger), norm.counts.rpk_edger)

colnames(df)[1] = "GeneID"

write.table(df, GeTMM_path, row.names = FALSE, sep="\t", quote = FALSE)
