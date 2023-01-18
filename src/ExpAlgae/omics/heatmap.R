










getmm = read.delim("C:/Users/Bisbii/PythonProjects/ExpAlgae/data/getmm.tsv", header = TRUE, row.names = "GeneID")


best = getmm[degs[0:50,]$gene_id,]

heatmap(as.matrix(getmm[degs[0:50,]$gene_id,]))
