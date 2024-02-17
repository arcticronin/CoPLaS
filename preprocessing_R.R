setwd("/home/ronin/Dev/notebooks/thesis_notebook/DaniDatasets")
load("relevant_genes.RData")

print(typeof(relevant_genes))

genes_df <- data.frame(Gene = relevant_genes)

write.csv(genes_df, "relevant_genes.csv", row.names = FALSE)