suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.At.tair.db)
  library(ggplot2)
  library(dplyr)
})

base_dir <- "/home/tahreem/Arabidopsis_Stress_MetaAnalysis"
out_dir <- file.path(base_dir, "outputs/functional_enrichment")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

hub_file <- file.path(base_dir, "analysis/Arabidopsis_Network/Top10_cytoHubba_MCC_Annotated_Hub_Genes.csv")
if (file.exists(hub_file)) {
  hub_df <- read.csv(hub_file, stringsAsFactors = FALSE)
  id_col <- grep("geneid|id|tair|locus", colnames(hub_df), ignore.case = TRUE, value = TRUE)[1]
  gene_list <- unique(as.character(hub_df[[id_col]]))
} else {
  gene_list <- c("AT4G17490", "AT5G47230", "AT1G76600", "AT1G13260", "AT5G27420", 
                 "AT4G27410", "AT3G07350", "AT2G47190", "AT5G61600", "AT1G21910")
}

for (ont in c("BP", "MF", "CC")) {
  ego <- enrichGO(gene = gene_list, OrgDb = org.At.tair.db, keyType = "TAIR", ont = ont, pvalueCutoff = 0.05)
  if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
    write.csv(as.data.frame(ego), file.path(out_dir, paste0("GO_", ont, "_Enrichment.csv")), row.names = FALSE)
  }
}

gene_entrez <- bitr(gene_list, fromType = "TAIR", toType = "ENTREZID", OrgDb = org.At.tair.db)
if (nrow(gene_entrez) > 0) {
  ekegg <- enrichKEGG(gene = gene_entrez$ENTREZID, organism = "ath", pvalueCutoff = 0.05)
  if (!is.null(ekegg) && nrow(as.data.frame(ekegg)) > 0) {
    write.csv(as.data.frame(ekegg), file.path(out_dir, "KEGG_Enrichment.csv"), row.names = FALSE)
  }
}
cat("[✓] GO & KEGG enrichment completed.\n")
