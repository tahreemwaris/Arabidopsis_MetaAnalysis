#!/usr/bin/env Rscript

# ==============================================================================
# Script Name: 04_six_way_intersection.R
# Author: Tahreem
# Description: Performs high-dimensional DEG intersection across 6 abiotic and
#              biotic stress conditions in Arabidopsis thaliana:
#              (1) Bacteria (Pseudomonas syringae)
#              (2) TCV Virus
#              (3) Drought Stress
#              (4) Salt Stress
#              (5) Heat Stress
#              (6) Whitefly Infestation
#
# Inputs: Standardized DESeq2 significant DEG tables (padj < 0.05, |log2FC| >= 1)
# Outputs:
#   - Six_Conditions_DEG_Presence_Matrix.csv (Binary membership table)
#   - Six_Conditions_DEGs_Merged_Summary.csv (Consolidated log2FC profiles)
#   - Overlap subsets (>= 3, >= 4, >= 5, 6-way core genes)
#   - Six_Conditions_UpSet_Plot.pdf / .png (Set intersection visualization)
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(readr)
})

# ------------------------------------------------------------------------------
# 1. Configuration and Path Definitions
# ------------------------------------------------------------------------------
base_dir   <- "/home/tahreem/Arabidopsis_Stress_MetaAnalysis"
output_dir <- file.path(base_dir, "analysis/intersection")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Input DEG file registry
deg_inputs <- list(
  Bacteria   = file.path(base_dir, "analysis/deseq2/bacteria/Bacteria_significant_DEGs.csv"),
  TCV        = file.path(base_dir, "analysis/deseq2/tcv/TCV_significant_DEGs.csv"),
  Drought    = file.path(base_dir, "analysis/deseq2/drought/Drought_significant_DEGs.csv"),
  Salt       = file.path(base_dir, "analysis/deseq2/salt/Salt_significant_DEGs.csv"),
  HeatStress = "/mnt/d/Arabidopsis_HeatStress_RNAseq/results/deseq2/DESeq2_significant_DEGs.csv",
  Whitefly   = "/mnt/d/white-fly-data/deseq2/DESeq2_significant_DEGs.csv"
)

# ------------------------------------------------------------------------------
# 2. Validation and Loading Functions
# ------------------------------------------------------------------------------
read_and_standardize_deg <- function(filepath, condition_name) {
  if (!file.exists(filepath)) {
    stop(sprintf("[-] Input file not found for %s: %s", condition_name, filepath))
  }
  
  df <- read.csv(filepath, stringsAsFactors = FALSE, check.names = FALSE)
  
  # Standardize primary gene identifier column
  if (!"GeneID" %in% colnames(df)) {
    if ("gene_id" %in% colnames(df)) {
      df <- rename(df, GeneID = gene_id)
    } else if ("X" %in% colnames(df)) {
      df <- rename(df, GeneID = X)
    } else {
      colnames(df)[1] <- "GeneID"
    }
  }
  
  # Ensure clean locus identifiers and clean numeric metrics
  df <- df %>%
    filter(!is.na(GeneID) & GeneID != "") %>%
    mutate(
      GeneID = trimws(GeneID),
      log2FoldChange = as.numeric(log2FoldChange),
      padj = as.numeric(padj),
      Condition = condition_name
    ) %>%
    select(GeneID, log2FoldChange, padj, Condition)
  
  message(sprintf("[+] %-12s: %d significant DEGs loaded", condition_name, nrow(df)))
  return(df)
}

# ------------------------------------------------------------------------------
# 3. Execution: Load Datasets
# ------------------------------------------------------------------------------
message("\n=== STEP 4: 6-WAY DEG INTERSECTION ANALYSIS ===")
message(sprintf("[*] Output directory: %s", output_dir))
message("[*] Ingesting condition-specific DEG lists...")

deg_list <- lapply(names(deg_inputs), function(nm) {
  read_and_standardize_deg(deg_inputs[[nm]], nm)
})
names(deg_list) <- names(deg_inputs)

# Extract gene sets per condition
gene_sets <- lapply(deg_list, function(df) unique(df$GeneID))
all_unique_genes <- sort(unique(unlist(gene_sets)))
total_genes <- length(all_unique_genes)

message(sprintf("[*] Total unique union of DEGs across all 6 conditions: %d genes", total_genes))

# ------------------------------------------------------------------------------
# 4. Binary Presence-Absence Matrix
# ------------------------------------------------------------------------------
message("[*] Constructing binary intersection presence-absence matrix...")

presence_df <- data.frame(GeneID = all_unique_genes, stringsAsFactors = FALSE)
for (cond in names(gene_sets)) {
  presence_df[[cond]] <- as.integer(presence_df$GeneID %in% gene_sets[[cond]])
}

presence_df$Intersection_Count <- rowSums(presence_df[, names(gene_sets)])
presence_df <- presence_df %>% arrange(desc(Intersection_Count), GeneID)

presence_csv <- file.path(output_dir, "Six_Conditions_DEG_Presence_Matrix.csv")
write.csv(presence_df, presence_csv, row.names = FALSE)
message(sprintf("[+] Written: %s", presence_csv))

# ------------------------------------------------------------------------------
# 5. Overlap Stratification and Subsets
# ------------------------------------------------------------------------------
message("\n--- INTERSECTION STRATIFICATION ---")
cat(sprintf("%-25s | %-12s\n", "Threshold", "Gene Count"))
cat("----------------------------------------------\n")
for (k in 6:1) {
  count_k <- sum(presence_df$Intersection_Count >= k)
  exact_k <- sum(presence_df$Intersection_Count == k)
  cat(sprintf("Present in >= %d conditions   | %-6d (Exact %d: %d)\n", k, count_k, k, exact_k))
}
cat("----------------------------------------------\n")

# Export stratified tables (>= 3, >= 4, >= 5, and all 6 conditions)
for (k in c(3, 4, 5, 6)) {
  sub_df <- presence_df %>% filter(Intersection_Count >= k)
  fname  <- file.path(output_dir, sprintf("DEGs_shared_in_at_least_%d_conditions.csv", k))
  write.csv(sub_df, fname, row.names = FALSE)
}

# ------------------------------------------------------------------------------
# 6. Comprehensive Multi-Condition Expression Merging
# ------------------------------------------------------------------------------
message("\n[*] Compiling multi-condition log2FoldChange profiles...")

lfc_matrix <- data.frame(GeneID = all_unique_genes, stringsAsFactors = FALSE)
for (cond in names(deg_list)) {
  sub_lfc <- deg_list[[cond]] %>% 
    select(GeneID, log2FoldChange) %>%
    rename(!!paste0(cond, "_log2FC") := log2FoldChange)
  lfc_matrix <- left_join(lfc_matrix, sub_lfc, by = "GeneID")
}

merged_summary <- left_join(presence_df, lfc_matrix, by = "GeneID")
merged_csv     <- file.path(output_dir, "Six_Conditions_DEGs_Merged_Summary.csv")
write.csv(merged_summary, merged_csv, row.names = FALSE)
message(sprintf("[+] Written: %s", merged_csv))

# ------------------------------------------------------------------------------
# 7. Visualization: UpSet Plot
# ------------------------------------------------------------------------------
if (requireNamespace("UpSetR", quietly = TRUE)) {
  library(UpSetR)
  message("[*] Generating 6-way UpSet intersection plot...")
  
  upset_pdf <- file.path(output_dir, "Six_Conditions_UpSet_Plot.pdf")
  upset_png <- file.path(output_dir, "Six_Conditions_UpSet_Plot.png")
  
  # PDF generation
  pdf(upset_pdf, width = 12, height = 8)
  print(upset(
    presence_df[, names(gene_sets)],
    nsets = 6,
    nintersects = 40,
    order.by = "freq",
    decreasing = TRUE,
    mainbar.y.label = "Shared Differential Genes",
    sets.x.label = "Total DEGs per Condition",
    point.size = 3.5,
    line.size = 1.2,
    text.scale = c(1.5, 1.3, 1.2, 1.2, 1.4, 1.2)
  ))
  dev.off()
  
  # High-resolution PNG generation
  png(upset_png, width = 12, height = 8, units = "in", res = 300)
  print(upset(
    presence_df[, names(gene_sets)],
    nsets = 6,
    nintersects = 40,
    order.by = "freq",
    decreasing = TRUE,
    mainbar.y.label = "Shared Differential Genes",
    sets.x.label = "Total DEGs per Condition",
    point.size = 3.5,
    line.size = 1.2,
    text.scale = c(1.5, 1.3, 1.2, 1.2, 1.4, 1.2)
  ))
  dev.off()
  
  message(sprintf("[+] Saved UpSet Plot (PDF): %s", upset_pdf))
  message(sprintf("[+] Saved UpSet Plot (PNG): %s", upset_png))
} else {
  message("[!] Package 'UpSetR' is not installed. To generate plots, run: install.packages('UpSetR')")
}

message("\n[✓] Step 4 Completed successfully.\n")
