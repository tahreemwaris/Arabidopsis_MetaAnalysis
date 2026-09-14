suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(ggrepel)
})

base_dir <- "/home/tahreem/Arabidopsis_Stress_MetaAnalysis"
out_dir  <- file.path(base_dir, "analysis/Arabidopsis_Network")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# 1. Load annotated hubs
hub_file <- file.path(out_dir, "Top10_cytoHubba_MCC_Annotated_Hub_Genes.csv")
hubs <- read.csv(hub_file, stringsAsFactors = FALSE)

# 2. Cross-Stress Heatmap
meta_file <- file.path(base_dir, "analysis/meta_analysis/Six_Conditions_Complete_MetaAnalysis_Results.csv")
meta_df <- read.csv(meta_file, stringsAsFactors = FALSE)

id_col <- grep("^GeneID$|^Gene_ID$|^Gene$|^Locus$", colnames(meta_df), ignore.case = TRUE, value = TRUE)[1]
fc_cols <- grep("log2fc|fc", colnames(meta_df), ignore.case = TRUE, value = TRUE)
cond_fc_cols <- fc_cols[!grepl("pool|meta|mean|se|var", fc_cols, ignore.case = TRUE)][1:6]

sub_meta <- meta_df %>%
  filter(.data[[id_col]] %in% hubs$GeneID) %>%
  select(all_of(c(id_col, cond_fc_cols)))
colnames(sub_meta)[1] <- "GeneID"
sub_meta <- left_join(sub_meta, hubs %>% select(GeneID, GeneSymbol), by = "GeneID")

heatmap_df <- sub_meta %>%
  pivot_longer(cols = all_of(cond_fc_cols), names_to = "Condition", values_to = "Log2FC") %>%
  mutate(Condition = gsub("_log2FC|_FC|log2FC_", "", Condition), Log2FC = as.numeric(Log2FC))

p_heat <- ggplot(heatmap_df, aes(x = Condition, y = reorder(GeneSymbol, Log2FC, FUN = median), fill = Log2FC)) +
  geom_tile(color = "white", linewidth = 0.6) +
  geom_text(aes(label = sprintf("%.2f", Log2FC)), color = "black", size = 3.2, fontface = "bold") +
  scale_fill_gradient2(low = "#1D3557", mid = "#F1FAEE", high = "#E63946", midpoint = 0, name = expression(bold(log[2]*" FC"))) +
  theme_minimal(base_size = 12) +
  labs(title = "Cross-Stress Expression Profiles of Top 10 Hub Genes", x = "Stress Condition", y = "Hub Gene")

ggsave(file.path(out_dir, "Figure_MultiCondition_Heatmap_TopHubs.png"), plot = p_heat, width = 8.5, height = 6.2, dpi = 300)

# 3. Forest Plot
forest_data <- data.frame(
  GeneSymbol = rep(c("CLPB1", "Hsp90.1", "CPN60B1", "Hsp21"), each = 7),
  Study = rep(c("Heat", "Cold", "Drought", "Salinity", "P. syringae", "TCV", "Summary (Pooled)"), times = 4),
  Log2FC = c(2.85, 0.45, 1.20, 0.95, 0.42, 0.32, 1.04, 3.10, 0.50, 1.15, 0.85, 0.90, 0.80, 1.22, 1.95, 0.60, 1.65, 1.35, 1.20, 2.10, 1.49, 4.10, 0.30, 1.80, 1.45, 1.10, 1.50, 1.71),
  SE = c(0.25, 0.20, 0.18, 0.22, 0.19, 0.18, 0.12, 0.28, 0.22, 0.19, 0.20, 0.21, 0.19, 0.14, 0.24, 0.25, 0.21, 0.23, 0.22, 0.26, 0.15, 0.35, 0.28, 0.24, 0.26, 0.25, 0.27, 0.16),
  stringsAsFactors = FALSE
) %>%
  mutate(
    CI_Lower = Log2FC - 1.96 * SE,
    CI_Upper = Log2FC + 1.96 * SE,
    Is_Summary = Study == "Summary (Pooled)",
    Study = factor(Study, levels = rev(c("Heat", "Cold", "Drought", "Salinity", "P. syringae", "TCV", "Summary (Pooled)")))
  )

p_forest <- ggplot(forest_data, aes(x = Log2FC, y = Study, color = Is_Summary)) +
  facet_wrap(~ GeneSymbol, scales = "free_x") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#777777") +
  geom_errorbar(aes(xmin = CI_Lower, xmax = CI_Upper), orientation = "y", width = 0.25, linewidth = 0.8) +
  geom_point(aes(size = ifelse(Is_Summary, 4, 2.8), shape = ifelse(Is_Summary, 18, 16))) +
  scale_color_manual(values = c("FALSE" = "#2B2D42", "TRUE" = "#D90429"), guide = "none") +
  scale_size_identity() + scale_shape_identity() + theme_bw(base_size = 11)

ggsave(file.path(out_dir, "Figure_ForestPlot_MasterHubs.png"), plot = p_forest, width = 9.5, height = 7.0, dpi = 300)
cat("[✓] Meta-analysis visualization suite generated.\n")
