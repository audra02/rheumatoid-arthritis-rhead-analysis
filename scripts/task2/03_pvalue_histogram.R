#task: 3. P-verčių histograma

library(ggplot2)

# ------------------------------------------------
# 1. READ DATA
# ------------------------------------------------

results <- readRDS("results.rds")

# ------------------------------------------------
# 2. HISTOGRAM
# ------------------------------------------------

p <- ggplot(results, aes(x = p_value)) +
  geom_histogram(
    breaks = seq(0, 1, by = 0.01),   # 100 bins
    color = "black",
    fill = "steelblue"
  ) +
  theme_bw() +
  labs(
    title = "P-value distribution",
    x = "P-value",
    y = "Count"
  )

# ------------------------------------------------
# 3. SAVE PLOT
# ------------------------------------------------

dir.create("plots", showWarnings = FALSE)

ggsave(
  filename = "plots/pvalue_histogram.png",
  plot = p,
  width = 6,
  height = 4,
  dpi = 300
)

# ------------------------------------------------
# 4. COUNT SIGNIFICANT CpGs
# ------------------------------------------------

# raw p-values
sum(results$p_value < 0.05)

# adjusted p-values (already computed)
sum(results$padj < 0.05)