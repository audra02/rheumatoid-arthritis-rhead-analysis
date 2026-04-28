#task: 2. Atvaizduoti 10 pačių patikimiausių citozinų grafiškai

library(ggplot2)

# ------------------------------------------------
# 1. NUSKAITYMAS
# ------------------------------------------------

beta <- readRDS("rhead_filtered.rds")
results <- readRDS("results.rds")

sample_anno <- attr(beta, ".annmatrix.cann")
group <- tolower(sample_anno$diagnosis)

# ------------------------------------------------
# 2. ATRANKA (patikimiausi CpG)
# ------------------------------------------------

significant <- results[results$padj < 0.05, ]

top10 <- significant[order(-abs(significant$effect_size)), ][1:10, ]

# ------------------------------------------------
# 3. GRAFIKAI
# ------------------------------------------------

for (i in 1:10) {
  
  cpg <- top10$CpG[i]
  
  values <- beta[cpg, ]
  
  plot_df <- data.frame(
    value = as.numeric(values),
    group = group
  )
  
  p <- ggplot(plot_df, aes(x = group, y = value, color = group)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.4) +
    geom_jitter(width = 0.2, alpha = 0.6) +
    stat_summary(
      fun = mean,
      geom = "point",
      shape = 18,      # diamond
      size = 4,
      color = "black"
    ) +
    theme_bw() +
    labs(
      title = cpg,
      x = "Group",
      y = "Methylation",
      caption = "Black diamond indicates group mean."
    ) +
    theme(
      legend.position = "none",
      plot.caption = element_text(size = 9, hjust = 0)
    )
  
  ggsave(
    filename = paste0("plots/cpg_", i, ".png"),
    plot = p,
    width = 6,
    height = 4,
    dpi = 300
  )
}
