# ------------------------------------------------
# 02_top10_plots.R
# ------------------------------------------------
# Atvaizduoja 10 patikimiausių CpG pozicijų:
# pirmiausia atrenkamos reikšmingos pozicijos (p_adj < 0.05),
# tada pasirenkamos 10 su didžiausiu absoliučiu efekto dydžiu.
# Grafikai išsaugomi aplanke plots/.
# ------------------------------------------------

library(annmatrix)
library(ggplot2)

# ------------------------------------------------
# 1. DUOMENŲ ĮKĖLIMAS
# ------------------------------------------------

obj <- readRDS("rhead_filtered.rds")
results <- readRDS("results.rds")

dir.create("plots", showWarnings = FALSE)


# ------------------------------------------------
# 2. PATIKIMIAUSIŲ CpG ATRANKA
# ------------------------------------------------

significant <- results[results$p_adj < 0.05, ]
top10 <- significant[order(-abs(significant$effect_size)), ][1:10, ]


# ------------------------------------------------
# 3. GRAFIKAI
# ------------------------------------------------

group <- tolower(obj$diagnosis)

for (i in seq_len(nrow(top10))) {
  
  cpg <- top10$CpG[i]
  
  plot_df <- data.frame(
    value = as.numeric(obj[cpg, ]),
    group = group
  )
  
  p <- ggplot(plot_df, aes(x = group, y = value, color = group)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.4) +
    geom_jitter(width = 0.2, alpha = 0.6) +
    stat_summary(
      fun = mean,
      geom = "point",
      shape = 18,
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