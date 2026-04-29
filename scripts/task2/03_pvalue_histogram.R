# ------------------------------------------------
# 03_pvalue_histogram.R
# ------------------------------------------------
# Atvaizduoja p-verčių histogramą (0–1 intervalas, 100 intervalų)
# ir suskaičiuoja statistiškai reikšmingų CpG skaičių.
# ------------------------------------------------

library(ggplot2)

# ------------------------------------------------
# 1. DUOMENŲ ĮKĖLIMAS
# ------------------------------------------------

results <- readRDS("results.rds")

# ------------------------------------------------
# 2. HISTOGRAMA
# ------------------------------------------------

p <- ggplot(results, aes(x = p_value)) +
  geom_histogram(
    breaks = seq(0, 1, by = 0.01),
    color = "black",
    fill = "steelblue"
  ) +
  theme_bw() +
  labs(
    title = "P-value distribution",
    x = "P-value",
    y = "Count"
  )

ggsave(
  "../../plots/task2/pvalue_histogram.png",
  p,
  width = 6,
  height = 4,
  dpi = 300
)


# ------------------------------------------------
# 3. REIKŠMINGŲ CpG SKAIČIUS
# ------------------------------------------------

n_raw <- sum(results$p_value < 0.05)
n_adj <- sum(results$p_adj < 0.05)

n_raw
n_adj