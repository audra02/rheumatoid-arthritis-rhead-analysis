#task: 7. Apžvalginė analizė

library(ggplot2)
library(annmatrix)

# Duomenys
obj <- readRDS("rhead_filtered.rds")
res <- readRDS("results.rds")

# Top 10 CpG pagal mažiausią p
top_cpg <- res[order(res$p_value), "CpG"][1:10]

#Metilinimo reiksmes
beta <- obj[top_cpg, ]

# Paverčiam į data frame
df <- as.data.frame(as.table(beta))
colnames(df) <- c("CpG", "Sample", "beta")

# Pridedam anotacijas
df$celltype <- obj$celltype[df$Sample]
df$diagnosis <- obj$diagnosis[df$Sample]


# Grafikas
p <- ggplot(df, aes(x = celltype, y = beta, colour = diagnosis)) +
  geom_jitter(width = 0.2, alpha = 0.6) +
  facet_wrap(~ CpG, scales = "free_y", ncol = 2) +
  labs(
    title = "Reikšmingiausių CpG metilinimo skirtumai pagal ląstelių tipą",
    x = "Ląstelių tipas",
    y = "Metilinimo lygis (beta)",
    colour = "Diagnozė"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )


ggsave(
  "../../plots/task2/cpg_celltype_analysis.png",
  p,
  width = 11,
  height = 15,
  dpi = 310
)