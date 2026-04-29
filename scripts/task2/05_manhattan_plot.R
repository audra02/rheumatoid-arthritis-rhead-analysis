#task: 5. Manhattan grafikas

library(ggplot2)
library(annmatrix)

obj <- readRDS("rhead_filtered.rds")
res <- readRDS("results.rds")

ann <- rowanns(obj)

ann <- ann[, c("id", "chr", "pos")]

dat <- merge(
  res,
  ann,
  by.x = "CpG",
  by.y = "id"
)

dat$chr <- gsub("chr", "", as.character(dat$chr))

chromosomes <- c(as.character(1:22), "X", "Y")
dat <- dat[dat$chr %in% chromosomes, ]

dat$chr <- factor(dat$chr, levels = chromosomes)

dat <- dat[order(dat$chr, dat$pos), ]

dat$log_p <- -log10(dat$p_value)

# Visų chromosomų Manhattan grafikas
chr_len <- aggregate(pos ~ chr, dat, max)
names(chr_len)[2] <- "length"

chr_len$length <- as.numeric(chr_len$length)
chr_len$offset <- cumsum(chr_len$length) - chr_len$length

dat <- merge(dat, chr_len, by = "chr")
dat$pos_all <- dat$pos + dat$offset

axis_pos <- aggregate(pos_all ~ chr, dat, mean)

p1 <- ggplot(dat, aes(pos_all, log_p, colour = chr)) +
  geom_point(size = 0.5, alpha = 0.6) +
  scale_x_continuous(
    breaks = axis_pos$pos_all,
    labels = axis_pos$chr
  ) +
  labs(
    title = "Manhattan grafikas: visos chromosomos",
    x = "Chromosoma",
    y = "-log10(p reikšmė)"
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )

ggsave(
  "../../plots/task2/manhattan_all_chromosomes.png",
  p1,
  width = 12,
  height = 6,
  dpi = 300
)

# Randame chromosomą su mažiausia p reikše 

chr_signal <- aggregate(log_p ~ chr, dat, max)
chr_signal <- chr_signal[order(chr_signal$log_p, decreasing = TRUE), ]

top_chr <- as.character(chr_signal$chr[1])

# Tos chromosomos Manhattan grafikas

dat_chr <- dat[dat$chr == top_chr, ]

p2 <- ggplot(dat_chr, aes(pos, log_p)) +
  geom_point(size = 0.7, alpha = 0.7) +
  labs(
    title = paste("Manhattan grafikas: chromosoma", top_chr),
    x = paste("Pozicija chromosomoje", top_chr),
    y = "-log10(p reikšmė)"
  ) +
  theme_bw()

ggsave(
  ("../../plots/task2/manhattan_chr.png"),
  p2,
  width = 10,
  height = 5,
  dpi = 300
)

# 20 CpG vietų su mažiausiomis p reikšmėmis

top20 <- dat[order(dat$p_value), ]
top20 <- top20[1:20, c("CpG", "chr", "pos", "p_value", "p_adj", "effect_size")]

write.csv(
  top20,
  "../../plots/task2/manhattan_top20_cpg.csv",
  row.names = FALSE
)
