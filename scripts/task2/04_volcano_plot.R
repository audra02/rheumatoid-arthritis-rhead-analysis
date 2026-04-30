# ------------------------------------------------
# 04_volcano_plot.R
# ------------------------------------------------
# Užduotis 2 – Volcano grafikai
# Sukuria du volcano grafikus:
# 1) y = -log10(p_value)
# 2) y = -log10(p_adj)
# ------------------------------------------------

library(ggplot2)

# ------------------------------------------------
# 1. DUOMENŲ ĮKĖLIMAS
# ------------------------------------------------

results <- readRDS("results.rds")

stopifnot(all(c("CpG", "p_value", "p_adj", "effect_size") %in% colnames(results)))

dir.create("../../plots/task2", showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------
# 2. PARAMETRAI
# ------------------------------------------------

fdr_cutoff <- 0.05
effect_cutoff <- 0.10
n_labels_per_side <- 3

# ------------------------------------------------
# 3. PARUOŠIMAS
# ------------------------------------------------

plot_dat <- results
plot_dat$p_value_for_plot <- pmax(plot_dat$p_value, .Machine$double.xmin)
plot_dat$p_adj_for_plot   <- pmax(plot_dat$p_adj,   .Machine$double.xmin)

plot_dat$neg_log10_p     <- -log10(plot_dat$p_value_for_plot)
plot_dat$neg_log10_p_adj <- -log10(plot_dat$p_adj_for_plot)

plot_dat$direction <- "Nereikšminga"
plot_dat$direction[plot_dat$p_adj < fdr_cutoff & plot_dat$effect_size > 0] <- "RA labiau metilinta"
plot_dat$direction[plot_dat$p_adj < fdr_cutoff & plot_dat$effect_size < 0] <- "Kontrolė labiau metilinta"

plot_dat$direction <- factor(
  plot_dat$direction,
  levels = c("Nereikšminga", "RA labiau metilinta", "Kontrolė labiau metilinta")
)

# ------------------------------------------------
# 4. TAŠKŲ ŽYMĖJIMAS
# ------------------------------------------------
# Pirmiausia žymime statistiškai reikšmingus ir biologiškai ryškesnius taškus
# (|effect_size| >= 0.10). Jei tokių per mažai, užpildome reikšmingais pagal p_adj.

sig_big <- subset(plot_dat, p_adj < fdr_cutoff & abs(effect_size) >= effect_cutoff)
sig_all <- subset(plot_dat, p_adj < fdr_cutoff)

pick_side <- function(dat, sign = c("positive", "negative"), n = 4) {
  sign <- match.arg(sign)
  if (nrow(dat) == 0) return(dat[0, , drop = FALSE])
  if (sign == "positive") {
    dat <- dat[dat$effect_size > 0, , drop = FALSE]
  } else {
    dat <- dat[dat$effect_size < 0, , drop = FALSE]
  }
  if (nrow(dat) == 0) return(dat[0, , drop = FALSE])
  dat <- dat[order(dat$p_adj, -abs(dat$effect_size)), , drop = FALSE]
  head(dat, n)
}

label_pos <- pick_side(sig_big, "positive", n_labels_per_side)
label_neg <- pick_side(sig_big, "negative", n_labels_per_side)

# Jei iš sig_big pusėje trūksta taškų, papildome iš visų reikšmingų
if (nrow(label_pos) < n_labels_per_side) {
  extra_pos <- pick_side(sig_all[!(sig_all$CpG %in% label_pos$CpG), , drop = FALSE],
                         "positive", n_labels_per_side - nrow(label_pos))
  label_pos <- rbind(label_pos, extra_pos)
}
if (nrow(label_neg) < n_labels_per_side) {
  extra_neg <- pick_side(sig_all[!(sig_all$CpG %in% label_neg$CpG), , drop = FALSE],
                         "negative", n_labels_per_side - nrow(label_neg))
  label_neg <- rbind(label_neg, extra_neg)
}

label_dat <- rbind(label_pos, label_neg)

# Lentelė interpretacijai / ataskaitai
label_export <- label_dat[, c("CpG", "p_value", "p_adj", "effect_size", "direction")]
write.csv(label_export, "../../plots/task2/volcano_labeled_cpg.csv", row.names = FALSE)

# Papildomai išsaugome top20 pagal p_adj, po to pagal |effect|
# Šio nenaudojame ataskaitoje, tiesiog dėl įdomumo peržvelgiame daugiau
top20 <- plot_dat[order(plot_dat$p_adj, -abs(plot_dat$effect_size)), ]
top20 <- head(top20[, c("CpG", "p_value", "p_adj", "effect_size", "direction")], 20)
write.csv(top20, "../../plots/task2/volcano_top20_cpg.csv", row.names = FALSE)

# ------------------------------------------------
# 5. STILIUS
# ------------------------------------------------

base_theme <- theme_bw() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold")
  )

volcano_colors <- c(
  "Nereikšminga" = "#F8766D",
  "RA labiau metilinta" = "#00BA38",
  "Kontrolė labiau metilinta" = "#619CFF"
)

# ------------------------------------------------
# 6. VOLCANO: -log10(p)
# ------------------------------------------------

p1 <- ggplot(plot_dat, aes(x = effect_size, y = neg_log10_p, colour = direction)) +
  geom_point(size = 0.45, alpha = 0.55) +
  geom_vline(xintercept = c(-effect_cutoff, effect_cutoff), linetype = "dashed", linewidth = 0.35) +
  geom_hline(yintercept = -log10(0.05), linetype = "dotted", linewidth = 0.35) +
  geom_text(
    data = label_dat,
    aes(label = CpG),
    size = 2.8,
    check_overlap = FALSE,
    show.legend = FALSE
  ) +
  scale_colour_manual(values = volcano_colors) +
  labs(
    title = "Volcano grafikas: RA ir kontrolės CpG metilinimo skirtumai",
    x = "Efekto dydis: vidutinis beta skirtumas (RA - kontrolė)",
    y = "-log10(p reikšmė)",
    colour = NULL
  ) +
  base_theme

ggsave(
  "../../plots/task2/volcano_plot.png",
  p1,
  width = 9,
  height = 6,
  dpi = 300
)

# ------------------------------------------------
# 7. VOLCANO: -log10(p_adj)
# ------------------------------------------------

p2 <- ggplot(plot_dat, aes(x = effect_size, y = neg_log10_p_adj, colour = direction)) +
  geom_point(size = 0.45, alpha = 0.55) +
  geom_vline(xintercept = c(-effect_cutoff, effect_cutoff), linetype = "dashed", linewidth = 0.35) +
  geom_hline(yintercept = -log10(fdr_cutoff), linetype = "dotted", linewidth = 0.35) +
  geom_text(
    data = label_dat,
    aes(y = neg_log10_p_adj, label = CpG),
    size = 2.8,
    check_overlap = FALSE,
    show.legend = FALSE
  ) +
  scale_colour_manual(values = volcano_colors) +
  labs(
    title = "Volcano grafikas su FDR koreguotomis p reikšmėmis",
    x = "Efekto dydis: vidutinis beta skirtumas (RA - kontrolė)",
    y = "-log10(koreguota p reikšmė)",
    colour = NULL
  ) +
  base_theme

ggsave(
  "../../plots/task2/volcano_plot_fdr.png",
  p2,
  width = 9,
  height = 6,
  dpi = 300
)
