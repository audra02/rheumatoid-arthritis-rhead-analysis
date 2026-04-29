# ------------------------------------------------
# PAPILDOMŲ OUTLIERIŲ ANALIZĖ
# ------------------------------------------------
# Iš originalaus rhead.rds pašalinami žinomi outlier mėginiai,
# sukuriama dendrograma su ID ir anotacijomis.
# Tuomet, pagal dendrogramą atrenkami papildomi outlieriai
# apie kuriuos minėjo dėstytojas 1 užd. ataskaitoje
# ir sukuriama nauja dendrograma po jų pašalinimo.
#
# Sukuriami failai:
# - ../../plots/task2/dendrogram_with_ids.png
# - ../../plots/task2/dendrogram_after_removal.png
# ------------------------------------------------

library(annmatrix)
library(WGCNA)

obj <- readRDS("rhead.rds")

dir.create("../../plots/task2", showWarnings = FALSE, recursive = TRUE)

known_outliers <- c(
  "GSM3833612_9704031135_R03C01",
  "GSM3833638_9704031135_R02C01",
  "GSM3833716_9259684070_R01C02",
  "GSM3833615_9704031135_R05C01",
  "GSM3833773_9274651074_R06C01"
)

obj1 <- obj[, !colnames(obj) %in% known_outliers]

hc1 <- hclust(as.dist(1 - cor(obj1)), method = "complete")

png("../../plots/task2/dendrogram_with_ids.png",
    width = 3000, height = 1600, res = 180)

plot(
  hc1,
  labels = colnames(obj1),
  cex = 0.35,
  main = "Dendrograma su mėginių ID",
  xlab = "",
  sub = ""
)

dev.off()


extra_outliers <- c(
  "GSM3833617_9704031129_R05C01",
  "GSM3833483_9406922147_R06C02",
  "GSM3833453_9406922031_R05C02",
  "GSM3833478_9407201070_R06C02",
  "GSM3833640_9704031140_R01C01"
)

obj2 <- obj1[, !colnames(obj1) %in% extra_outliers]

hc2 <- hclust(as.dist(1 - cor(obj2)), method = "complete")

colors <- data.frame(
  Celltype  = labels2colors(as.numeric(factor(obj2$celltype))),
  Diagnosis = labels2colors(as.numeric(factor(obj2$diagnosis))),
  Plate     = labels2colors(as.numeric(factor(obj2$plate)))
)

png("../../plots/task2/dendrogram_after_removal.png",
    width = 2600, height = 1600, res = 180)

plotDendroAndColors(
  dendro = hc2,
  colors = colors,
  groupLabels = colnames(colors),
  dendroLabels = FALSE,
  main = "Dendrograma po outlier pašalinimo"
)

dev.off()