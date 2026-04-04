data <- readRDS("data/rhead.rds")

str(data)
dim(data)
names(data)
head(data)
summary(data)

sample_annot <- attr(data, ".annmatrix.cann")
probe_annot  <- attr(data, ".annmatrix.rann")


str(sample_annot)
table(sample_annot$celltype)
table(sample_annot$diagnosis)
table(sample_annot$celltype, sample_annot$diagnosis)
summary(sample_annot$age)
table(sample_annot$sex)
table(sample_annot$smoking)


table(sample_annot$diagnosis)
table(sample_annot$celltype)
table(sample_annot$celltype, sample_annot$diagnosis)

library(ggplot2)


hist(sample_annot$age,
     main = "Donorų amžiaus pasiskirstymas",
     xlab = "Amžius",
     ylab = "Dažnis")


# 1. Amžius pagal diagnozę
ggplot(sample_annot, aes(x = diagnosis, y = age)) +
  geom_boxplot() +
  geom_jitter(width = 0.2) +
  labs(title = "Amžius pagal diagnozę",
       x = "Diagnozė",
       y = "Amžius")

# 2. Ląstelių tipai pagal diagnozę
ggplot(sample_annot, aes(x = celltype, fill = diagnosis)) +
  geom_bar(position = "dodge") +
  labs(title = "Mėginių skaičius pagal ląstelių tipą ir diagnozę",
       x = "Ląstelių tipas",
       y = "Mėginių skaičius")

# 3. Chromosomų lentelė
chr_tab <- as.data.frame(table(probe_annot$chr))
colnames(chr_tab) <- c("chr", "count")

chr_tab$chr <- factor(chr_tab$chr,
                      levels = c(paste0("chr", 1:22), "chrX", "chrY"))

# 4. CpG pasiskirstymas chromosomose
ggplot(chr_tab, aes(x = chr, y = count)) +
  geom_col() +
  labs(title = "CpG zondų pasiskirstymas chromosomose",
       x = "Chromosoma",
       y = "CpG zondų skaičius") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))


## Paketai
install.packages("flextable")
install.packages("officer")
install.packages("dplyr")

library(flextable)
library(officer)
library(dplyr)

# Duomenys
data <- readRDS("data/rhead.rds")
sample_annot <- attr(data, ".annmatrix.cann")

# Pasitikrinam, kaip vadinasi diagnoziu grupes
unique(sample_annot$diagnosis)

# Jei reikia, pasikeisk "control" i tikra savo kontrolines grupes pavadinima
ra <- sample_annot %>% filter(diagnosis == "ra")
ctrl <- sample_annot %>% filter(diagnosis == "control")

# Pagalbines funkcijos
mean_sd <- function(x, digits = 1) {
  paste0(round(mean(x, na.rm = TRUE), digits), " ± ", round(sd(x, na.rm = TRUE), digits))
}

n_pct <- function(x, value, digits = 1) {
  n <- sum(x == value, na.rm = TRUE)
  p <- round(100 * n / sum(!is.na(x)), digits)
  paste0(n, " (", p, "%)")
}

# p reiksmes
p_age <- wilcox.test(age ~ diagnosis, data = sample_annot)$p.value
p_smoking <- chisq.test(table(sample_annot$diagnosis, sample_annot$smoking))$p.value
p_sex <- chisq.test(table(sample_annot$diagnosis, sample_annot$sex))$p.value
p_race <- chisq.test(table(sample_annot$diagnosis, sample_annot$race))$p.value

# Aprašomoji statistika
desc_table <- data.frame(
  Rodiklis = c(
    "Mėginių skaičius",
    "Amžius, vidurkis ± SD",
    "Rūkymas: taip, n (%)",
    "Rūkymas: ne, n (%)",
    "Moterų, n (%)",
    "Caucasian, n (%)",
    "QC detection, vidurkis ± SD",
    "QC purity, vidurkis ± SD",
    "QC IAC, vidurkis ± SD"
  ),
  RA = c(
    nrow(ra),
    mean_sd(ra$age),
    n_pct(ra$smoking, "yes"),
    n_pct(ra$smoking, "no"),
    n_pct(ra$sex, "F"),
    n_pct(ra$race, "caucasian"),
    mean_sd(ra$qc_detection, 3),
    mean_sd(ra$qc_purity, 3),
    mean_sd(ra$qc_iac, 3)
  ),
  Kontrole = c(
    nrow(ctrl),
    mean_sd(ctrl$age),
    n_pct(ctrl$smoking, "yes"),
    n_pct(ctrl$smoking, "no"),
    n_pct(ctrl$sex, "F"),
    n_pct(ctrl$race, "caucasian"),
    mean_sd(ctrl$qc_detection, 3),
    mean_sd(ctrl$qc_purity, 3),
    mean_sd(ctrl$qc_iac, 3)
  ),
  P_reiksme = c(
    "",
    round(p_age, 4),
    round(p_smoking, 4),
    "",
    round(p_sex, 4),
    round(p_race, 4),
    "",
    "",
    ""
  ),
  stringsAsFactors = FALSE
)

desc_table

# Flextable
ft <- flextable(desc_table)

ft <- set_header_labels(
  ft,
  Rodiklis = "Rodiklis",
  RA = "RA",
  Kontrole = "Kontrolė",
  P_reiksme = "p reikšmė"
)

ft <- bold(ft, part = "header")
ft <- align(ft, align = "center", part = "all")
ft <- align(ft, j = 1, align = "left", part = "body")
ft <- font(ft, fontname = "Times New Roman", part = "all")
ft <- fontsize(ft, size = 11, part = "all")
ft <- theme_box(ft)
ft <- autofit(ft)

ft <- set_caption(
  ft,
  caption = "Lentelė 1. Mėginių anotacijų aprašomoji statistika pagal diagnozės grupę"
)

ft

# Išsaugojimas į Word failą
doc <- read_docx()
doc <- body_add_par(
  doc,
  "Mėginių anotacijų aprašomoji statistika",
  style = "heading 1"
)
doc <- body_add_flextable(doc, ft)

print(doc, target = "reports/aprasomoji_statistika_meginiu_lygiu.docx")

#heatmap
install.packages("pheatmap")
library(pheatmap)

data <- readRDS("data/rhead_filtered.rds")
sample_annot <- attr(data, ".annmatrix.cann")

cpg_var <- apply(data, 1, var, na.rm = TRUE)
top_idx <- order(cpg_var, decreasing = TRUE)[1:500]
heatmap_data <- data[top_idx, ]

annotation_col <- data.frame(
  CellType = sample_annot$celltype,
  Diagnosis = sample_annot$diagnosis
)

rownames(annotation_col) <- colnames(heatmap_data)

pheatmap(
  heatmap_data,
  scale = "row",
  show_rownames = FALSE,
  show_colnames = FALSE,
  annotation_col = annotation_col,
  main = "500 labiausiai kintančių CpG pozicijų heatmap"
)



help("heatmap")
help("apply")


# rukymo grafikas 
table(sample_annot$diagnosis)
table(sample_annot$smoking)

#grafikas 
library(ggplot2)
library(dplyr)


# grafikas
ggplot(sample_annot, aes(x = diagnosis, fill = smoking)) +
  geom_bar(position = "fill") +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(
    title = "Rūkymo pasiskirstymas pagal diagnozę",
    x = "Diagnozė",
    y = "Proporcija",
    fill = "Rūkymas"
  ) +theme_classic()
