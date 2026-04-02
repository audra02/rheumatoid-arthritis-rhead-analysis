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


#aprosamoji statistika 
# duomenys
data <- readRDS("data/rhead.rds")
sample_annot <- attr(data, ".annmatrix.cann")

# funkcija vienam lasteliu tipui
make_summary_for_celltype <- function(beta_mat, annot, celltype_name) {
  
  idx <- annot$celltype == celltype_name
  sub_beta <- beta_mat[, idx, drop = FALSE]
  sub_annot <- annot[idx, , drop = FALSE]
  
  ra_idx <- sub_annot$diagnosis == "ra"
  control_idx <- sub_annot$diagnosis == "control"
  
  pvals <- numeric(nrow(sub_beta))
  med_diff <- numeric(nrow(sub_beta))
  
  for (i in 1:nrow(sub_beta)) {
    x_ra <- sub_beta[i, ra_idx]
    x_control <- sub_beta[i, control_idx]
    
    # Wilcoxon rank sum test, two-sided
    pvals[i] <- wilcox.test(x_ra, x_control, alternative = "two.sided")$p.value
    
    # kryptis pagal medianu skirtuma
    med_diff[i] <- median(x_ra) - median(x_control)
  }
  
  # BH korekcija, kad butu arciau straipsnio logikos
  qvals <- p.adjust(pvals, method = "BH")
  
  hypo <- med_diff < 0
  hyper <- med_diff > 0
  
  result <- data.frame(
    Cell_type = c(celltype_name, celltype_name),
    Direction = c("Hypomethylated", "Hypermethylated"),
    
    Raw_P = c(
      sum(pvals < 0.05 & hypo),
      sum(pvals < 0.05 & hyper)
    ),
    
    Diff_more_10 = c(
      sum(abs(med_diff) > 0.10 & pvals < 0.05 & hypo),
      sum(abs(med_diff) > 0.10 & pvals < 0.05 & hyper)
    ),
    
    Diff_1_to_10 = c(
      sum(abs(med_diff) >= 0.01 & abs(med_diff) <= 0.10 & pvals < 0.05 & hypo),
      sum(abs(med_diff) >= 0.01 & abs(med_diff) <= 0.10 & pvals < 0.05 & hyper)
    ),
    
    FDR_q = c(
      sum(qvals < 0.05 & hypo),
      sum(qvals < 0.05 & hyper)
    ),
    
    FDR_q_diff_more_1 = c(
      sum(qvals < 0.05 & abs(med_diff) > 0.01 & hypo),
      sum(qvals < 0.05 & abs(med_diff) > 0.01 & hyper)
    )
  )
  
  return(result)
}

# 4 lasteliu tipu lentele

final_table <- rbind(
  make_summary_for_celltype(data, sample_annot, "mono"),
  make_summary_for_celltype(data, sample_annot, "bcell"),
  make_summary_for_celltype(data, sample_annot, "tcd4m"),
  make_summary_for_celltype(data, sample_annot, "tcd4n")
)

final_table$Cell_type <- factor(
  final_table$Cell_type,
  levels = c("mono", "bcell", "tcd4m", "tcd4n"),
  labels = c("CD14+ monocytes", "CD19+ B cells", "CD4+ memory T cells", "CD4+ naive T cells")
)

final_table

install.packages("flextable")
install.packages("officer")

library(flextable)
library(officer)

ft <- flextable(final_table)

ft <- set_header_labels(
  ft,
  Cell_type = "Cell type",
  Direction = "",
  Raw_P = "CpGs with raw P < 0.05",
  Diff_more_10 = "CpGs with absolute median methylation difference of >10%",
  Diff_1_to_10 = "CpGs with absolute median methylation difference between 1% and 10%",
  FDR_q = "CpGs with FDR q < 0.05",
  FDR_q_diff_more_1 = "CpGs with FDR q < 0.05 and median methylation difference >1%"
)

ft <- merge_v(ft, j = "Cell_type")
ft <- valign(ft, j = "Cell_type", valign = "top")

ft <- align(ft, align = "center", part = "all")
ft <- align(ft, j = c("Cell_type", "Direction"), align = "left", part = "body")

ft <- bold(ft, part = "header")
ft <- font(ft, fontname = "Times New Roman", part = "all")
ft <- fontsize(ft, size = 11, part = "all")

ft <- border_remove(ft)
ft <- hline_top(ft, border = fp_border(width = 1.2), part = "all")
ft <- hline(ft, i = 2, border = fp_border(width = 0.8), part = "body")
ft <- hline(ft, i = 4, border = fp_border(width = 0.8), part = "body")
ft <- hline(ft, i = 6, border = fp_border(width = 0.8), part = "body")
ft <- hline_bottom(ft, border = fp_border(width = 1.2), part = "all")

ft <- autofit(ft)

doc <- read_docx()
doc <- body_add_par(doc, "Differential methylation summary by cell type", style = "Normal")
doc <- body_add_flextable(doc, ft)

print(doc, target = "reports/differential_methylation_summary_wilcox.docx")

#heatmap
install.packages("pheatmap")
library(pheatmap)

data <- readRDS("data/rhead.rds")
sample_annot <- attr(data, ".annmatrix.cann")

# 1. Suskaiciuojame kiekvienos CpG pozicijos dispersija
cpg_var <- apply(data, 1, var)

# 2. Atsirenkame į00 labiausiai kintanciu CpG
top_idx <- order(cpg_var, decreasing = TRUE)[1:500]

heatmap_data <- data[top_idx, ]

# 3. Anotacija stulpeliams
annotation_col <- data.frame(
  CellType = sample_annot$celltype,
  Diagnosis = sample_annot$diagnosis
)

rownames(annotation_col) <- colnames(heatmap_data)

# 4. Heatmap
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