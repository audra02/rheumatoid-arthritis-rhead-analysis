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