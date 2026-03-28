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

ggplot(sample_annot, aes(x = diagnosis, y = age)) +
  geom_boxplot() +
  labs(title = "Amžius pagal diagnozę",
       x = "Diagnozė",
       y = "Amžius")

