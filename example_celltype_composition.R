devtools::install_github('satijalab/seurat-data')
library(SeuratData)
options(timeout = max(300, getOption("timeout")))
InstallData("panc8")
data("panc8")
panc8 = UpdateSeuratObject(object = panc8)

source("R/main.R")

# Calculate overall composition
celltype.compositions.overall <- calc_CTcomp(object = panc8, annot.cols = "celltype")

par(mar = c(8.1, 4.1, 4.1, 2.1)) 
barplot(unlist(celltype.compositions.overall[["celltype"]][["freq"]]),
        ylab = "Relative abundance (%)",
        las = 2)

# Calculate sample-wise composition
celltype.compositions.sample_wise <- calc_CTcomp(object = panc8, annot.cols = "celltype", sample.col = "orig.ident")

df <- celltype.compositions.sample_wise[["celltype"]][["freq"]]
dflong <- reshape2::melt(t(df))

library(ggplot2)
ggplot(dflong, aes(x=Var1, y=value, color=Var1)) +
  geom_boxplot() + geom_point() + theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + labs(x = "", y = "Relative abundance (%)") + NoLegend()
