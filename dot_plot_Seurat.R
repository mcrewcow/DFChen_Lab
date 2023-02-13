#Upon data transfer from Python to R, we can now show the gene dot plot with x = datasets and y = populations

library(Seurat)

WT <- microgliaall[['RNA']]@data["Trem2",] * (microgliaall$stage == "WT")
KO <- microgliaall[['RNA']]@data["Trem2",] * (microgliaall$stage == "KO")
PBS <- microgliaall[['RNA']]@data["Trem2",] * (microgliaall$stage == "1168")
IGFBPL1 <- microgliaall[['RNA']]@data["Trem2",] * (microgliaall$stage == "1167")
microgliaall[['NEW']] <- CreateAssayObject(data = rbind(WT, KO, PBS, IGFBPL1))
DefaultAssay(microgliaall) <- "NEW"

DotPlot(microgliaall, features = c('WT', 'KO', 'PBS', 'IGFBPL1'), split.by = 'annotation', cols = c("#82b232","#FFB200", "#E23D28", "#B9D9EB"), dot.scale = 10)

#Same for Tlr4
