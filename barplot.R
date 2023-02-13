#Numbers received after DATASET.obs.annotation.value_counts() in Jupyter
#Dataset	Population	Number	Percentage
#WT	Inflammatory	85	4.216269841
#WT	IFN-response	251	12.45039683
#WT	Proliferative	148	7.341269841
#WT	Resident	1532	75.99206349
#IGFBPL1 KO	Inflammatory	914	69.24242424
#IGFBPL1 KO	IFN-response	212	16.06060606
#IGFBPL1 KO	Proliferative	118	8.939393939
#IGFBPL1 KO	Resident	76	5.757575758
#Glaucoma PBS	Inflammatory	542	24.07818747
#Glaucoma PBS	IFN-response	477	21.19058196
#Glaucoma PBS	Proliferative	58	2.576632608
#Glaucoma PBS	Resident	1174	52.15459796
#Glaucoma IGFBPL1	Inflammatory	72	8.737864078
#Glaucoma IGFBPL1	IFN-response	54	6.553398058
#Glaucoma IGFBPL1	Proliferative	0	0
#Glaucoma IGFBPL1	Resident	698	84.70873786


numbers <- read.table(file = "clipboard", sep = "\t", header=TRUE)
View(numbers)

numbers$Dataset <- factor(numbers$Dataset , levels = c('WT','IGFBPL1 KO','Glaucoma PBS','Glaucoma IGFBPL1'))
numbers$Population <- factor(numbers$Population , levels = c('Resident','Proliferative','IFN-response','Inflammatory'))

library(ggplot2)
ggplot(numbers, aes(x=Dataset, y=Percentage, fill=Population))+
    geom_bar(stat="identity", color="black", position = 'dodge') + scale_fill_manual(values = c("#228B22", "#DC143C","#FFAA1D", "#00FFFF")) + theme_bw()

