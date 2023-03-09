library(Seurat)

#Read the counts
message("Reading counts...")
x <- read.csv("raw-counts.csv",header=TRUE)
rownames(x) <- x[,1]
x[,1] <- NULL
print(dim(x))
print(x[1:5,1:5])

#Read the metadata
message("Reading metadata...")
m <- read.csv("metadata.csv",header=TRUE)
rownames(m) <- m[,1]
colnames(m)[1] <- "sample"
print(dim(m))
print(head(m))

#Create and save Seurat object, could be changed to .h5Seurat
message("Writing seurat object...")
saveRDS(
  CreateSeuratObject(counts=t(x),meta.data=m,project="seurat",min.cells=0,min.features=0),
  "seurat.Rds"
)
