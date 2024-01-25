library(Seurat)
#library(BPCells)
library(ggplot2)
#library(rlang)
library(biomaRt)
library(org.Hs.eg.db)
#library(AnnotationDbi)

L1 <- read.csv('RawData/Obese - Human Adipose/CD45P-L1.csv', row.names = 1)
L1 <- as.matrix(L1)
L1 <- CreateSeuratObject(counts = L1, project = 'L1')

L2 <- read.csv('RawData/Obese - Human Adipose/CD45P-L2.csv', row.names = 1)
L2 <- as.matrix(L2)
L2 <- CreateSeuratObject(counts = L2, project = 'L2')

L3 <- read.csv('RawData/Obese - Human Adipose/CD45P-L3.csv', row.names = 1)
L3 <- as.matrix(L3)
L3 <- CreateSeuratObject(counts = L3, project = 'L3')

O4 <- read.csv('RawData/Obese - Human Adipose/CD45P-O1.csv', row.names = 1)
O4 <- as.matrix(O4)
O4 <- CreateSeuratObject(counts = O4, project = 'O1')

O5 <- read.csv('RawData/Obese - Human Adipose/CD45P-O2.csv', row.names = 1)
O5 <- as.matrix(O5)
O5 <- CreateSeuratObject(counts = O5, project = 'O2')

O6 <- read.csv('RawData/Obese - Human Adipose/CD45P-O3.csv', row.names = 1)
O6 <- as.matrix(O6)
O6 <- CreateSeuratObject(counts = O6, project = 'O3')

# Merge the Seurat objects
obj <- merge(L1, y = c(L2, L3, O4, O5, O6), add.cell.ids = c('L1', 'L2', 'L3', 'O1', 'O2', 'O3'), project = "Weight_Loss_Cohort")

obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")

obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 3000)
obj <- ScaleData(obj)
obj <- RunPCA(obj, features = VariableFeatures(object = obj))
ElbowPlot(obj, ndims = 30)
obj <- FindNeighbors(obj, dims = 1:30)
obj <- FindClusters(obj, resolution = 0.7)
obj <- RunUMAP(obj, dims = 1:30)

obj <- JoinLayers(obj)
m <- FindAllMarkers(obj, logfc.threshold = 1)

m$Gene <- mapIds(org.Hs.eg.db,
       keys = m$gene,
       keytype = 'ENSEMBL',
       column = 'SYMBOL', multiVals = 'CharacterList')

FeaturePlot(obj, features = 'ENSG00000153563', label = T) + NoLegend()
VlnPlot(obj, features = 'ENSG00000153563', group.by = 'orig.ident', pt.size = 0)
DotPlot(obj, features = 'ENSG00000153563', group.by = 'orig.ident')

saveRDS(obj, 'Processed Data/obese_human_adipose.rds')
