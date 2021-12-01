library(Seurat)

baron <- readRDS("baron-human.rds")
segerstolpe <- readRDS("segerstolpe.rds")

###remove uncertain cell types
segerstolpe.idx <- segerstolpe$cell_type1 %in% c("unclear","co-expression","not applicable","unclassified",
                                     "unclassified endocrine","dropped","alpha.contaminated","beta.contaminated",
									 "delta.contaminated","gamma.contaminated") == FALSE
segerstolpe <- segerstolpe[,segerstolpe.idx]


####generate expression profile
baron_exp <- counts(baron)
segerstolpe_exp <- counts(segerstolpe)

####normalize 
baron_exp <- CreateSeuratObject(baron_exp)
baron_exp <- NormalizeData(baron_exp, normalization.method = "LogNormalize", scale.factor = 10000)
baron_exp <- FindVariableFeatures(object = baron_exp, selection.method = 'vst', nfeatures =10000)
baron_feature <- VariableFeatures(baron_exp)
baron_exp <- as.matrix(GetAssayData(baron_exp)[baron_feature,])

segerstolpe_exp <- CreateSeuratObject(segerstolpe_exp)
segerstolpe_exp <- NormalizeData(segerstolpe_exp, normalization.method = "LogNormalize", scale.factor = 10000)
segerstolpe_exp <- FindVariableFeatures(object = segerstolpe_exp, selection.method = 'vst', nfeatures =10000)
segerstolpe_feature <- VariableFeatures(segerstolpe_exp)
segerstolpe_exp <- as.matrix(GetAssayData(segerstolpe_exp)[segerstolpe_feature,])
#####
#compute correlation matrix and perform SVD
library(rARPACK)
genes <- intersect(rownames(baron_exp),rownames(segerstolpe_exp))
cor <- cor(baron_exp[genes, ], segerstolpe_exp[genes, ])
cor_svd <- rARPACK::svds(cor, k = 30)
cell_embed <- rbind(cor_svd$u, cor_svd$v)

##generate new seurat object
pancreas <- CreateSeuratObject(cbind(baron_exp[genes, ], segerstolpe_exp[genes, ]))
pancreas$celltype <- c(as.character(as.matrix(baron$cell_type1)),segerstolpe$cell_type1)
pancreas$batch <- c(rep("baron",ncol(baron_exp)),rep("segerstolpe",ncol(segerstolpe_exp)))
colnames(cell_embed) <- paste0("PC-",1:30)
rownames(cell_embed) <- colnames(pancreas)
pancreas[['svd']] <- CreateDimReducObject(cell_embed,key="PC_")
pancreas <- RunUMAP(pancreas, dims = 1:30,reduction="svd")
p1 <- DimPlot(pancreas,reduction="umap",group.by = "batch",label=1)
p2 <- DimPlot(pancreas,reduction="umap",group.by = "celltype",label=1)
p_bc <- p1+p2

##Raw reductions
pancreas <- ScaleData(pancreas)
pancreas <- RunPCA(pancreas,npcs = 30,verbose = F,features=rownames(pancreas)) 
pancreas <- RunUMAP(pancreas, dims = 1:30,reduction="pca",reduction.name = "umap_raw")
p1 <- DimPlot(pancreas,reduction="umap_raw",group.by = "batch",label=1)
p2 <- DimPlot(pancreas,reduction="umap_raw",group.by = "celltype",label=1)
p_raw <- p1+p2