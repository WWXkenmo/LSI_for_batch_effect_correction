PATH <- "/PATH/to/bcc_data"
mat0 = data.table::fread(file = paste0(PATH, "/Raw_Data/bcc_scRNA_counts.txt"), sep = "\t")
geneID <- mat0$V1
mat0 <- as.data.frame(mat0[,-1])
rownames(mat0) <- geneID

# The cell information
coldata0 = read.table(
file = paste0(PATH, "/Raw_Data/bcc_all_metadata.txt"),
sep = "\t", header = T, stringsAsFactors = F
)
coldata0 = coldata0[c("cell.id", "patient", "treatment", "sort")]
# su001 pre-treatment
cell1 = coldata0$cell.id[coldata0$treatment == "pre" & coldata0$patient == "su001"]
coldata1 = coldata0[coldata0$cell.id %in% cell1, c("cell.id", "patient", "treatment", "sort")]
rownames(coldata1) = coldata1$cell.id
mat1 = mat0[,coldata1$cell.id]
keep = rowSums(mat1 > 0) > 0
mat1 = mat1[keep,]
rowdata0 = data.frame(Symbol = rownames(mat1), row.names = rownames(mat1))
coldata1 = coldata1[,-1]
## Create RISC object
# make sure the row-names of mat1 equal to row-names of rowdata0,
# and the col-names of mat1 equal to row-names of coldata1
dat1 = readscdata(mat1, coldata1, rowdata0, is.filter = F)
print(dat1)


cell1 = coldata0$cell.id[coldata0$treatment == "post" & coldata0$patient == "su001"]
coldata1 = coldata0[coldata0$cell.id %in% cell1, c("cell.id", "patient", "treatment", "sort")]
rownames(coldata1) = coldata1$cell.id
mat1 = mat0[,coldata1$cell.id]
keep = rowSums(mat1 > 0) > 0
mat1 = mat1[keep,]
rowdata0 = data.frame(Symbol = rownames(mat1), row.names = rownames(mat1))
coldata1 = coldata1[,-1]
dat2 = readscdata(mat1, coldata1, rowdata0, is.filter = F)

# su005 pre-treatment
cell1 = coldata0$cell.id[coldata0$treatment == "pre" & coldata0$patient == "su005"]
coldata1 = coldata0[coldata0$cell.id %in% cell1, c("cell.id", "patient", "treatment", "sort")]
rownames(coldata1) = coldata1$cell.id
mat1 = mat0[,coldata1$cell.id]
keep = rowSums(mat1 > 0) > 0
mat1 = mat1[keep,]
rowdata0 = data.frame(Symbol = rownames(mat1), row.names = rownames(mat1))
coldata1 = coldata1[,-1]
dat3 = readscdata(mat1, coldata1, rowdata0, is.filter = F)

# su005 post-treatment
cell1 = coldata0$cell.id[coldata0$treatment == "post" & coldata0$patient == "su005"]
coldata1 = coldata0[coldata0$cell.id %in% cell1, c("cell.id", "patient", "treatment", "sort")]
rownames(coldata1) = coldata1$cell.id
mat1 = mat0[,coldata1$cell.id]
keep = rowSums(mat1 > 0) > 0
mat1 = mat1[keep,]
rowdata0 = data.frame(Symbol = rownames(mat1), row.names = rownames(mat1))
coldata1 = coldata1[,-1]
dat4 = readscdata(mat1, coldata1, rowdata0, is.filter = F)

# su006 pre-treatment
cell1 = coldata0$cell.id[coldata0$treatment == "pre" & coldata0$patient == "su006"]
coldata1 = coldata0[coldata0$cell.id %in% cell1, c("cell.id", "patient", "treatment", "sort")]
rownames(coldata1) = coldata1$cell.id
mat1 = mat0[,coldata1$cell.id]
keep = rowSums(mat1 > 0) > 0
mat1 = mat1[keep,]
rowdata0 = data.frame(Symbol = rownames(mat1), row.names = rownames(mat1))
coldata1 = coldata1[,-1]
dat5 = readscdata(mat1, coldata1, rowdata0, is.filter = F)

# su006 post-treatment
cell1 = coldata0$cell.id[coldata0$treatment == "post" & coldata0$patient == "su006"]
coldata1 = coldata0[coldata0$cell.id %in% cell1, c("cell.id", "patient", "treatment", "sort")]
rownames(coldata1) = coldata1$cell.id
mat1 = mat0[,coldata1$cell.id]
keep = rowSums(mat1 > 0) > 0
mat1 = mat1[keep,]
rowdata0 = data.frame(Symbol = rownames(mat1), row.names = rownames(mat1))
coldata1 = coldata1[,-1]
dat6 = readscdata(mat1, coldata1, rowdata0, is.filter = F)

# su008 pre-treatment
cell1 = coldata0$cell.id[coldata0$treatment == "pre" & coldata0$patient == "su008"]
coldata1 = coldata0[coldata0$cell.id %in% cell1, c("cell.id", "patient", "treatment", "sort")]
rownames(coldata1) = coldata1$cell.id
mat1 = mat0[,coldata1$cell.id]
keep = rowSums(mat1 > 0) > 0
mat1 = mat1[keep,]
rowdata0 = data.frame(Symbol = rownames(mat1), row.names = rownames(mat1))
coldata1 = coldata1[,-1]
dat7 = readscdata(mat1, coldata1, rowdata0, is.filter = F)

# su008 post-treatment
cell1 = coldata0$cell.id[coldata0$treatment == "post" & coldata0$patient == "su008"]
coldata1 = coldata0[coldata0$cell.id %in% cell1, c("cell.id", "patient", "treatment", "sort")]
rownames(coldata1) = coldata1$cell.id
mat1 = mat0[,coldata1$cell.id]
keep = rowSums(mat1 > 0) > 0
mat1 = mat1[keep,]
rowdata0 = data.frame(Symbol = rownames(mat1), row.names = rownames(mat1))
coldata1 = coldata1[,-1]
dat8 = readscdata(mat1, coldata1, rowdata0, is.filter = F)


process0 <- function(obj0){
# Filter cells and genes
obj0 = scFilter(obj0, min.UMI = 1000, max.UMI = 40000, min.gene = 200, min.cell = 3)
# Normalize the raw counts
obj0 = scNormalize(obj0)
# Find highly variable genes
obj0 = scDisperse(obj0)
# print(length(obj0@vargene))
return(obj0)
}
dat1 = process0(dat1)
dat2 = process0(dat2)
dat3 = process0(dat3)
dat4 = process0(dat4)
dat5 = process0(dat5)
dat6 = process0(dat6)
dat7 = process0(dat7)
dat8 = process0(dat8)

var0 = Reduce(
intersect, list(
dat1@rowdata$Symbol, dat2@rowdata$Symbol, dat3@rowdata$Symbol, dat4@rowdata$Symbol,
dat5@rowdata$Symbol, dat6@rowdata$Symbol, dat7@rowdata$Symbol, dat8@rowdata$Symbol
)
)
# Choose a good reference dataset
data0 = list(dat6, dat1, dat2, dat3, dat4, dat5, dat7, dat8)
data0 = scMultiIntegrate(
objects = data0, eigens = 18, add.Id = NULL, var.gene = var0,
method = "RPCI", align = 'OLS', npc = 50, adjust = FALSE,
ncore = 4, do.fast = "AUTO"
)