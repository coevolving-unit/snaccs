library(Seurat)
library(WGCNA)
library(Matrix)
library(DescTools)
library(glmGamPoi)
library(scuttle)
library(Azimuth)
library(scCustomize)

## Allen M1

expr = data.table::fread('allen_M1_matrix.csv', sep=',', check.names = T)
dim(expr)
samples = expr$sample_name
genes = colnames(expr)[2:length(colnames(expr))]
expr = as.matrix(expr[,2:length(colnames(expr))])
exprt = transposeBigData(expr, blocksize = 5000)
colnames(exprt) = samples

meta = data.table::fread('allen_M1_metadata.csv', sep = ',')
meta = as.data.frame(meta)
dim(meta)
rownames(meta) = meta$sample_name
table(rownames(meta) == colnames(exprt))

allen = CreateSeuratObject(counts = exprt, meta.data = meta, min.cells = 10, min.features = 200)
saveRDS(allen, file = 'allen_M1_seurat.rds')

allen = readRDS('allen_M1_seurat.rds')
allen = subset(allen, subset = outlier_call == FALSE)

myprocess <- function( input ) {
  output <- SCTransform(input, method = "glmGamPoi", verbose = FALSE)
  output <- RunPCA(output, verbose = FALSE)
  output <- RunUMAP(output, dims = 1:50, verbose = FALSE,return.model = TRUE)
  
output.ref <- AzimuthReference(
    output,
    refUMAP = "umap",
    refDR = "pca",
    refAssay = "SCT",
    dims = 1:50,
    k.param = 31,
    plotref = "umap",
    plot.metadata = NULL,
    ori.index = NULL,
    colormap = NULL,
    assays = NULL,
    metadata = c('cluster_label', 'subclass_label', 'class_label', 'cortical_layer_label'),
    reference.version = "0.0.0",
    verbose = FALSE
  )
  
  return(output.ref)
}

allen = myprocess(allen)
ValidateAzimuthReference(allen)

saveRDS(allen, file = 'allen_M1_ref.rds')

## Allen MCA

expr = data.table::fread('allen_MCA_matrix.csv', sep=',', check.names = T)
dim(expr)
samples = expr$sample_name
genes = colnames(expr)[2:length(colnames(expr))]
expr = as.matrix(expr[,2:length(colnames(expr))])
exprt = transposeBigData(expr, blocksize = 5000)
colnames(exprt) = samples
dim(exprt)
exprt[1:10, 1:10]

meta = data.table::fread('allen_MCA_metadata.csv', sep = ',')
meta = as.data.frame(meta)
dim(meta)
rownames(meta) = meta$sample_name
table(rownames(meta) == colnames(exprt))

allen = CreateSeuratObject(counts = exprt, meta.data = meta, min.cells = 10, min.features = 200)
saveRDS(allen, file = 'allen_MCA_seurat.rds')

allen = subset(allen, subset = outlier_call == FALSE)
myprocess <- function( input ) {
  output <- SCTransform(input, method = "glmGamPoi", verbose = FALSE)
  output <- RunPCA(output, verbose = FALSE)
  output <- RunUMAP(output, dims = 1:50, verbose = FALSE,return.model = TRUE)
  
output.ref <- AzimuthReference(
    output,
    refUMAP = "umap",
    refDR = "pca",
    refAssay = "SCT",
    dims = 1:50,
    k.param = 31,
    plotref = "umap",
    plot.metadata = NULL,
    ori.index = NULL,
    colormap = NULL,
    assays = NULL,
    metadata = c('cluster_label', 'subclass_label', 'class_label', 'cortical_layer_label'),
    reference.version = "0.0.0",
    verbose = FALSE
  )
  
  return(output.ref)
}

allen = myprocess(allen)
ValidateAzimuthReference(allen)

saveRDS(allen, file = 'allen_MCA_ref.rds')

## BRAIN

# merge other cell types

allen1 = readRDS('../res_alevin/h5ad/BRAIN-oligo.rds')
allen2 = readRDS('../res_alevin/h5ad/BRAIN-OPC.rds')
#allen3 = readRDS('../res_alevin/h5ad/BRAIN-COPs.rds')
allen4 = readRDS('../res_alevin/h5ad/BRAIN-Endo.rds')
allen5 = readRDS('../res_alevin/h5ad/BRAIN-VLMC.rds')
allen6 = readRDS('../res_alevin/h5ad/BRAIN-MicroPVM.rds')
allen7 = readRDS('../res_alevin/h5ad/BRAIN-Astro.rds')

allen1 = CreateSeuratObject(counts = allen1@assays$RNA@data, meta.data = allen1@meta.data)
allen2 = CreateSeuratObject(counts = allen2@assays$RNA@data, meta.data = allen2@meta.data)
#allen3 = CreateSeuratObject(counts = allen3@assays$RNA@data, meta.data = allen3@meta.data)
allen4 = CreateSeuratObject(counts = allen4@assays$RNA@data, meta.data = allen4@meta.data)
allen5 = CreateSeuratObject(counts = allen5@assays$RNA@data, meta.data = allen5@meta.data)
allen6 = CreateSeuratObject(counts = allen6@assays$RNA@data, meta.data = allen6@meta.data)
allen7 = CreateSeuratObject(counts = allen7@assays$RNA@data, meta.data = allen7@meta.data)

obj_list = list(allen1, allen2, allen4, allen5, allen6, allen7)
allen = Merge_Seurat_List(list_seurat = obj_list, merge.data = TRUE, project = "SeuratProject")
allen[['RNA']] = JoinLayers(allen[['RNA']])
str(allen)
table(allen@meta.data$Subclass)
table(allen@meta.data$Supertype)

saveRDS(allen, file = 'BRAIN-other.rds')

# merge EXN

allen1 = readRDS('../res_alevin/h5ad/BRAIN-L23IT.rds')
allen2 = readRDS('../res_alevin/h5ad/BRAIN-L4IT.rds')
allen3 = readRDS('../res_alevin/h5ad/BRAIN-L5IT.rds')
allen4 = readRDS('../res_alevin/h5ad/BRAIN-L5ET.rds')
allen5 = readRDS('../res_alevin/h5ad/BRAIN-L6ITCar3.rds')
allen6 = readRDS('../res_alevin/h5ad/BRAIN-L56NP.rds')
allen7 = readRDS('../res_alevin/h5ad/BRAIN-L6CT.rds')
allen8 = readRDS('../res_alevin/h5ad/BRAIN-L6b.rds')
allen9 = readRDS('../res_alevin/h5ad/BRAIN-L6IT.rds')

allen1 = CreateSeuratObject(counts = allen1@assays$RNA@data, meta.data = allen1@meta.data)
allen2 = CreateSeuratObject(counts = allen2@assays$RNA@data, meta.data = allen2@meta.data)
allen3 = CreateSeuratObject(counts = allen3@assays$RNA@data, meta.data = allen3@meta.data)
allen4 = CreateSeuratObject(counts = allen4@assays$RNA@data, meta.data = allen4@meta.data)
allen5 = CreateSeuratObject(counts = allen5@assays$RNA@data, meta.data = allen5@meta.data)
allen6 = CreateSeuratObject(counts = allen6@assays$RNA@data, meta.data = allen6@meta.data)
allen7 = CreateSeuratObject(counts = allen7@assays$RNA@data, meta.data = allen7@meta.data)
allen8 = CreateSeuratObject(counts = allen8@assays$RNA@data, meta.data = allen8@meta.data)
allen9 = CreateSeuratObject(counts = allen9@assays$RNA@data, meta.data = allen9@meta.data)

obj_list = list(allen1, allen2, allen3, allen4, allen5, allen6, allen7, allen8, allen9)
allen = Merge_Seurat_List(list_seurat = obj_list, merge.data = TRUE, project = "SeuratProject")
allen[['RNA']] = JoinLayers(allen[['RNA']])
str(allen)
table(allen@meta.data$Subclass)
table(allen@meta.data$Supertype)

saveRDS(allen, file = 'BRAIN-exn.rds')

# merge INN

allen1 = readRDS('../res_alevin/h5ad/BRAIN-LAMP5.rds')
allen2 = readRDS('../res_alevin/h5ad/BRAIN-SstChodl.rds')
allen3 = readRDS('../res_alevin/h5ad/BRAIN-Pax6.rds')
allen4 = readRDS('../res_alevin/h5ad/BRAIN-Chandelier.rds')
allen5 = readRDS('../res_alevin/h5ad/BRAIN-Sncg.rds')
allen6 = readRDS('../res_alevin/h5ad/BRAIN-Lamp5Lhx6.rds')
allen7 = readRDS('../res_alevin/h5ad/BRAIN-PVALB.rds')
allen8 = readRDS('../res_alevin/h5ad/BRAIN-VIP.rds')
allen9 = readRDS('../res_alevin/h5ad/BRAIN-SST.rds')

allen1 = CreateSeuratObject(counts = allen1@assays$RNA@data, meta.data = allen1@meta.data)
allen2 = CreateSeuratObject(counts = allen2@assays$RNA@data, meta.data = allen2@meta.data)
allen3 = CreateSeuratObject(counts = allen3@assays$RNA@data, meta.data = allen3@meta.data)
allen4 = CreateSeuratObject(counts = allen4@assays$RNA@data, meta.data = allen4@meta.data)
allen5 = CreateSeuratObject(counts = allen5@assays$RNA@data, meta.data = allen5@meta.data)
allen6 = CreateSeuratObject(counts = allen6@assays$RNA@data, meta.data = allen6@meta.data)
allen7 = CreateSeuratObject(counts = allen7@assays$RNA@data, meta.data = allen7@meta.data)
allen8 = CreateSeuratObject(counts = allen8@assays$RNA@data, meta.data = allen8@meta.data)
allen9 = CreateSeuratObject(counts = allen9@assays$RNA@data, meta.data = allen9@meta.data)

obj_list = list(allen1, allen2, allen3, allen4, allen5, allen6, allen7, allen8, allen9)
allen = Merge_Seurat_List(list_seurat = obj_list, merge.data = TRUE, project = "SeuratProject")
str(allen)
allen[['RNA']] = JoinLayers(allen[['RNA']])
table(allen@meta.data$Subclass)
table(allen@meta.data$Supertype)

saveRDS(allen, file = 'BRAIN-inn.rds')

# run 

myprocess <- function( input ) {
  output <- SCTransform(input, method = "glmGamPoi", verbose = FALSE)
  output <- RunPCA(output, verbose = FALSE)
  output <- RunUMAP(output, dims = 1:50, verbose = FALSE,return.model = TRUE)
  
output.ref <- AzimuthReference(
    output,
    refUMAP = "umap",
    refDR = "pca",
    refAssay = "SCT",
    dims = 1:50,
    k.param = 31,
    plotref = "umap",
    plot.metadata = NULL,
    ori.index = NULL,
    colormap = NULL,
    assays = NULL,
    metadata = c('Supertype','Subclass'),
    reference.version = "0.0.0",
    verbose = FALSE
  )
  
  return(output.ref)
}

allen = readRDS('BRAIN-inn.rds')
allen = myprocess(allen)
ValidateAzimuthReference(allen)
saveRDS(allen, file = 'BRAIN-inn-ref.rds')

allen = readRDS('BRAIN-other.rds')
allen = myprocess(allen)
ValidateAzimuthReference(allen)
saveRDS(allen, file = 'BRAIN-other-ref.rds')

allen = readRDS('BRAIN-exn.rds')
allen = myprocess(allen)
ValidateAzimuthReference(allen)
saveRDS(allen, file = 'BRAIN-exn-ref.rds')

