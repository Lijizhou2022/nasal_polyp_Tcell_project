##scRNA-seq data analysis
#Part1 For basic transcriptome analysis and TCR analysis(use Seurat V3.0.2) #####
## Step1: treat all single cd45+ sample data as following
library(Seurat)
sample.data<-Read10X(data.dir = paste0(data_dir,"/filtered_feature_bc_matrix"))
sample <- CreateSeuratObject(counts = sample.data,
                                     min.cells = 0,
                                     min.features = 0, 
                                     project = "sample_CD45")
sample <- NormalizeData(sample, verbose = FALSE)
sample <- FindVariableFeatures(sample, selection.method = "vst", 
                                       nfeatures = 2000, verbose = FALSE)
## Step2: merge samples (for CD45 dataset)
sample.merge<-merge(sample1,y=c(sample2,...),
                    add.cell.ids = c("sample1","sample2",...),project = "nasal")
sample.merge <- NormalizeData(sample.merge, verbose = FALSE)
sample.merge<- FindVariableFeatures(sample.merge, selection.method = "vst", 
                                    nfeatures = 2000, verbose = FALSE)
sample.merge <- ScaleData(sample.merge, verbose = FALSE)
sample.merge <- RunPCA(sample.merge, npcs = 30, verbose = FALSE)
sample.merge <- RunUMAP(sample.merge, reduction = "pca", dims = 1:30)
DimPlot(sample.merge, reduction = "tsne",group.by = "seurat_clusters",cols = all_color,label = T)+coord_fixed() # demo umap and clustering result
# quality control for each cluster
sample.merge[["percent.mt"]] <- PercentageFeatureSet(sample.merge, pattern = "^MT-")
VlnPlot(sample.merge,features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),cols =  all_color, ncol = 3,pt.size = 0)
# demo gene expression
FeaturePlot(object = sample.merge,
            features = c(gene-list
            ), ncol = 4
)
DotPlot(object = sample.merge,
        features = c(gene-list
        ))+RotatedAxis()

## Step3 subset and integrate data
#subset ab-T cells
Tcell_merge<-subset(sample.merge,idents = c(0:4,6,14,17,18,21,27,29))
#subset genes
all_genes<-rownames(Tcell_merge)
subset_genes<-all_genes[-grep('TR[AB]V',perl = T,rownames(Tcell_merge))]
Tcell_merge_subgene<-subset(Tcell_merge,features = subset_genes)
#integrate data
Tcell_merge_subgene_list<-SplitObject(Tcell_merge_subgene,split.by = 'orig.ident')
integr_anchor<-FindIntegrationAnchors(object.list = Tcell_merge_subgene_list)
Tcell_merge_subgene_integr<-IntegrateData(integr_anchor)
## Step4 analyze ab-T cell dataset(integrate blood data is the same as this)
Tcell_merge_subgene_integr<- FindVariableFeatures(Tcell_merge_subgene_integr, selection.method = "vst", 
                                                  nfeatures = 2000, verbose = FALSE)
Tcell_merge_subgene_integr <- ScaleData(Tcell_merge_subgene_integr, verbose = FALSE)
Tcell_merge_subgene_integr <- RunPCA(Tcell_merge_subgene_integr, npcs = 40, verbose = FALSE)
ElbowPlot(Tcell_merge_subgene_integr)
Tcell_merge_subgene_integr <- RunTSNE(Tcell_merge_subgene_integr, reduction = "pca", dims = 1:30,check_duplicates=FALSE)
Tcell_merge_subgene_integr <- RunUMAP(Tcell_merge_subgene_integr, reduction = "pca", dims = 1:10)
DimPlot(Tcell_merge_subgene_integr, reduction = "pca")
Tcell_merge_subgene_integr <- FindNeighbors(Tcell_merge_subgene_integr ,dims = 1:10)
Tcell_merge_subgene_integr <- FindClusters(Tcell_merge_subgene_integr,graph.name ='integrated_snn' ,resolution =0.7)
## Step5 TCR analysis and demo
#generate meta matrix with TCR info
nasal_total_TCR <- read.table('TCR igblast result')
meta <- merge(x = data.frame(colnames(Tcell_merge_subgene_integ),row.names = 1),y = nasal_total_TCR,
               by.x = 'row.names',by.y = 'sequence_id_H',all.x = T,sort = T)
meta$TCR_clone<-paste0(meta$germline_cdr3_H,'-',meta$germline_cdr3_LK)  #define TCR clone. 
# TCR and demo results using ggplot2 and seurat package

#Part2 For cell cluster prediction (use Seurat V4.3.0.1)#####
nasal.anchors <- FindTransferAnchors(reference =Tcell_merge_subgene_integr, query = query_data,
                                     dims = 1:30, reference.reduction = "pca")
predictions <- TransferData(anchorset = nasal.anchors, refdata = Tcell_merge_subgene_integr$integrated_snn_res.0.7,
                            dims = 1:30)
query_data <- AddMetaData(query_data, metadata = predictions)
Tcell_merge_subgene_integr<-RunUMAP(Tcell_merge_subgene_integr, dims = 1:10, reduction = "pca", return.model = TRUE,
                                    verbose = TRUE)

query_data <- MapQuery(anchorset = nasal.anchors, reference =  Tcell_merge_subgene_integr, query = query_data,
                                     refdata = list(celltype = "integrated_snn_res.0.7"), reference.reduction = "pca", reduction.model = "umap")

query_data <- IntegrateEmbeddings(anchorset = nasal.anchors, reference = Tcell_merge_subgene_integr,
                                                query = query_data, new.reduction.name = "ref.pca")


