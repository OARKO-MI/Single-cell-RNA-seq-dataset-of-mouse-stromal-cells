
# Background

Single cell analysis, start making the tsne using the normalized counts within seurat 

#The normalized counts are downloaded from [here](https://singlecell.broadinstitute.org/single_cell/study/SCP361/mouse-bone-marrow-stroma-in-homeostasis)


#The analysis was performed using Seurat/Seurat_4.1.0 within R/4.1.3. 



library(Seurat)
#setwd("/mnt/WorkingDir/P23-007/Intermediate/Seurat")



normCounts<-read.csv("stroma.TP4K.txt",sep="\t")
rownames(normCounts)<-normCounts$GENE
normCounts$GENE<-NULL


colSums(normCounts)

seuratObject <- CreateSeuratObject(normCounts)



#Find the variable features


seuratObject <- FindVariableFeatures(seuratObject, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(seuratObject), 10)

top10

plot1 <- VariableFeaturePlot(seuratObject)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

all.genes <- rownames(seuratObject)
seuratObject <- ScaleData(seuratObject, features = all.genes)


seuratObject <- RunPCA(seuratObject, features = VariableFeatures(object = seuratObject))
DimPlot(seuratObject, reduction = "pca")

seuratObject <- FindNeighbors(seuratObject, dims = 1:10)
seuratObject <- FindClusters(seuratObject, resolution = 0.86)


#UMAP


seuratObject <- RunUMAP(seuratObject, dims = 1:10)

#grDevices::tiff(filename = "umap.tiff", units="in", width=5, height=5 ,res=300)
DimPlot(seuratObject, reduction = "umap")
#dev.off()

#grDevices::tiff(filename = "Ar_umap.tiff", units="in", width=5, height=5 ,res=300)
FeaturePlot(object = seuratObject, features = 'Ar', reduction = "umap",cols = c("lightgrey", "red"))
#dev.off()

#grDevices::tiff(filename = "Cxcl12_umap.tiff", units="in", width=5, height=5 ,res=300)
FeaturePlot(object = seuratObject, features = 'Cxcl12', reduction = "umap",cols = c("lightgrey", "red"))
#dev.off()

#grDevices::tiff(filename = "Angpt1_umap.tiff", units="in", width=5, height=5 ,res=300)
FeaturePlot(object = seuratObject, features = 'Angpt1', reduction = "umap",cols = c("lightgrey", "red"))
#dev.off()


#tsnee (Online they have 17 clusters)



seuratObject<-RunTSNE(seuratObject, dims = 1:10)


#grDevices::tiff(filename = "tsne.tiff", units="in", width=5, height=5 ,res=300)
DimPlot(seuratObject, reduction = "tsne")
#dev.off()

#grDevices::tiff(filename = "Ar_tsne.tiff", units="in", width=5, height=5 ,res=300)
FeaturePlot(object = seuratObject, features = 'Ar', reduction = "tsne",cols = c("lightgrey", "red"))
#dev.off()

#grDevices::tiff(filename = "Cxcl12_tsne.tiff", units="in", width=5, height=5 ,res=300)
FeaturePlot(object = seuratObject, features = 'Cxcl12', reduction = "tsne",cols = c("lightgrey", "red"))
#dev.off()

#grDevices::tiff(filename = "Angpt1_tsne.tiff", units="in", width=5, height=5 ,res=300)
FeaturePlot(object = seuratObject, features = 'Angpt1', reduction = "tsne",cols = c("lightgrey", "red"))
#dev.off()





#Dotplot for Ar, Cxcl12, Angpt1


#grDevices::tiff(filename = "DotPlot.tiff", units="in", width=5, height=5 ,res=300)
DotPlot(seuratObject, features = c("Ar", "Cxcl12", "Angpt1"), cols = c("lightgrey", "red"))
#dev.off()



saveRDS(seuratObject, file = "stroma.TP4K_final.rds")

sessionInfo()

writeLines(capture.output(sessionInfo()), "sessionInfo.txt")



seuratObject<-readRDS("stroma.TP4K_final.rds", refhook = NULL)

# Plotting

#grDevices::tiff(filename = "Ar_violin.tiff", units="in", width=10, height=10 ,res=300)
VlnPlot(seuratObject, features = "Ar",)
#dev.off()

#grDevices::tiff(filename = "Cxcl12_violin.tiff", units="in", width=10, height=10 ,res=300)
VlnPlot(seuratObject, features = "Cxcl12")
#dev.off()

#grDevices::tiff(filename = "Angpt1_violin.tiff", units="in", width=10, height=10 ,res=300)
VlnPlot(seuratObject, features = "Angpt1")
#dev.off()



cluster11<-subset(x = seuratObject, idents = "11")
hist(as.vector(cluster11@assays$RNA["Cxcl12",]))

head(table(as.vector(cluster11@assays$RNA["Cxcl12",])),n=1)
#Cluster0 have 228
228/length(as.vector(cluster11@assays$RNA["Cxcl12",])) # 22% is 0 


cluster0<-subset(x = seuratObject, idents = "0")
hist(as.vector(cluster0@assays$RNA["Cxcl12",]))
head(table(as.vector(cluster0@assays$RNA["Cxcl12",])),n=1)
#Cluster0 have 1739
1739/length(as.vector(cluster0@assays$RNA["Cxcl12",])) # 76% is 0 


cluster8<-subset(x = seuratObject, idents = "8")
hist(as.vector(cluster8@assays$RNA["Cxcl12",]))
head(table(as.vector(cluster8@assays$RNA["Cxcl12",])),n=1)

#Cluster8 have 884
884/length(as.vector(cluster8@assays$RNA["Cxcl12",])) # 69% is 0 

grDevices::tiff(filename = "Bglap_umap.tiff", units="in", width=5, height=5 ,res=300)
FeaturePlot(object = seuratObject, features = 'Bglap', reduction = "umap",cols = c("lightgrey", "red"))
dev.off()

grDevices::tiff(filename = "Lepr_umap.tiff", units="in", width=5, height=5 ,res=300)
FeaturePlot(object = seuratObject, features = 'Lepr', reduction = "umap",cols = c("lightgrey", "red"))
dev.off()



# Ar vs Cxcl12

#grDevices::tiff(filename = "Correlation_Ar_vs_Cxcl12_allClusters.tiff", units="in", width=5, height=5 ,res=300)
FeatureScatter(object = seuratObject, feature1 = 'Ar', feature2 = 'Cxcl12')
#dev.off()

#grDevices::tiff(filename = "Correlation_Ar_vs_Cxcl12_Cluster10.tiff", units="in", width=5, height=5 ,res=300)
FeatureScatter(object = seuratObject, feature1 = 'Ar', feature2 = 'Cxcl12', cells=rownames(seuratObject@meta.data[which(seuratObject@meta.data$seurat_clusters==10),]))
#dev.off()

#grDevices::tiff(filename = "Correlation_Ar_vs_Cxcl12_Cluster2.tiff", units="in", width=5, height=5 ,res=300)
FeatureScatter(object = seuratObject, feature1 = 'Ar', feature2 = 'Cxcl12', cells=rownames(seuratObject@meta.data[which(seuratObject@meta.data$seurat_clusters==2),]))
#dev.off()

# Angpt1 vs Cxcl12

#grDevices::tiff(filename = "Correlation_Angpt1_vs_Cxcl12_allClusters.tiff", units="in", width=5, height=5 ,res=300)
FeatureScatter(object = seuratObject, feature1 = 'Ar', feature2 = 'Angpt1')
#dev.off()

#grDevices::tiff(filename = "Correlation_Angpt1_vs_Cxcl12_Cluster10.tiff", units="in", width=5, height=5 ,res=300)
FeatureScatter(object = seuratObject, feature1 = 'Ar', feature2 = 'Angpt1', cells=rownames(seuratObject@meta.data[which(seuratObject@meta.data$seurat_clusters==10),]))
#dev.off()

#grDevices::tiff(filename = "Correlation_Angpt1_vs_Cxcl12_Cluster2.tiff", units="in", width=5, height=5 ,res=300)
FeatureScatter(object = seuratObject, feature1 = 'Ar', feature2 = 'Angpt1', cells=rownames(seuratObject@meta.data[which(seuratObject@meta.data$seurat_clusters==2),]))
#dev.off()




all.markers<-FindAllMarkers(object=seuratObject)

#print(all.markers)

write.table(all.markers, file = "DE_allMarkers.txt", quote=F, sep="\t")


