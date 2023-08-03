library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)
library(gridExtra)
library(stringi)
library(ROCit)
library(png)
library(SeuratData)
library(patchwork)
####Data Preprocessing####
#RI model
slide_list <- c("c_2.h5",
                "c_3.h5",
                "c_4.h5",
                "c_7.h5",
                "c_20.h5",
                "c_34.h5",
                "c_36.h5",
                "c_39.h5",
                "c_45.h5",
                "c_51.h5")
#NRI model
slide_list <- c("b_1.h5",
                "b_18.h5",
                "a_3.h5",
                "a_15.h5")
spatial_list <- sapply(slide_list,function(slide){
  print(slide)
  expr.url <- paste0("data/matrix/",name[slide],".h5")
  spatial_object <- Seurat::Read10X_h5(filename =  expr.url )# Collect all genes coded on the mitochondrial genome
  spatial_object <- Seurat::CreateSeuratObject(counts = spatial_object, project = name[slide], assay = 'Spatial')
  mt.genes <- grep(pattern = "^MT-", x = rownames(spatial_object), value = TRUE)
  spatial_object$percent.mito <- (Matrix::colSums(spatial_object@assays$Spatial@counts[mt.genes, ])/Matrix::colSums(spatial_object@assays$Spatial@counts))*100
  #remove mt genes
  genes_to_keep <- setdiff(names(which(Matrix::rowSums(spatial_object@assays$Spatial@counts )>5)),mt.genes)
  spatial_object_subset <- subset(spatial_object,features =genes_to_keep, subset = nFeature_Spatial > 300 & percent.mito < 30)
  cat("Spots removed: ", ncol(spatial_object) - ncol(spatial_object_subset), "\n")
  cat("Genes kept: ", length(genes_to_keep),"from",nrow(spatial_object), "\n") 
  spatial_object_subset <- SCTransform(spatial_object_subset, assay = "Spatial", verbose = T)
  
})

####obtain the interaction set of genes in the expression matrix####
#take the samples used in RI model as an example
intersection <- intersect(rownames(spatial_list$c_3.h5),rownames(spatial_list$c_4.h5))
for (i in slide_list){
  intersection <- intersect(rownames(spatial_list[[i]]),intersection)
}
for (x in slide_list){
  spatial_list[[x]] <- spatial_list[[x]][intersection,]
}

####Add annotation####
#RI model
slide_list <- c("c_3.h5",
                "c_4.h5",
                "c_36.h5",
                "c_51.h5")
csv_list <- c("c_3.csv",
              "c_4.csv",
              "c_36.csv",
              "c_51.csv")
#NRI model
slide_list <- c("b_1.h5",
                "b_18.h5",
                "a_3.h5",
                "a_15.h5")
csv_list <- c("b_1.csv",
              "b_18.csv",
              "a_3.csv",
              "a_15.csv")
for(slide in slide_list){
  spatial_list[[slide]]$orig.ident <- slide
  for(x in csv_list){
    setwd("C:\\Users\\HUAWEI\\Documents\\Differntial Expression\\data\\flitered_matiexs_ffpe\\annotation_ffpe/")
    #the folder of annotation csv, use the data in 'annotation' folder
    annot_table <- read.csv(x)
    rownames(annot_table) <- annot_table[,1]
    annot_table$Barcode <- NULL
    spatial_list[[slide]] <- AddMetaData(object= spatial_list[[slide]],
                                         metadata = annot_table,
                                         col.name = paste0(colnames(annot_table),"_annot"))
  }
}

####Integrate the training dataset, takeRI model as an example####
train_list <- c(spatial_list[['c_3.h5']],spatial_list[['c_4.h5']],spatial_list[['c_36.h5']])
features <- SelectIntegrationFeatures(object.list = train_list,nfeatures = 14214)
train.anchors <- FindIntegrationAnchors(object.list = train_list, anchor.features = features)
train.combined <- IntegrateData(anchorset = train.anchors)

####Identify the DEG####
TLS_list_all_samples <- FindMarkers(train.combined,group.by = "TLS_2_cat_annot",ident.1 = "TLS",test.use = "MAST",logfc.threshold = 0.01)

####Figure 2: the Volcano Plot####
library(ggplot2)
library(ggrepel)
df <- TLS_list_all_samples
# 注：pthreshold和fcthreshold为自己设置的p值和fc值阈值
fcthreshold <- 2
pthreshold <- 0.05
# 设置点的分类
df$change <- as.factor(ifelse(df$p_val_adj<pthreshold & abs(df$avg_log2FC)>log2(fcthreshold),
                              ifelse(df$avg_log2FC>log2(fcthreshold),"Up","Down"),"Non"))
rownames(df)<-c
# 样本标签
df$label <- ifelse(df$p_val_adj<pthreshold & abs(df$avg_log2FC)>log2(fcthreshold),as.character(rownames(df)),"")
# 绘制火山图
p.vol <- ggplot(data = df,
                aes(x = avg_log2FC,y = -log10(p_val_adj),colour = change,fill = change))+
  scale_color_manual(values = c('green','grey','red'))+
  geom_point(alpha = 0.4,size = 3.5)+
  # 标签
  geom_text_repel(aes(x = avg_log2FC,y = -log10(p_val_adj),label = label),size = 3,
                  box.padding = unit(0.6,"lines"),point.padding = unit(0.7,"lines"),
                  segment.color = "black",show.legend = FALSE)+
  # 辅助线
  geom_vline(xintercept = c(-log2(fcthreshold),log2(fcthreshold)),lty = 4,col = "black",lwd = 0.8)+
  geom_hline(yintercept = -log10(pthreshold),lty = 4,col = "black",lwd = 0.8)+
  theme_bw()+
  labs(x = "log2(FoldChange)",y = "-log10(p-adjust)",title = "Volcano Plot of  Different Expression Proteins")+
  # 坐标轴标题、标签和图例相关设置
  theme(axis.text = element_text(size = 11),axis.title = element_text(size = 13), # 坐标轴标签和标题
        plot.title = element_text(hjust = 0.5,size = 15,face = "bold"), # 标题
        legend.text = element_text(size = 11),legend.title = element_text(size = 13), # 图例标签和标题
        plot.margin = unit(c(0.5,0.5,0.5,0.5),"cm")) # 图边距
setwd("C:\\Users\\HUAWEI\\Desktop/")
ggsave(p.vol,filename = paste("Volcano.pdf"))