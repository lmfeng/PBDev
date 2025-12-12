suppressPackageStartupMessages(library(scrabbitr))
suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(miloR))
suppressPackageStartupMessages(library(DelayedArray))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(ggraph))
suppressPackageStartupMessages(library(igraph))
suppressPackageStartupMessages(library(viridis))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(jsonlite))
suppressPackageStartupMessages(library(ggrastr))
suppressPackageStartupMessages(library(ggridges))
suppressPackageStartupMessages(library(ggalluvial))
suppressPackageStartupMessages(library(ggrepel))
library(Seurat)





set.seed(1234)
m_species="Macaque"
r_species="Human"
file=paste0("MF2_related/scrabbitr_Lister/",m_species,"2",r_species,"/")



marmoset=readRDS("5.macaque_pcd127.147_rmlowQ.rds")
DimPlot(marmoset,group.by = "subclass",label = T)
sort(unique(marmoset$subclass))
unique(marmoset$class)
dim(marmoset)

m_seurat=marmoset
m_seurat$species=m_species
m_data=as.SingleCellExperiment(m_seurat)



r_seurat=readRDS("Lister.rds")
r_seurat$species=r_species
r_data=as.SingleCellExperiment(r_seurat)


p=DimPlot(r_seurat,group.by = "class",label = T,label.size = 3,repel = T)+theme(aspect.ratio = 1/1)+labs(title = r_species)
print(p)
# p=DimPlot(r_seurat,group.by = "subclass",label = T,label.size = 3,repel = T)+theme(aspect.ratio = 1/1)+labs(title = r_species)
# print(p)
# p=DimPlot(r_seurat,group.by = "subtype",label = T,label.size = 3,repel = T)+theme(aspect.ratio = 1/1)+labs(title = r_species)
# print(p)
p=DimPlot(m_seurat,group.by = "class",label = T,label.size = 3,repel = T)+theme(aspect.ratio = 1/1)+labs(title = m_species)
print(p)
p=DimPlot(m_seurat,group.by = "subclass",label = T,label.size = 3,repel = T)+theme(aspect.ratio = 1/1)+labs(title = m_species)
print(p)
p=DimPlot(m_seurat,group.by = "subtype",label = T,label.size = 3,repel = T)+theme(aspect.ratio = 1/1)+labs(title = m_species)
print(p)




head(reducedDim(m_data, "PCA"))
class(reducedDim(m_data, "PCA"))
head(reducedDim(r_data, "PCA"))
class(reducedDim(r_data, "PCA"))


r_milo <- Milo(r_data)
r_milo <- buildGraph(r_milo, k=30, d=50, reduced.dim="PCA")
r_milo <- makeNhoods(r_milo, prop=0.05, k=30, d=50,refined=T, reduced_dims="PCA")

r_milo <- buildNhoodGraph(r_milo)
r_milo


p1 <- scrabbitr::plotNhoodSizeHist(r_milo, colour="blue")+theme(aspect.ratio = 1/1)

p2 <- plotNhoodGraph(r_milo, size_range=c(0.1,3), node_stroke=0.1) + 
  scale_fill_viridis(name = "Nhood size", option = "viridis", direction = 1) +theme(aspect.ratio = 1/1)+labs(title = r_species)


grid.arrange(p1, p2, nrow=1,widths=c(1, 3))


m_milo <- Milo(m_data)
m_milo <- buildGraph(m_milo, k=30, d=50, reduced.dim="PCA")
m_milo <- makeNhoods(m_milo, prop=0.05, k=30, d=50,refined=T, reduced_dims="PCA")
m_milo <- makeNhoods(m_milo, prop=0.05, k=30, d=50,refined=T, reduced_dims="PCA")
m_milo <- buildNhoodGraph(m_milo)
m_milo

p1 <- plotNhoodSizeHist(m_milo)+theme(aspect.ratio = 1/1)


p2 <- plotNhoodGraph(m_milo,size_range=c(0.1,3) ,node_stroke=0.1) + 
  scale_fill_viridis(name = "Nhood size", option = "viridis", direction=1)+theme(aspect.ratio = 1/1)+labs(title = m_species)

grid.arrange(p1, p2, nrow=1,widths=c(1, 3))



rm_orthologs <- data.frame(ref=intersect(rownames(m_seurat),rownames(r_seurat)),query=intersect(rownames(m_seurat),rownames(r_seurat)))
##excluding mitochondrial and cell-cycle related genes
mt.genes <- grep(pattern = "^MT-",rownames(m_seurat),ignore.case = T,value = T)
CaseMatch(c(cc.genes$s.genes,cc.genes$g2m.genes),rownames(m_seurat))
g2m_genes = cc.genes$g2m.genes
g2m_genes = CaseMatch(search = g2m_genes, match = rownames(m_seurat))
s_genes = cc.genes$s.genes
s_genes = CaseMatch(search = s_genes, match = rownames(m_seurat))
m_exclude=c(mt.genes,g2m_genes,s_genes)


mt.genes <- grep(pattern = "^MT-",rownames(r_seurat),ignore.case = T,value = T)
CaseMatch(c(cc.genes$s.genes,cc.genes$g2m.genes),rownames(r_seurat))
g2m_genes = cc.genes$g2m.genes
g2m_genes = CaseMatch(search = g2m_genes, match = rownames(r_seurat))
s_genes = cc.genes$s.genes
s_genes = CaseMatch(search = s_genes, match = rownames(r_seurat))
r_exclude=c(mt.genes,g2m_genes,s_genes)


save(r_milo,m_milo,m_data,r_data,m_seurat,r_seurat,file=paste0(file,r_species,"2",m_species,".milo.Rdata"))
out <- scrabbitr::calcNhoodSim(r_milo, m_milo, 
                               orthologs = rm_orthologs, 
                               sim_preprocessing="gene_spec", sim_measure="pearson",
                               hvg_join_type="intersection", max_hvgs=2000, 
                               r_exclude = r_exclude, m_exclude = m_exclude,
                               export_dir = file, 
                               verbose = TRUE)





nhood_sim <- as.matrix(fread(file = paste0(file,"nhood_sim.tsv"), sep="\t"), rownames=1)
dim(nhood_sim)
r_vals <- fread(paste0(file,"r_vals.tsv"), sep="\t")
m_vals <- fread(paste0(file,"m_vals.tsv"), sep="\t")
out <- list(r_vals = r_vals, m_vals = m_vals, nhood_sim = nhood_sim)
save(out,file =paste0(file,"nhood_vals.Rdata"))








library("limma")
library(ggrepel)
library(patchwork)
library(RAPToR)
library(Seurat)
library(dplyr)
library(openxlsx)
library(ggplot2)
library(glmnet)
library(purrr)
library(ComplexHeatmap)
library(dplyr)
library(tidyr)
library(gridExtra)
genes = intersect(rownames(Lister), rownames(gexpr.pig))
Lister_ALL = Lister[genes, ]
gexpr2.pig_ALL = gexpr.pig[genes, ]


for(i in 1:3){
  class.use=c("ExN","InN","Glia")[i]
  cat(class.use,"\t","\n")
  Lister_1=subset(Lister_ALL,class %in% class.use)
  gexpr2.pig=subset(gexpr2.pig_ALL,class %in% class.use)
  unique(Lister_1$class)
  unique(gexpr2.pig$class)
  Idents(Lister_1)="orig.ident"
  unique(Lister_1$orig.ident)
  gene_cell_exp <- AverageExpression(Lister_1,
                                     features = rownames(Lister_1),
                                     group.by = 'orig.ident',
                                     slot = 'data') 
  gene_cell_exp <- as.data.frame(gene_cell_exp$RNA)
  
  meta.hsa=unique(FetchData(Lister_1,vars = c("orig.ident","age","Age_Range","logDays","labsite")))
  rownames(meta.hsa)=NULL
  rownames(meta.hsa)=meta.hsa$orig.ident
  meta.hsa=meta.hsa[colnames(gene_cell_exp),]
  
  
  mat_hsa = ge_im(X = gene_cell_exp, p = meta.hsa,
                  formula = "X ~ s(logDays, bs = 'ts')")
  
  ref_hsa <-make_ref(
    m=mat_hsa,
    n.inter=100, #interpolationresolution
    t.unit ="hpastegg-laying", #timeunit
    metadata=list("organism"="hsa", #anymetadata
                  "profiling"="whole-organism,bulk",
                  "technology"="RNAseq")
  )
  
  saveRDS(ref_hsa,file = paste0("ref_hsa_Lister_",class.use,".rds"))
  
  
  
  
  Idents(gexpr2.pig)="orig.ident"
  # sort(table(gexpr2.pig$age_subclass))
  gene_cell_exp_pig <- AverageExpression(gexpr2.pig,
                                         features = rownames(gexpr2.pig),
                                         group.by = 'orig.ident',
                                         slot = 'data') 
  gene_cell_exp_pig <- as.data.frame(gene_cell_exp_pig$RNA)
  
  
  
  meta.pig_sc=unique(FetchData(gexpr2.pig,vars = c("orig.ident","time","labsite")))
  rownames(meta.pig_sc)=NULL
  rownames(meta.pig_sc)=meta.pig_sc$orig.ident
  meta.pig_sc=meta.pig_sc[colnames(gene_cell_exp_pig),]
  pig_age_predict_1 = ae(gene_cell_exp_pig, ref_hsa)
  
  
  gene_cell_exp <- AverageExpression(ListetKriegstein2023,
                                     features = intersect(rownames(Lister_1),rownames(ListetKriegstein2023)),
                                     group.by = 'orig.ident',
                                     slot = 'data') 
  gene_cell_exp <- as.data.frame(gene_cell_exp$RNA)
  
  
  
  meta.hsa=unique(FetchData(ListetKriegstein2023,vars = c("orig.ident","age","Age_Range","logDays","labsite")))
  rownames(meta.hsa)=NULL
  rownames(meta.hsa)=meta.hsa$orig.ident
  meta.hsa=meta.hsa[colnames(gene_cell_exp),]
  
  dataMat = cbind(gene_cell_exp, gene_cell_exp_pig[rownames(gene_cell_exp),])
  pheno = data.frame(logDays = c(meta.hsa$logDays, pig_age_predict_1$age.estimates[,1]),
                     Age_Range = c(meta.hsa$Age_Range, meta.pig_sc$time),
                     Age = c(meta.hsa$age, meta.pig_sc$time),
                     Species = c(rep("Human", nrow(meta.hsa)), rep("Pig", nrow(meta.pig_sc))),
                     orig.ident=c(meta.hsa$orig.ident, meta.pig_sc$orig.ident),
                     labsite=c(meta.hsa$labsite, meta.pig_sc$labsite)
  )
  idx.sort = order(pheno$logDays, decreasing = F)
  dataMat = dataMat[, idx.sort]
  pheno = pheno[idx.sort, ]
  pheno$Days=round(10^(pheno$logDays),1)
  
  save(dataMat,pheno,file = paste0("dataMat.pheno_Lister_",class.use,".RData"))
  
}





pheno_1=pheno[!(pheno$Age_Range %in% c("3rd trimester")),]
dataMat_1=dataMat[,pheno_1$orig.ident]

colorRef=c("#a9be7b","#007947","#ed1941")
names(colorRef)=c("Listerlab","Kriegstein2023","Lilab")

genesRef=c("DCX",
           "CAMK2A","MAP1A","MAPT",
           "SYN1","SYP","SYPL1","SYPL2")

genesRef=intersect(genesRef,rownames(ListetKriegstein2023))

p_list=list()
for (j in 1:length(genesRef)) {
  gene.use=genesRef[j]
  dcx_expr <- dataMat_1[gene.use, ,drop = FALSE] %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column("orig.ident") %>%
    rename(expression = gene.use)
  
  
  plot_data <- pheno_1 %>%
    inner_join(dcx_expr, by = "orig.ident") %>%
    select(orig.ident, labsite, Days, logDays, expression)
  plot_data$labsite=factor(plot_data$labsite,levels = c("Listerlab","Kriegstein2023","Lilab"))
  
  
  custom_days <- c(100, 200, 500, 2000, 10000, 20000)
  
  custom_logdays <- log10(custom_days)  
  
  
  p=ggplot(plot_data, aes(x = logDays, y = expression)) +
    geom_point(aes(color = labsite),size = 1) +
    geom_smooth(
      method = "loess", 
      formula = y ~ x, 
      se = TRUE, 
      size = 1, 
      alpha = 0.2
    ) +
    geom_vline(
      xintercept = log10(266),  
      color = "black",          
      linetype = "solid",      
      linewidth = 0.7           
    )+
    geom_vline(
      xintercept = log10(c(98,196,645,1010,1740,3930,7580)-14),  
      color = "grey",          
      linetype = "dashed",      
      linewidth = 0.7           
    )+
    
    scale_x_continuous(
      breaks = custom_logdays,       
      labels = custom_days,          
      name = "Days",                
      limits = c(min(custom_logdays)-0.08,max(custom_logdays)) 
    ) +
    labs(
      y = paste0("Expression"),
      title = gene.use,
      color = "Labsite"
    ) +
    scale_color_manual(values = colorRef)+
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 12),
      legend.title = element_text(size = 11),
      panel.grid = element_blank()
    )
  
  p_list[[j]]=p
}





wrap_plots(p_list,guides = "collect",ncol=4)






scRNAlist=SplitObject(gexpr, split.by = "orig.ident")
for (i in 1:length(scRNAlist)) {
  scRNAlist[[i]] <- NormalizeData(scRNAlist[[i]])
  scRNAlist[[i]] <- FindVariableFeatures(scRNAlist[[i]], selection.method = "vst", nfeatures = 2000)
}

features <- SelectIntegrationFeatures(object.list = scRNAlist)


scRNA.anchors <- FindIntegrationAnchors(object.list = scRNAlist,dims = 1:20,anchor.features = features)

scRNA4 <- IntegrateData(anchorset = scRNA.anchors,dims = 1:20)
dim(scRNA4)
testAB.integrated <- scRNA4


DefaultAssay(testAB.integrated) <- "integrated" 
testAB.integrated <- ScaleData(testAB.integrated, features = rownames(testAB.integrated))
testAB.integrated <- RunPCA(testAB.integrated, npcs = 50, verbose = FALSE)  
DimPlot(testAB.integrated, reduction = "pca", group.by="orig.ident")
ElbowPlot(testAB.integrated, ndims=30, reduction="pca") 
testAB.integrated <- FindNeighbors(testAB.integrated, dims = 1:30)
testAB.integrated <- FindClusters(testAB.integrated, resolution = 0.5)
testAB.integrated <- RunUMAP(testAB.integrated, dims = 1:30)
testAB.integrated <- RunTSNE(testAB.integrated, dims = 1:30)
p1<- DimPlot(testAB.integrated,label = T,split.by = 'orig.ident') 

DefaultAssay(testAB.integrated) <- "RNA"
testAB.integrated <- ScaleData(testAB.integrated, features = rownames(testAB.integrated))
testAB.integrated <- RunPCA(testAB.integrated, npcs = 50, verbose = FALSE)
DimPlot(testAB.integrated, reduction = "pca", group.by="orig.ident")
ElbowPlot(testAB.integrated, ndims=30, reduction="pca") 
testAB.integrated <- FindNeighbors(testAB.integrated, dims = 1:30)  
testAB.integrated <- FindClusters(testAB.integrated, resolution = c(1,2,3))
testAB.integrated <- RunUMAP(testAB.integrated, dims = 1:30)
testAB.integrated <- RunTSNE(testAB.integrated, dims = 1:30)
DimPlot(testAB.integrated,label = T,group.by = 'orig.ident')+theme(aspect.ratio = 1/1)
DimPlot(testAB.integrated,label = T,group.by = 'species')+theme(aspect.ratio = 1/1)
DimPlot(testAB.integrated,label = T,group.by = "subclass")+theme(aspect.ratio = 1/1)
DimPlot(testAB.integrated,label = T,group.by = "integrated_snn_res.0.5")+theme(aspect.ratio = 1/1)
DimPlot(testAB.integrated,label = T,group.by = "integrated_snn_res.1")+theme(aspect.ratio = 1/1)
DimPlot(testAB.integrated,label = T,group.by = "integrated_snn_res.2")+theme(aspect.ratio = 1/1)
DimPlot(testAB.integrated,label = T,group.by = "integrated_snn_res.3")+theme(aspect.ratio = 1/1)
DimPlot(testAB.integrated,label = T,group.by = "subclass",split.by = "species")&theme(aspect.ratio = 1/1)
DimPlot(testAB.integrated,label = T,group.by = "integrated_snn_res.1",split.by = "species")&theme(aspect.ratio = 1/1)
saveRDS(testAB.integrated,file = "CCA.fourSpeciesKriegstein2023_2nd.rds")

