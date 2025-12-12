A <- prop.table(table(gexpr$subtype_ord,gexpr$time_ord), margin = 2)
A <- as.data.frame(A)
colnames(A) <- c("subtype","age", "Freq")

p4=ggplot(A,aes(x = age,y =Freq,
                group=subtype))+
  stat_summary(geom = 'line',fun='mean',cex=1,col='white')+
  geom_area(data = A,aes(fill=subtype))+
  scale_fill_manual(values=subtypecol)+
  labs(x=NULL,y=NULL)+
  theme_bw()+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        axis.text = element_text(color = "black",size = 10))+
  geom_vline(aes(xintercept ="E55"),linetype="dashed", size=1, colour="white")+
  geom_vline(aes(xintercept ="E66"),linetype="dashed", size=1, colour="white")+
  geom_vline(aes(xintercept ="E76"),linetype="dashed", size=1, colour="white")+
  geom_vline(aes(xintercept ="E85"),linetype="dashed", size=1, colour="white")+
  geom_vline(aes(xintercept ="E94"),linetype="dashed", size=1, colour="white")+
  geom_vline(aes(xintercept ="E104"),linetype="dashed", size=1, colour="white")+
  geom_vline(aes(xintercept ="E109"),linetype="dashed", size=1, colour="white")










library(MetaNeighbor)
refdata=Adult
refsce=as.SingleCellExperiment(refdata)
refsce$study_id <- "ref"
sce=as.SingleCellExperiment(gexpr1)
sce$study_id="test"
hvgs = variableGenes(refsce, exp_labels = refsce$study_id)[1:1000]
pretrained_model <- MetaNeighbor::trainModel(var_genes = unique(c(hvgs,recMarkers$gene)),
                                             dat = refsce,
                                             study_id = refsce$study_id, 
                                             cell_type = refsce$subclass 
)

aurocs <- MetaNeighborUS(trained_model=pretrained_model,
                         dat = sce,
                         study_id = sce$study_id,
                         cell_type = sce$subtype,
                         fast_version=TRUE)
colnames(aurocs)=sapply(colnames(aurocs),function(x) strsplit(x,fixed = T,split = "|")[[1]][2])
rownames(aurocs)=sapply(rownames(aurocs),function(x) strsplit(x,fixed = T,split = "|")[[1]][2])
plotHeatmapPretrained(t(aurocs),margins = c(15,15))





pretrained_model <- MetaNeighbor::trainModel(var_genes = unique(c(hvgs,recMarkers$gene)),
                                             dat = refsce,
                                             study_id = refsce$study_id, 
                                             cell_type = refsce$subtype 
)

aurocs <- MetaNeighborUS(trained_model=pretrained_model,
                         dat = sce,
                         study_id = sce$study_id,
                         cell_type = sce$subtype,
                         fast_version=TRUE)
colnames(aurocs)=sapply(colnames(aurocs),function(x) strsplit(x,fixed = T,split = "|")[[1]][2])
rownames(aurocs)=sapply(rownames(aurocs),function(x) strsplit(x,fixed = T,split = "|")[[1]][2])
plotHeatmapPretrained(t(aurocs),margins = c(15,15))







library(AUCell)
exprMatrix <- as.matrix(GetAssayData(gexpr, slot = "data"))
geneSets <- lapply(term_names, function(genes) {
  intersect(genes, rownames(exprMatrix))
})


cells_rankings <- AUCell_buildRankings(exprMatrix, plotStats = FALSE)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings,nCores =10)


auc_matrix <- getAUC(cells_AUC)
gexpr[["AUC"]] <- CreateAssayObject(data = auc_matrix)


gene_cell_exp <- AverageExpression(gexpr, 
                                   features = names(term_names),
                                   assays = "AUC",
                                   group.by = c("subtype_ord"),  # 元数据中的时间字段
                                   slot = "data")[[1]]


gene_cell_exp <- t(scale(t(gene_cell_exp),scale = T,center = T))

p=Heatmap(as.matrix(gene_cell_exp),
          name = "Expression",
          col = colorRampPalette(c("#a0d8ef","white","#d9333f"))(100),
          na_col = "grey90",
          cluster_columns = F,
          cluster_rows = F,
          row_gap = unit(2, "mm"),
          row_title_rot = 0,
          width = unit(ncol(gene_cell_exp)*0.5, "cm"),  
          height = unit(nrow(gene_cell_exp)*0.5, "cm"),
          
          column_gap = unit(0.5, "mm"),
          show_column_names = TRUE,
          
          rect_gp = gpar(col = "white", lwd = 0.5),
          border = "black",
          heatmap_legend_param = list(
            title_position = "leftcenter-rot",
            legend_width = unit(6, "cm"))
)





gexpr_merge$age_subtype=paste(gexpr_merge$age,gexpr_merge$subtype_ord,sep = ".")
pseudo_expr <- AverageExpression(gexpr_merge,
                                 group.by = "age_subtype",
                                 assays = "RNA",
                                 features = unique(c("PVALB","SST","NKX2-1")),
                                 slot = "data",
                                 return.seurat = FALSE)$RNA


expr_long <- as.data.frame(as.table(pseudo_expr)) %>%
  dplyr::rename(gene = Var1, sample_celltype = Var2, expression = Freq) %>%
  tidyr::separate(sample_celltype, into = c("age", "subtype"), sep = "\\.") 

plot_gene_gam_all <- function(gene_name) {
  dat <- expr_long %>% filter(gene == gene_name)
  
  
  
  all_subtypes <- names(subtypecol)[sort(match(unique(dat$subtype),names(subtypecol)))]
  dat <- dat %>% 
    mutate(subtype = factor(subtype, levels = all_subtypes))
  
  dat=dat[dat$expression!=0,]
  ggplot(dat, aes(x = age, y = expression, group = subtype)) +
    
    geom_line(aes(color = subtype), linewidth = 0.6) +
    
    geom_point(aes(color = subtype), size = 1.5) +
    
    scale_color_manual(
      values = subtypecol,
      drop = FALSE  
    ) +
    scale_linetype_manual(
      values = c(2,1)
    ) +
    
    labs(title = gene_name) +
    theme_bw() +
    theme(
      aspect.ratio = 1/1,
      # legend.position = "none", 
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
}



plots <- purrr::map(unique(c("PVALB","SST","NKX2-1")), ~plot_gene_gam_all(.x) )
print(plot_grid(plotlist = plots[1:9], ncol = 3))





FeaturePlot(gexpr3,features =c("PVALB","SST","NKX2-1") ,
            max.cutoff = 4,min.cutoff = 0,ncol = 3 )&
  theme(aspect.ratio = 1/1,axis.line = element_blank(),axis.title = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),legend.key.width=unit(3,'mm'),plot.title = element_text(size = 8))&
  scale_color_gradientn(colours = c("lightgrey" ,"#e41a1c"))



p1=DimPlot(gexpr3,group.by = "subtype_ord",label = T,cols = subtypecol)+
  theme(legend.position = "right",panel.border = element_blank(),aspect.ratio = 1/1,
        axis.title = element_blank(),  #轴标题
        axis.text = element_blank(), # 文本
        axis.ticks = element_blank(),axis.line = element_blank(),plot.title = element_blank())

p2=DimPlot(gexpr3,group.by = "region",label = T,cols = region_color,raster = T)+
  theme(legend.position = "right",panel.border = element_blank(),aspect.ratio = 1/1,
        axis.title = element_blank(),  #轴标题
        axis.text = element_blank(), # 文本
        axis.ticks = element_blank(),axis.line = element_blank(),plot.title = element_blank())








gexpr=subset(STR,subclass%in% c("mgeInN"))
gexpr2=readRDS("result/PFC/second_PFC_InN/PFC_first_second_InN_rmlowQ.rds")
gexpr2=subset(gexpr2,subclass %in% c("mgeInN"))
dim(gexpr2)

aa=subset(gexpr,PVALB>0)
bb=subset(gexpr2,PVALB>0)

mm=as.data.frame(table(bb$time))
mm$region="PFC"

nn=as.data.frame(table(aa$time))
nn$region="STR"

mmnn=rbind(mm,nn)
mmnn$Var1=factor(mmnn$Var1,levels = names(age_color)[sort(match(unique(mmnn$Var1),names(age_color)))])


p3=ggplot(mmnn, aes(x=region, y=Freq,fill=Var1)) + 
  geom_bar(stat="identity", position = "fill",width = 0.9)+
  scale_fill_manual(values = age_color)+
  xlab("")+
  ylab("")+
  # ggtitle(species.use)+
  theme(panel.background=element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.ticks.x = element_blank(),
        aspect.ratio = 1/1
  )+
  theme(legend.position = "right",legend.title=element_text(size=12),legend.text=element_text(size=8))





mydata <- FetchData(aa,vars = c("subtype"))
mydata <- mydata %>% group_by(subtype) %>% summarise(counts=n())
mydata$Freq <- round(mydata$counts/sum(mydata$counts),4)

library(scales)
p1=ggplot(mydata,aes(x="",y=Freq,fill=subtype))+geom_bar(stat = "identity")+
  scale_fill_manual(values =subtypecol)+coord_polar(theta = "y")+
  theme(axis.text = element_blank(),axis.ticks = element_blank(),panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5,face = "bold"),panel.background = element_blank(),
  )+
  geom_text(aes(y=cumsum(rev(Freq))-rev(Freq)/2,x=1.5,label=rev(percent(Freq))),size=2)




mydata <- FetchData(bb,vars = c("subtype"))
mydata <- mydata %>% group_by(subtype) %>% summarise(counts=n())
mydata$Freq <- round(mydata$counts/sum(mydata$counts),4)
mydata$subtype=factor(mydata$subtype,levels = levels(gexpr2$subtype_ord)[sort(match(unique(mydata$subtype),levels(gexpr2$subtype_ord)))])

p2=ggplot(mydata,aes(x="",y=Freq,fill=subtype))+geom_bar(stat = "identity")+
  scale_fill_manual(values =subtypecol)+coord_polar(theta = "y")+
  theme(axis.text = element_blank(),axis.ticks = element_blank(),panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5,face = "bold"),panel.background = element_blank(),
  )+
  geom_text(aes(y=cumsum(rev(Freq))-rev(Freq)/2,x=1.5,label=rev(percent(Freq))),size=2)







gene_cell_exp <- AverageExpression(gexpr, 
                                   features = rownames(auc_matrix),
                                   assays = "AUC",
                                   group.by = c("time_ord","subclass_ord"),  
                                   slot = "data")[[1]]

gene_cell_exp <- t(scale(t(gene_cell_exp),scale = T,center = T))
table(gene_cell_exp>2)
table(gene_cell_exp< -2)
gene_cell_exp[gene_cell_exp > 2]=2
gene_cell_exp[gene_cell_exp < -2]=-2

gene_cell_exp=as.data.frame(gene_cell_exp) %>% rownames_to_column(var = "ion")

gene_cell_exp=gene_cell_exp %>% pivot_longer(cols = !ion,names_to = "Age_subclass", values_to = "expression") %>% 
  tidyr::separate(col =Age_subclass,into = c("age","subclass"),sep = "_" )


gene_cell_exp <- gene_cell_exp %>%
  mutate(ion_subclass = paste(ion, subclass, sep = "_")) %>%
  select(-ion, -subclass)


gene_cell_exp=gene_cell_exp %>% pivot_wider(id_cols=ion_subclass,names_from = age, values_from = expression) %>% 
  column_to_rownames(var = "ion_subclass")


row_meta <- data.frame(full_name = rownames(gene_cell_exp)) %>%
  tidyr::separate(full_name, 
                  into = c("ion", "subclass"),
                  sep = "_", 
                  extra = "merge")



library(ComplexHeatmap)
library(circlize)
p3=Heatmap(as.matrix(gene_cell_exp),
           name = "Expression",
           col = colorRampPalette(c("#a0d8ef","white","#d9333f"))(100),
           na_col = "grey90",
           cluster_columns = F,
           cluster_rows = F,
           
           
           row_split = factor(row_meta$ion, 
                              levels = c("CL", "Na", "Glu", "GABA", "K", "Ca","Kca")),
           row_gap = unit(2, "mm"),
           row_title_rot = 0,
           width = unit(ncol(gene_cell_exp)*0.5, "cm"),  
           height = unit(nrow(gene_cell_exp)*0.5, "cm"),
           
           column_gap = unit(0.5, "mm"),
           show_column_names = TRUE,
           
           rect_gp = gpar(col = "white", lwd = 0.5),
           border = "black",
           heatmap_legend_param = list(
             title_position = "leftcenter-rot",
             legend_width = unit(6, "cm"))
)






gene_cell_exp <- AverageExpression(gexpr,
                                   features = rownames(auc_matrix),
                                   group.by = c("time_ord","subtype_ord"),
                                   assays = "AUC",
                                   slot = 'data')
gene_cell_exp <- as.data.frame(gene_cell_exp$AUC) 
gene_cell_exp <- t(scale(t(gene_cell_exp),scale = T,center = T))
table(gene_cell_exp>2)
table(gene_cell_exp< -2)
gene_cell_exp[gene_cell_exp > 2]=2
gene_cell_exp[gene_cell_exp < -2]=-2

gene_cell_exp=as.data.frame(gene_cell_exp) %>% rownames_to_column(var = "ion")

gene_cell_exp=gene_cell_exp %>% pivot_longer(cols = !ion,names_to = "Age_subtype", values_to = "expression") %>% 
  tidyr::separate(col =Age_subtype,into = c("age","subtype"),sep = "_" )
gene_cell_exp_all=gene_cell_exp



gene_cell_exp=gene_cell_exp_all[gene_cell_exp_all$subtype %in% subtype.use,]
gene_cell_exp <- gene_cell_exp %>%
  mutate(ion_subtype = paste(ion, subtype, sep = "_")) %>%
  select(-ion, -subtype)
gene_cell_exp=gene_cell_exp %>% pivot_wider(id_cols=ion_subtype,names_from = age, values_from = expression) %>% 
  column_to_rownames(var = "ion_subtype")


rownames.use=paste(c("CL", "Na", "Glu", "GABA", "K", "Ca","Kca"),rep(subtype.use,each=7),sep="_")
gene_cell_exp=gene_cell_exp[rownames.use,levels(gexpr$time_ord)]


row_meta <- data.frame(full_name = rownames(gene_cell_exp)) %>%
  tidyr::separate(full_name, 
                  into = c("ion", "subtype"),
                  sep = "_", 
                  extra = "merge")

p=Heatmap(as.matrix(gene_cell_exp),
          name = "Expression",
          col = colorRampPalette(c("#a0d8ef","white","#d9333f"))(100),
          na_col = "grey90",
          cluster_columns = F,
          cluster_rows = F,
          
          
          row_split = factor(row_meta$ion, 
                             levels = c("CL", "Na", "Glu", "GABA", "K", "Ca","Kca")),
          
          width = unit(ncol(gene_cell_exp)*0.5, "cm"),  
          height = unit(nrow(gene_cell_exp)*0.5, "cm"), 
          row_gap = unit(0.2, "mm"),
          row_title_rot = 0,
          
          column_gap = unit(0.5, "mm"),
          show_column_names = TRUE,
          
          rect_gp = gpar(col = "white", lwd = 0.5),
          border = "black",
          heatmap_legend_param = list(
            title_position = "leftcenter-rot",
            legend_width = unit(6, "cm"))
          
          
