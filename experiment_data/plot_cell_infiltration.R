#SDY144 

#CITMIC
GSE52005<-read.csv("GSE52005.csv",row.names = 1)
lnscore_GSE52005<-CITMIC(GSE52005)

#observed
cell_protation_SDY144<-read.csv("cell_protation_SDY144.csv")
colnames(cell_protation_SDY144)<-c("cell","patient","predict")
compare_method_prediction<-merge(result_predction,cell_protation_SDY144,by=c("cell","patient"))
cor.test(as.numeric(compare_method_prediction[,3]),as.numeric(compare_method_prediction[,4]),method="spearman")

#ssGSEA
library(GSVA)
GSE52005_zscore<-t(scale(t(GSE52005)))
SKCM_ssgsea_2000<-gsva(as.matrix(GSE52005_zscore),TMEcell_list,method="ssgsea",kcdf="Gaussian")

#xCell
library(xCell)
xcell_out<-xCellAnalysis(GSE52005)

#CIBERSORT
lm22<-read.delim("LM22.txt",row.names = 1)
cibersort<-my_CIBERSORT(GSE52005,lm22, perm=10, QN=TRUE, cores = 3)

#quanTIseq
tcga_quantiseq<-deconvolute(as.matrix(GSE52005),"quantiseq")

#EPIC
EPIC_out <- EPIC(bulk = GSE52005)
library(MCPcounter)
MCPcounter_out<-MCPcounter.estimate(GSE52005,featuresType=c("affy133P2_probesets","HUGO_symbols","ENTREZ_ID","ENSEMBL_ID")[2],
                                    probesets=read.table("probesets.txt",sep="\t",stringsAsFactors=FALSE,colClasses="character"),
                                    genes=read.table("genes.txt",sep="\t",stringsAsFactors=FALSE,header=TRUE,colClasses="character",check.names=FALSE)
                                    
)

#TIMER  website tools

timer<-read.csv("timer_new.csv",row.names = 1)

#predict
lnscore_GSE52005<-t(lnscore_GSE52005)
result_predction<-c()

for(i in colnames(lnscore_GSE52005)){
  
  x<-cbind(rownames(lnscore_GSE52005),lnscore_GSE52005[,i],i)
  result_predction<-rbind(result_predction,x)
  
}
colnames(result_predction)<-c("cell","predict","patient")

#observed
cell_protation_SDY144<-read.csv("SDY144.csv")
colnames(cell_protation_SDY144)<-c("cell","patient","predict")

#compare
compare_method_prediction<-merge(result_predction,cell_protation_SDY144,by=c("cell","patient"))
cor.test(as.numeric(compare_method_prediction[,3]),as.numeric(compare_method_prediction[,4]),method="spearman")

#single compare
comm_cell<-intersect(unique(result_predction[,3]),unique(cell_protation_SDY144[,1]))
single_result<-c()
for(i in comm_cell){
  x<-result_predction[which(result_predction[,3]==i),]
  y<-cell_protation_SDY144[which(cell_protation_SDY144[,1]==i),]
  result<-merge(x,y,by=c("cell","patient"))
  re<-cor.test(as.numeric(result[,3]),as.numeric(result[,4]),method="spearman")
  re.pval<-re$p.value
  re.cor<-as.numeric(re$estimate)
  single_result<-rbind(single_result,cbind(i,re.cor,re.pval))
}
write.csv(single_result,"result.csv")


#scatter plot
mcor<-read.csv("compare_cell_SDY67.csv",row.names = 1)
mcor<-read.csv("compare_cell_SDY144.csv",row.names = 1)
library(corrplot)
col3 <- colorRampPalette(c("#2b821d", "white", "#c12e34"))
corrplot(as.matrix(mcor),tl.col = 'black',
         method = 'pie', 
         col = col3(100),
         addrect = 2,
         na.label = "-",
         rect.col = 'black',
         rect.lwd = 2)

#bar graph
colnames(compare_method_prediction)<-c("cell","patient","predict","observed")
compare_method_prediction_single<-compare_method_prediction[which(compare_method_prediction$cell=="Activated CD4+ T cells"),]
ggplot(compare_method_prediction_single, aes(x = predict, y = observed))+geom_point(size = 6, color = "black")+
  geom_smooth(method = "lm",aes(x =predict , y = observed),
              se = FALSE,inherit.aes = FALSE,color="#d0648a",linewidth=2)+
  geom_text(aes(x = -19, y = 0.55,label = "R = 0.591"),size = 6)+
  geom_text(aes(x = -19, y = 0.57,label = "p = 0.041"),size = 6)



#single cell infiltration
lnScore_GSE86363<-CITMIC(GSE86363_symbol_log2_20000,cl.cores = 8)
write.csv(lnScore_GSE86363,"cell_infiltration_GSE86363.csv")
cell_infiltration_GSE86363<-read.csv("cell_infiltration_GSE86363.csv",row.names = 1)
GSE86363_harmonized_annotation133A<-read.delim("GSE86363_harmonized_annotation133A.txt")
nu_immune<-c("Kidney","Breast","Skin","Lung","Ovary","Large intestine","Central Nervous System",
  "Prostate","Pre-GC B cells","Ewing's sarcoma","Upper aerodigestive tract","Cervix",
  "Liver","Osteosarcoma","Soft Tissue Sarcoma","Peripheral Nervous System","Cell_Lines")
GSE86363_harmonized_annotation133A<-GSE86363_harmonized_annotation133A[-which(GSE86363_harmonized_annotation133A[,4]%in%nu_immune),]
GSE86363_harmonized_annotation133A<-GSE86363_harmonized_annotation133A[,c(1,4)]
common<-intersect(GSE86363_harmonized_annotation133A[,1],colnames(cell_infiltration_GSE86363))
GSE86363_harmonized_annotation133A<-GSE86363_harmonized_annotation133A[which(GSE86363_harmonized_annotation133A[,1]%in%common),]
                                       
rownames(GSE86363_harmonized_annotation)<-GSE86363_harmonized_annotation133A[,1]
GSE86363_harmonized_annotation<-data.frame(GSE86363_harmonized_annotation133A[,-1])
colnames(GSE86363_harmonized_annotation)<-"group"

GSE86363_harmonized_annotation[,1]<-factor(GSE86363_harmonized_annotation[,1],
                                                levels = c(
  "Activated CD4+ T cells","Activated Memory CD4 T cells","CD4 T cells","White Blood Cells",
  "CD8 T cells","Memory CD4 T cells","Resting Memory CD4 T cells","T cells","T gamma delta",
  "NK cells","Th1","Th2","Canonical CD4 Treg cells", "Plasma cells",
  "B cells","CD40 activated B cells","Germinal center B cells","Memory B cells","Naive B cells",
  "Resting B cells","post Germinal-center B cells",
  "Endothelium","Fibroblast",
   "Monocytes","Neutrophils",
  "Eosinophils","Mast cell",
  "Activated Macrophages","Alternatively Activated Macrophages","Macrophages",
  "Immature Dendritic cells","Mature dendritic cells","Myeloid Dendritic Cells","Plasmacytoid Dendritic Cells"
  ))

my_colors <- colorRampPalette(c("#5394cd","white", "#c12e34"))(100)
                                    
min_value <- -2
max_value <- 2
breaks <- seq(min_value, max_value, length.out = 101)
pheatmap(cell_infiltration_GSE86363[,rownames(GSE86363_harmonized_annotation)[order(GSE86363_harmonized_annotation[,1],decreasing = F)]],
         show_colnames = F,show_rownames = T,cluster_cols = F,cluster_rows = F,breaks=breaks,color=my_colors,
         scale="row",annotation_col=GSE86363_harmonized_annotation)
