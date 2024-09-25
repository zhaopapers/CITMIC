#Dumbbell diagram
pancancer_AUC<-read.csv("AUC_Adenocarcinoma.csv")
pancancer_AUC<-cbind(pancancer_AUC[which(pancancer_AUC[,4]=="all"),1:3],pancancer_AUC[which(pancancer_AUC[,4]!="all"),2])
colnames(pancancer_AUC)[4]<-"high_auc"
library(ggplot2)
library(ggalt)
library(ggprism)

ssgsea_auc<-c(0.6335379,0.6465963,0.6354688,0.6838575,0.6405544)
BRCA_AUC<-pancancer_AUC[which(pancancer_AUC[,1]=="BRCA"),]
p1 <- ggplot(BRCA_AUC) +
  geom_segment(aes(x=year, xend=year, y=auc, yend=iii_auc), color="grey",size=1) +
  geom_point( aes(x=year, y=auc), color='#005eaa', size=4 ) +
  geom_point( aes(x=year, y=iii_auc), color='#c12e34', size=4 ) +
  geom_line(aes(x = as.numeric(year), y = ssgsea_auc), color = "black", size = 1.5)+
  theme_prism(palette = "pearl", 
              base_fontface = "plain",
              base_family = "serif", 
              base_size = 14, 
              base_line_size = 0.8,
            ) +ylim(0,1)+
  theme(legend.position = "none") + 
  xlab("Year") +
  ylab("AUC") +
  ggtitle("BRCA")
#LUAD
ssgsea_auc<-c(0.6416234,0.5753193,0.5825402,0.5916086,0.6041333)
test_data_LUAD<-cbind(test_data_LUAD[1:5,],ssgsea_auc)
test_data_LUAD[,5]<-as.numeric(test_data_LUAD[,5])
colnames(test_data_LUAD)[4]<-c("iii_auc")
p2 <- ggplot(test_data_LUAD) 
  geom_segment(aes(x=year, xend=year, y=auc, yend=iii_auc), color="grey",size=1) +
  geom_point( aes(x=year, y=auc), color='#005eaa', size=4 ) +
  geom_point( aes(x=year, y=iii_auc), color='#c12e34', size=4 ) +
  geom_line(aes(x = as.numeric(year), y = ssgsea_auc), color = "black", size = 1.5, linetype = "dashed")+
  theme_prism(palette = "pearl",  
              base_fontface = "plain",
              base_family = "serif", 
              base_size = 14, 
              base_line_size = 0.8,
  ) +ylim(0,1)+
  theme(legend.position = "none") + 
  xlab("Year") +
  ylab("AUC") +
  ggtitle("LUAD")

#READ
ssgsea_auc<-c(0.6703816,0.6578412,0.6674999,0.4769310,0.6525333)
test_data_READ<-cbind(test_data_READ[1:5,],ssgsea_auc)
test_data_READ[,5]<-as.numeric(test_data_READ[,5])
colnames(test_data_READ)[4]<-c("iii_auc")
p3 <- ggplot(test_data_READ) +
  geom_segment(aes(x=year, xend=year, y=auc, yend=iii_auc), color="grey",size=1) +
  geom_point( aes(x=year, y=auc), color='#005eaa', size=4 ) +
  geom_point( aes(x=year, y=iii_auc), color='#c12e34', size=4 ) +
  geom_line(aes(x = as.numeric(year), y = ssgsea_auc), color = "black", size = 1.5, linetype = "dashed")+
  theme_prism(palette = "pearl",  
              base_fontface = "plain",
              base_family = "serif", 
              base_size = 14, 
              base_line_size = 0.8,
  ) +ylim(0,1)+
  theme(legend.position = "none") + 
  xlab("Year") +
  ylab("AUC") +
  ggtitle("READ")

#STAD
test_data_STAD<-new[which(new[,1]=="STAD"),]
ssgsea_auc<-c(0.7149579,0.6062693,0.6559163,0.5692674,0.6998115)
test_data_STAD<-cbind(test_data_STAD[1:5,],ssgsea_auc)
test_data_STAD[,5]<-as.numeric(test_data_STAD[,5])
colnames(test_data_STAD)[4]<-c("iii_auc")
p4 <- ggplot(test_data_STAD) +
  geom_segment(aes(x=year, xend=year, y=auc, yend=iii_auc), color="grey",size=1) +
  geom_point( aes(x=year, y=auc), color='#005eaa', size=4 ) +
  geom_point( aes(x=year, y=iii_auc), color='#c12e34', size=4 ) +
  geom_line(aes(x = as.numeric(year), y = ssgsea_auc), color = "black", size = 1.5, linetype = "dashed")+
  theme_prism(palette = "pearl",  
              base_fontface = "plain",
              base_family = "serif", 
              base_size = 14, 
              base_line_size = 0.8,
  )  +ylim(0,1)+
  theme(legend.position = "none") + 
  xlab("Year") +
  ylab("AUC") +
  ggtitle("STAD")

library(patchwork)
(p1+p2)/(p3+p4)

 
#BRCA
library(TCGAbiolinks)
query<-GDCquery(
  project = "TCGA-BRCA",
  data.category = "Simple Nucleotide Variation",
  data.type = "Masked Somatic Mutation",
  access="open"
  
)

GDCdownload(query)
GDCprepare(query, save=T, save.filename="TCGA-BRCA_SNP.Rdata")


 
library(CITMIC)
library(parallel)
load("BRCA_exp.Rdata")#BRCA exp Standardized methodology consistent with SKCM
lnScore_BRCA<-CITMIC(BRCA_exp,cl.cores=8)
survival<-read.delim("BRCA_survival.txt",header=T,sep = "\t")
survival<-survival[,c("sample","OS","OS.time")]
survival<-survival[which(survival$OS.time!=""),]
survival$sample<-paste0(survival$sample,"A")
survival$sample<-gsub("-",".",survival$sample)

stage<-read.delim("TCGA.BRCA.sampleMap_BRCA_clinicalMatrix",header=T)
stage_all<-read.delim("TCGA.BRCA.subtype.sampleMap_BRCA_clinicalMatrix",header=T)
stage<-stage[,c("sampleID","pathologic_stage")]
stage<-stage[which(stage$pathologic_stage!=""),]
stage$sampleID<-paste0(stage$sampleID,"A")
stage$sampleID<-gsub("-",".",stage$sampleID)

status <- stage_all[,c(8,10,14,16,20,22,92)]
rownames(status) <- stage_all[,1]
rownames(status)<-paste0(rownames(status),"A")
rownames(status)<-gsub("-",".",rownames(status))

var_stage_early<-c("Stage I","Stage II")
stage_I_II<-stage[which(stage$pathologic_stage%in%var_stage_early),1]

survival<-merge(survival,stage,by.x='sample',by.y='sampleID')
survival<-survival[,c("sample","OS","OS.time")]


cell_interact_survival<-cbind(survival[which(survival[,1]%in%intersect(survival[,1],colnames(lnScore_BRCA))),],t(lnScore_BRCA)[intersect(survival[,1],colnames(lnScore_BRCA)),])
cell_interact_survival_unI_II<-cell_interact_survival[-which(cell_interact_survival[,1]%in%stage_I_II),]

single_cox_cell_interact_late<-factor_sing_multi(cell_interact_survival_unI_II)
 
single_cox_cell_interact_late_new<-single_cox_cell_interact_late[,c(1,9)]
BRCA.BQ<-single_cox_cell_interact_late_new
library(maftools)
load("TCGA-BRCA_SNP.Rdata")
maf.BRCA<-data
maf<-read.maf(maf.BRCA)
data<-maf@data
PIK3CA_mut<-data.frame(data[which(data[,1]=="TP53"),])
PIK3CA_mut_sam<-substr(PIK3CA_mut$Tumor_Sample_Barcode,1,16)
PIK3CA_mut_sam<-gsub("-",".",PIK3CA_mut_sam)

FAT_mut<-data.frame(data[which(data[,2]=="FAT3"),])
FAT_mut_sam<-substr(FAT_mut$Tumor_Sample_Barcode,1,16)
FAT_mut_sam<-gsub("-",".",FAT_mut_sam)

BRCA.BQ$TP53 <- 0
BRCA.BQ$FAT3 <- 0
BRCA.BQ[which(rownames(BRCA.BQ)%in%PIK3CA_mut_sam),3]<-"mutant"
BRCA.BQ[-which(rownames(BRCA.BQ)%in%PIK3CA_mut_sam),3]<-"wt"
BRCA.BQ[which(rownames(BRCA.BQ)%in%PIK3CA_mut_sam1),4]<-"mutant"
BRCA.BQ[-which(rownames(BRCA.BQ)%in%PIK3CA_mut_sam1),4]<-"wt"

library(estimate)
write.table(exp_BRCA_fpkm_tumor_aggregate_log2_20000,file="BRCA_exp.txt",quote=F,sep="\t")
filterCommonGenes(input.f="BRCA_exp.txt", output.f="BRCA_exp.gct", id="GeneSymbol")
estimateScore("BRCA_exp.gct", "BRCA_exp_estimate_score.gct", platform="illumina")
BRCA_score<-read.delim("BRCA_exp_estimate_score.gct",fill=TRUE)
BRCA_score <- BRCA_score[-1,]
BRCA_score <- BRCA_score[,-2]
BRCA_score <- t(BRCA_score)
colnames(BRCA_score) <- BRCA_score[1,]
BRCA_score <- BRCA_score[-1,]
rownames(BRCA_score) <- BRCA_score[,1]
BRCA_score <- BRCA_score[,-1]
BRCA.BQ <- merge(BRCA.BQ,BRCA_score,by.x = 0, by.y = 0)
rownames(BRCA.BQ) <- BRCA.BQ[,1]
BRCA.BQ <- BRCA.BQ[,-c(1,2)]

#stem score 
library(xlsx)
TCGA_score_Stemness<-read.xlsx("TCGA_score.xlsx",sheetIndex  =1)
TCGA_score_Stemness<-TCGA_score_Stemness[which(TCGA_score_Stemness[,2]=="BRCA"),c(1,4)]
TCGA_score_Stemness[,1]<-substr(TCGA_score_Stemness[,1],1,16)
TCGA_score_Stemness[,1]<-gsub("-",".",TCGA_score_Stemness[,1])
rownames(TCGA_score_Stemness)<-TCGA_score_Stemness[,1]
BRCA.BQ<-cbind(BRCA.BQ,stemness_Score=TCGA_score_Stemness[rownames(BRCA.BQ),2])
BRCA.BQ <- BRCA.BQ[,-c(1,2)]
BRCA.BQ <- merge(BRCA.BQ,status,by.x = 0,by.y = 0,all.x=TRUE)
rownames(BRCA.BQ) <- BRCA.BQ[,1]
BRCA.BQ <- BRCA.BQ[,-1]
colnames(BRCA.BQ) <- c("CTTMEScore","TP53","FAT3","StromalScore","ImmuneScore","ESTIMATEScore","stemness_Score","ER","HER2","sample_type","Node","Subtype","PR","histological_type")
BRCA.BQ[which(BRCA.BQ$histological_type=="Infiltrating Lobular Carcinoma"),14] <- "ILC"
BRCA.BQ[which(BRCA.BQ$histological_type=="Infiltrating Ductal Carcinoma"),14] <- "IDC"
BRCA.BQ[which(BRCA.BQ$histological_type=="Mixed Histology (please specify)"),14] <- "Mixed"
BRCA.BQ[which(BRCA.BQ$histological_type=="Other, specify"),14] <- "Other"
BRCA.BQ[which(BRCA.BQ$histological_type=="Metaplastic Carcinoma"),14] <- "Other"
BRCA.BQ[which(BRCA.BQ$histological_type=="Mucinous Carcinoma"),14] <- "Other"
BRCA.BQ[which(BRCA.BQ$histological_type=="Medullary Carcinoma"),14] <- "Other"
BRCA.BQ[,4] <- as.numeric(BRCA.BQ[,4])
BRCA.BQ[,5] <- as.numeric(BRCA.BQ[,5])
BRCA.BQ[,6] <- as.numeric(BRCA.BQ[,6])
BRCA.BQ[,7] <- as.numeric(BRCA.BQ[,7])
BRCA.BQ[which(is.na(BRCA.BQ[,7])),7] <- ""

library(RColorBrewer)
colors <- list(
  TP53 = c("wt"="grey", "mutant"="#516b91","white"),
  Node = c("Negative"="grey", "Positive"="#516b91","white"),
  PR = c("Negative"="grey", "Positive"="#516b91","white","Indeterminate"="white"),
  ER = c("Negative"="grey", "Positive"="#516b91","Indeterminate"="white","white"),
  HER2 = c("Negative"="grey", "Positive"="#516b91","white","Equivocal"="white"),
  sample_type = c("Negative"="grey", "Positive"="#516b91","white"),
 #PIK3CA = c("wt"="grey", "mutant"="#516b91","white"),
  #KRAS = c("wt"="grey", "mutant"="#516b91","white"),
  #ARID1A = c("wt"="grey", "mutant"="#516b91","white"),
 #RHOA = c("wt"="grey", "mutant"="#516b91","white"),
  FAT3 = c("wt"="grey", "mutant"="#516b91","white"),
  StromalScore = c("#f6efa6","#bf444c"),
 stemness_Score = c("#f6efa6","#bf444c"),
  ImmuneScore = c("#f6efa6","#bf444c"),
  ESTIMATEScore = c("#f6efa6","#bf444c"),
  Lauren.Class = c("Intestinal"="#c12e34", "Diffuse"="#e6b600","Mixed"="#0098d9","NA"="white","white"),
 histological_type= c("ILC"="#d0648a","IDC"="#96dee8","Mixed"="#bebebe","Other"="#fad860"),
 Subtype= c("Basal"="#4ea397","Her2"="#7bd9a5","LumA"="#d0648a","LumB"="#f58db2","Normal"="#ffffff","white"),
 pathologic_stage = c("Stage III"="#008acd","Stage IV"="#d87a80","Stage X"="#ff0000","white"),  
  #MSI.status = c("MSS"="#2b821d", "MSI-H"="#005eaa","MSI-L"="#339ca8","white"),
  #EBV.positive = c("1"="#516b91", "0"="grey","white"),
  #Molecular.Subtype =c("EBV"="#4ea397", "MSI"="#7bd9a5","CIN"="#d0648a","GS"="#f58db2","white"),
 n = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100)
)

library(pheatmap)
min_value <- -2
max_value <- 2

CTTMEScore<-BRCA.BQ$CTTMEScore
names(CTTMEScore)<-rownames(BRCA.BQ)
CTTMEScore<-CTTMEScore[order(CTTMEScore,decreasing = T)]
annotation_col<-BRCA.BQ[,-1]
pheatmap(t(CTTMEScore),cluster_rows = F,cluster_cols = F,annotation_col =annotation_col,
         annotation_colors = colors,show_colnames = F,show_rownames = F)



library(ggpubr)
ggscatter(BRCA.BQ, x = "CTTMEScore", y = "ImmuneScore", 
          color = "black",fill = "lightgray",
          add = "reg.line", conf.int = TRUE, 
          add.params = list( color = "#D0648A",fill = "lightgray",fill = "lightgray"),
          cor.coef = T,
          cor.coeff.args = list(),
          cor.method = "pearson",
          cor.coef.coord = c(-83,-75),
          cor.coef.size = 8)



library(ggstatsplot)
ggbetweenstats(data=BRCA.BQ,TP53,CTTMEScore)
ggbetweenstats(data=BRCA.BQ,histological_type,CTTMEScore)


#STAD
library(TCGAbiolinks)
query<-GDCquery(
  project = "TCGA-STAD",
  data.category = "Simple Nucleotide Variation",
  data.type = "Masked Somatic Mutation",
  access="open"
  
)

GDCdownload(query)
GDCprepare(query, save=TRUE, save.filename="TCGA-STAD_SNP.Rdata")

library(CITMIC)
library(parallel)
load("STAD_exp.Rdata")#STAD exp Standardized methodology consistent with SKCM
lnScore_STAD<-CITMIC(exp_fpkm_tumor_aggregate_log2_20000,cl.cores=8)
survival<-read.delim("STAD_survival.txt",header=T,sep = "\t")
survival<-survival[,c("sample","OS","OS.time")]
survival<-survival[which(survival$OS.time!=""),]
survival$sample<-paste0(survival$sample,"A")
survival$sample<-gsub("-",".",survival$sample)

stage<-read.delim("TCGA.STAD.sampleMap_STAD_clinicalMatrix",header=T)
sub_type<-read.delim("STAD_subtype.csv",header=T,sep = ",",row.names=1)
stage<-stage[,c("sampleID","pathologic_stage")]
stage<-stage[which(stage$pathologic_stage!=""),]
stage$sampleID<-paste0(stage$sampleID,"A")
stage$sampleID<-gsub("-",".",stage$sampleID)

status <- sub_type[,c(3,22,28)]
rownames(status) <- sub_type[,1]
rownames(status)<-paste0(rownames(status),"-01A")
rownames(status)<-gsub("-",".",rownames(status))


var_stage_early<-c("Stage I","Stage II")
stage_I_II<-stage[which(stage$pathologic_stage%in%var_stage_early),1]

survival<-merge(survival,stage,by.x='sample',by.y='sampleID')
survival<-survival[,c("sample","OS","OS.time")]


cell_interact_survival<-cbind(survival[which(survival[,1]%in%intersect(survival[,1],colnames(lnScore_STAD))),],t(lnScore_STAD)[intersect(survival[,1],colnames(lnScore_STAD)),])
cell_interact_survival_unI_II<-cell_interact_survival[-which(cell_interact_survival[,1]%in%stage_I_II),]

single_cox_cell_interact_late<-factor_sing_multi(cell_interact_survival_unI_II)
 
single_cox_cell_interact_late_new<-single_cox_cell_interact_late[,c(1,7)]
STAD.BQ<-single_cox_cell_interact_late_new
library(maftools)
load("TCGA-STAD_SNP.Rdata")
maf.BRCA<-data
maf<-read.maf(maf.BRCA)
data<-maf@data
TP53_mut<-data.frame(data[which(data[,1]=="TP53"),])
TP53_mut_sam<-substr(PIK3CA_mut$Tumor_Sample_Barcode,1,16)
TP53_mut_sam<-gsub("-",".",PIK3CA_mut_sam)

PIK3CA_mut<-data.frame(data[which(data[,2]=="PIK3CA"),])
PIK3CA_mut_sam<-substr(PIK3CA_mut$Tumor_Sample_Barcode,1,16)
PIK3CA_mut_sam<-gsub("-",".",PIK3CA_mut_sam)

RHOA_mut<-data.frame(data[which(data[,2]=="RHOA"),])
RHOA_mut_sam<-substr(RHOA_mut$Tumor_Sample_Barcode,1,16)
RHOA_mut_sam<-gsub("-",".",RHOA_mut_sam)

KRAS_mut<-data.frame(data[which(data[,2]=="KRAS"),])
KRAS_mut_sam<-substr(KRAS_mut$Tumor_Sample_Barcode,1,16)
KRAS_mut_sam<-gsub("-",".",KRAS_mut_sam)

ARID1A_mut<-data.frame(data[which(data[,2]=="ARID1A"),])
ARID1A_mut_sam<-substr(ARID1A_mut$Tumor_Sample_Barcode,1,16)
ARID1A_mut_sam<-gsub("-",".",ARID1A_mut_sam)

STAD.BQ$TP53 <- 0
STAD.BQ$PIK3CA <- 0
STAD.BQ$RHOA <- 0
STAD.BQ$KRAS <- 0
STAD.BQ$ARID1A <- 0
STAD.BQ[which(STAD.BQ[,1]%in%TP53_mut_sam),3]<-"mutant"
STAD.BQ[-which(STAD.BQ[,1]%in%TP53_mut_sam),3]<-"wt"
STAD.BQ[which(STAD.BQ[,1]%in%PIK3CA_mut_sam),4]<-"mutant"
STAD.BQ[-which(STAD.BQ[,1]%in%PIK3CA_mut_sam),4]<-"wt"
STAD.BQ[which(STAD.BQ[,1]%in%RHOA_mut_sam),5]<-"mutant"
STAD.BQ[-which(STAD.BQ[,1]%in%RHOA_mut_sam),5]<-"wt"
STAD.BQ[which(STAD.BQ[,1]%in%KRAS_mut_sam),6]<-"mutant"
STAD.BQ[-which(STAD.BQ[,1]%in%KRAS_mut_sam),6]<-"wt"
STAD.BQ[which(STAD.BQ[,1]%in%ARID1A_mut_sam),7]<-"mutant"
STAD.BQ[-which(STAD.BQ[,1]%in%ARID1A_mut_sam),7]<-"wt"

library(estimate)
write.table(exp_BRCA_fpkm_tumor_aggregate_log2_20000,file="STAD_exp.txt",quote=F,sep="\t")
filterCommonGenes(input.f="STAD_exp.txt", output.f="STAD_exp.gct", id="GeneSymbol")
estimateScore("STAD_exp.gct", "STAD_exp_estimate_score.gct", platform="illumina")
STAD_score<-read.delim("STAD_exp_estimate_score.gct",fill=TRUE)
STAD_score <- STAD_score[-1,]
STAD_score <- STAD_score[,-2]
STAD_score <- t(STAD_score)
colnames(STAD_score) <- STAD_score[1,]
STAD_score <- STAD_score[-1,]
rownames(STAD_score) <- STAD_score[,1]
STAD_score <- STAD_score[,-1]
STAD.BQ <- merge(STAD.BQ,STAD_score,by.x = 0, by.y = 0)
rownames(STAD.BQ) <- STAD.BQ[,1]
STAD.BQ <- STAD.BQ[,-c(1,2)]

#stem score 
library(xlsx)
TCGA_score_Stemness<-read.xlsx("TCGA_score.xlsx",sheetIndex  =1)
TCGA_score_Stemness<-TCGA_score_Stemness[which(TCGA_score_Stemness[,2]=="STAD"),c(1,4)]
TCGA_score_Stemness[,1]<-substr(TCGA_score_Stemness[,1],1,16)
TCGA_score_Stemness[,1]<-gsub("-",".",TCGA_score_Stemness[,1])
rownames(TCGA_score_Stemness)<-TCGA_score_Stemness[,1]
STAD.BQ<-cbind(STAD.BQ,stemness_Score=TCGA_score_Stemness[rownames(STAD.BQ),2])

STAD.BQ <- merge(STAD.BQ,status,by.x = 0,by.y = 0,all.x = TRUE)
STAD.BQ <- merge(STAD.BQ,stage,by.x = 1,by.y = 1,all.x = TRUE)
rownames(STAD.BQ) <- STAD.BQ[,1]
STAD.BQ <- STAD.BQ[,-1]
colnames(STAD.BQ) <- c("CTTMEScore","TP53","PIK3CA","RHOA","KRAS","ARID1A","StromalScore","ImmuneScore","ESTIMATEScore","stemness_Score","Lauren.Class","Molecular.Subtype","MSI.status","EBV.positive","pathologic_stage")
STAD.BQ1$EBV.positive[which(STAD.BQ1$EBV.positive==0)]<-"Negative"
STAD.BQ1$EBV.positive[which(STAD.BQ1$EBV.positive==1)]<-"Positive"
STAD.BQ[,7] <- as.numeric(STAD.BQ[,7])
STAD.BQ[,8] <- as.numeric(STAD.BQ[,8])
STAD.BQ[,9] <- as.numeric(STAD.BQ[,9])
STAD.BQ[,10] <- as.numeric(STAD.BQ[,10])

library(RColorBrewer)
colors <- list(
  TP53 = c("wt"="grey", "mutant"="#516b91","white"),
  #Node = c("Negative"="grey", "Positive"="#516b91","white"),
  #PR = c("Negative"="grey", "Positive"="#516b91","white","Indeterminate"="white"),
  #ER = c("Negative"="grey", "Positive"="#516b91","Indeterminate"="white","white"),
  #HER2 = c("Negative"="grey", "Positive"="#516b91","white","Equivocal"="white"),
  #sample_type = c("Negative"="grey", "Positive"="#516b91","white"),
 PIK3CA = c("wt"="grey", "mutant"="#516b91","white"),
  KRAS = c("wt"="grey", "mutant"="#516b91","white"),
  ARID1A = c("wt"="grey", "mutant"="#516b91","white"),
 RHOA = c("wt"="grey", "mutant"="#516b91","white"),
  #FAT3 = c("wt"="grey", "mutant"="#516b91","white"),
  StromalScore = c("#f6efa6","#bf444c"),
 stemness_Score = c("#f6efa6","#bf444c"),
  ImmuneScore = c("#f6efa6","#bf444c"),
  ESTIMATEScore = c("#f6efa6","#bf444c"),
  Lauren.Class = c("Intestinal"="#c12e34", "Diffuse"="#e6b600","Mixed"="#0098d9","NA"="white","white"),
 #histological_type= c("ILC"="#d0648a","IDC"="#96dee8","Mixed"="#bebebe","Other"="#fad860"),
 Subtype= c("Basal"="#4ea397","Her2"="#7bd9a5","LumA"="#d0648a","LumB"="#f58db2","Normal"="#ffffff","white"),
 pathologic_stage = c("Stage III"="#008acd","Stage IV"="#d87a80","Stage X"="#ff0000","white"),  
  MSI.status = c("MSS"="#2b821d", "MSI-H"="#005eaa","MSI-L"="#339ca8","white"),
  EBV.positive = c("Positive"="#516b91", "Negative"="grey","white"),
  Molecular.Subtype =c("EBV"="#4ea397", "MSI"="#7bd9a5","CIN"="#d0648a","GS"="#f58db2","white"),
 n = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100)
)

library(pheatmap)
min_value <- -2
max_value <- 2
CTTMEScore<-STAD.BQ$CTTMEScore
names(CTTMEScore)<-rownames(STAD.BQ)
CTTMEScore<-CTTMEScore[order(CTTMEScore,decreasing = T)]
annotation_col<-STAD.BQ[,-1]
pheatmap(t(CTTMEScore),cluster_rows = F,cluster_cols = F,annotation_col =annotation_col,
         annotation_colors = colors,show_colnames = F,show_rownames = F)


#cor plot
library(ggpubr)
ggscatter(STAD.BQ, x = "CTTMEScore", y = "ImmuneScore", 
          color = "black",fill = "lightgray",
          add = "reg.line", conf.int = TRUE, 
          add.params = list( color = "#D0648A",fill = "lightgray",fill = "lightgray"),
          cor.coef = T,
          cor.coeff.args = list(),
          cor.method = "pearson",
          cor.coef.coord = c(-83,-75),
          cor.coef.size = 8)


#box plot
library(ggstatsplot)
ggbetweenstats(data=STAD.BQ,ARID1A,CTTMEScore)
ggbetweenstats(data=STAD.BQ,Molecular.Subtype,CTTMEScore)
