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
(p1+p2/(p3+p4)

 
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
lnScore_BRCA<-CITMIC(BRCA_exp,cl.core=8)
survival<-read.delim("BRCA_survival.txt",header=T,sep = "\t")
survival<-survival[,c("sample","OS","OS.time")]
survival<-survival[which(survival$OS.time!=""),]
survival$sample<-paste0(survival$sample,"A")
survival$sample<-gsub("-",".",survival$sample)

stage<-read.delim("TCGA.BRCA.sampleMap_BRCA_clinicalMatrix",header=T)
stage<-stage[,c("sampleID","pathologic_stage")]
stage<-stage[which(stage$pathologic_stage!=""),]
stage$sampleID<-paste0(stage$sampleID,"A")
stage$sampleID<-gsub("-",".",stage$sampleID)

var_stage_early<-c('Stage IA',"Stage I","Stage IIA","Stage IIC","Stage IIB","I/II NOS","Stage IB","Stage 0","Stage II")
stage_I_II<-stage[which(stage$pathologic_stage%in%var_stage_early),1]

survival<-merge(survival,stage,by.x='sample',by.y='sampleID')
survival<-survival[,c("sample","OS","OS.time")]
survival<-survival[-which(survival$OS.time==0),]


cell_interact_survival<-cbind(survival[which(survival[,1]%in%intersect(survival[,1],colnames(lnScore_BRCA))),],t(lnScore_BRCA)[intersect(survival[,1],colnames(lnScore_BRCA)),])
cell_interact_survival_unI_II<-cell_interact_survival[-which(cell_interact_survival[,1]%in%stage_I_II),]

single_cox_cell_interact_late<-factor_sing_multi(cell_interact_survival_unI_II)
 
single_cox_cell_interact_late<-single_cox_cell_interact_late[,c(1,5)]
BRCA.BQ<-single_cox_cell_interact_late
library(maftools)
load("TCGA-BRCA_SNP.Rdata")
maf.BRCA<-data
maf<-read.maf(maf.BRCA)
data<-maf@data
PIK3CA_mut<-data.frame(data[which(data[,1]=="TP53"),])
PIK3CA_mut_sam<-substr(PIK3CA_mut[,16],1,16)
PIK3CA_mut_sam<-gsub("-",".",PIK3CA_mut_sam)

BRCA.BQ[which(rownames(BRCA.BQ)%in%PIK3CA_mut_sam),6]<-"mutant"
BRCA.BQ[-which(rownames(BRCA.BQ)%in%PIK3CA_mut_sam),6]<-"wt"




library(estimate)

write.table(exp_BRCA_fpkm_tumor_aggregate_log2_20000,file="BRCA_exp.txt",quote=F,sep="\t")
filterCommonGenes(input.f="D:/Users/89800/Desktop/CellRankScore/BRCA_exp.txt", output.f="BRCA_exp.gct", id="GeneSymbol")
estimateScore("BRCA_exp.gct", "BRCA_exp_estimate_score.gct", platform="illumina")
BRCA_score<-read.delim("D:/Users/89800/Desktop/CellRankScore/BRCA_exp_estimate_score.gct",fill=TRUE)

#stem score 
library(xlsx)
TCGA_score_Stemness<-read.xlsx("TCGA_score.xlsx",sheetIndex  =1)
TCGA_score_Stemness<-TCGA_score_Stemness[which(TCGA_score_Stemness[,2]=="BRCA"),c(1,4)]
TCGA_score_Stemness[,1]<-substr(TCGA_score_Stemness[,1],1,16)
TCGA_score_Stemness[,1]<-gsub("-",".",TCGA_score_Stemness[,1])
rownames(TCGA_score_Stemness)<-TCGA_score_Stemness[,1]
rownames(TCGA_score_Stemness)<-substr(rownames(TCGA_score_Stemness),0,12)

BRCA.BQ<-cbind(BRCA.BQ,stemness_Score=TCGA_score_Stemness[rownames(BRCA.BQ),2])
sample<-intersect(TCGA_score_Stemness[,1],single_cox_cell_interact_late[,1])

TCGA_score_Stemness1<-TCGA_score_Stemness[which(TCGA_score_Stemness[,1]%in%sample),]
TCGA_score_Stemness1<-TCGA_score_Stemness1[order(TCGA_score_Stemness1[,2],decreasing = T),]

min_value <- -2
max_value <- 2
rownames(TCGA_score_Stemness2)<-TCGA_score_Stemness1[,1]
colnames(TCGA_score_Stemness2)<-"mRNAsi"
TCGA_score_Stemness2<-data.frame(TCGA_score_Stemness1[,-1])
breaks <- seq(min_value, max_value, length.out = 101)
pheatmap(stem_cell_CTscore[,rownames(TCGA_score_Stemness2)],scale="row",breaks=breaks,show_colnames = F,show_rownames = T,cluster_cols = T,annotation_col=TCGA_score_Stemness2)


BRCA.BQ<-cbind(BRCA.BQ,TCGA_score_Stemness[rownames(BRCA.BQ),2])
colnames(BRCA.BQ)[14]<-"Stemness Score"


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
 #Subtype= c("Basal"="#4ea397","Her2"="#7bd9a5","LumA"="#d0648a","LumB"="#f58db2","Normal"="#ffffff","white"),
 pathologic_stage = c("Stage III"="#008acd","Stage IV"="#d87a80","Stage X"="#ff0000","white"),  
  MSI.status = c("MSS"="#2b821d", "MSI-H"="#005eaa","MSI-L"="#339ca8","white"),
  EBV.positive = c("1"="#516b91", "0"="grey","white"),
  Molecular.Subtype =c("EBV"="#4ea397", "MSI"="#7bd9a5","CIN"="#d0648a","GS"="#f58db2","white"),
 n = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100)
)
STAD.BQ<-STAD.BQ[,-1]
rownames(survival)<-survival[,1]
survival<-survival[,-1]
BRCA.BQ<-cbind(survival[rownames(BRCA.BQ),3],BRCA.BQ)
colnames(BRCA.BQ)[1]<-"pathologic_stage"
p1<-pheatmap(t(PP),cluster_rows = F,cluster_cols = F,annotation_col =STAD.BQ,
         annotation_colors = colors,show_colnames = F,show_rownames = F )
pp<-data.frame(single_cox_cell_interact_late$RiskScore)
rownames(pp)<-single_cox_cell_interact_late[,1]
PP<-data.frame(pp[order(pp[,1],decreasing = T),])
rownames(PP)<-rownames(pp)[order(pp[,1],decreasing = T)]
BRCA.BQ[,1]
library(patchwork)
pheatmap::pheatmap(t(PP),cluster_rows = F,cluster_cols = F,annotation_col =STAD.BQ,
            annotation_colors = colors,show_colnames = F,show_rownames = F )


STAD_score1<-data.frame(STAD_score1)

#柱状图
BRCA.BQ
a<-t(exp_BRCA_fpkm_tumor_aggregate_log2_20000["TP53",rownames(READ.TB)])
comm<-intersect(rownames(READ.TB),rownames(TCGA_score_Stemness))
result<-data.frame(cbind(RiskScore= READ.TB[comm,2],Score=TCGA_score_Stemness[comm,2]))
result<-data.frame(cbind(RiskScore= READ.TB[,2],Score=a))
colnames(kirp_ri)<-c("RiskScore","StromalScore","ImmuneScore","ESTIMATEScore","mRNAsi",
                     "MET","CDKN2A","FGFR3","TERT","SETD2",
                     "CDKN2A_e","FGFR3_e","TERT_e","KIT_e","MET_e","SETD2_e")
ggplot(result, aes(x = RiskScore, y = Score)) +
  geom_point(size = 6, color = "")+
  geom_smooth(method = "lm",aes(x =RiskScore , y = Score),
              se = FALSE,inherit.aes = FALSE,color="#d0648a",linewidth=2)+
  geom_text(aes(x = -19, y = 0.55,label = "r = 0.2974966"),size = 6)+
  geom_text(aes(x = -19, y = 0.57,label = "p = 1.195e-06"),size = 6)


cor.test(result[,1],result[,2])



#箱式图
single_cox_cell_interact_late
library(ggstatsplot)
result[which(result[,2]!="Lifelong Non-smoker"),2]<-"smoker"
result<-data.frame(cbind(RiskScore=READ.TB[rownames(READ.TB),2],Group=as.character(READ.TB[rownames(READ.TB),6])))
result[,1]<-as.numeric(result[,1])
result<-result[-which(result[,2]=="Indeterminate"),]
ggbetweenstats(data=result,Group,RiskScore)
result<-result[which(result[,2]!=''),]
ggplot(data=result,aes(x=Group,y=RiskScore,fill=Group))+
  stat_boxplot(geom = "errorbar",width=0.15)+
  geom_boxplot(fill=c("#ef6d6d","","","","#5470c6"),width=0.3)+
  ggtitle("TP53")+theme(plot.title = element_text(hjust = 0.5))+
