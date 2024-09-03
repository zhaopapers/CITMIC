#Dumbbell diagram
pancancer_AUC<-read.csv("D:/Users/89800/Desktop/数据整理/D.csv")
pancancer_AUC<-cbind(pancancer_AUC[which(pancancer_AUC[,4]=="all"),1:3],pancancer_AUC[which(pancancer_AUC[,4]!="all"),2])
colnames(pancancer_AUC)[4]<-"iii_auc"
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
#加TP53突变 PD-L1突变
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
 
load("BRCA_exp.Rdata")
lnScore_BRCA<-CITMIC(BRCA_exp,cl.core=8)
single_cox
single_cox_cell_interact_late<-single_cox_cell_interact_late[,c(1,5)]
LIHC.BQ<-single_cox_cell_interact_late
library(maftools)
load("D:/Users/89800/Desktop/CellRankScore/TCGA-BRCA_SNP.Rdata")
maf.KIRP<-data
maf<-read.maf(maf.KIRP)
data<-maf@data
PIK3CA_mut<-data.frame(data[which(data[,1]=="PIK3CA"),])

PIK3CA_mut_sam<-substr(PIK3CA_mut[,16],1,16)
PIK3CA_mut_sam<-gsub("-",".",PIK3CA_mut_sam)

READ.TB[which(rownames(READ.TB)%in%PIK3CA_mut_sam),6]<-"mutant"
READ.TB[-which(rownames(READ.TB)%in%PIK3CA_mut_sam),6]<-"wt"

LUAD.TB<-TB
READ.TB<-single_cox_cell_interact_late[,c(1,5)]
TB

CTNNB1_mut<-data.frame(data[which(data[,2]=="CTNNB1"),])

CTNNB1_mut_sam<-substr(CTNNB1_mut[,17],1,16)
CTNNB1_mut_sam<-gsub("-",".",CTNNB1_mut_sam)

LIHC.BQ[which(rownames(LIHC.BQ)%in%CTNNB1_mut_sam),4]<-"mutant"
LIHC.BQ[-which(rownames(LIHC.BQ)%in%CTNNB1_mut_sam),4]<-"wt"


AXIN1_mut<-data.frame(data[which(data[,2]=="AXIN1"),])

AXIN1_mut_sam<-substr(AXIN1_mut[,17],1,16)
AXIN1_mut_sam<-gsub("-",".",AXIN1_mut_sam)

LIHC.BQ[which(rownames(LIHC.BQ)%in%AXIN1_mut_sam),5]<-"mutant"
LIHC.BQ[-which(rownames(LIHC.BQ)%in%AXIN1_mut_sam),5]<-"wt"

PIK3CA_mut<-data.frame(data[which(data[,2]=="PIK3CA"),])

PIK3CA_mut_sam<-substr(PIK3CA_mut[,17],1,16)
PIK3CA_mut_sam<-gsub("-",".",PIK3CA_mut_sam)

STAD.BQ2[which(rownames(STAD.BQ2)%in%PIK3CA_mut_sam),10]<-"mutant"
STAD.BQ2[-which(rownames(STAD.BQ2)%in%PIK3CA_mut_sam),10]<-"wt"


TP53_mut<-data.frame(data[which(data[,2]=="TP53"),])

TP53_mut_sam<-substr(TP53_mut[,17],1,16)
TP53_mut_sam<-gsub("-",".",TP53_mut_sam)

STAD.BQ2[which(rownames(STAD.BQ2)%in%TP53_mut_sam),11]<-"mutant"
STAD.BQ2[-which(rownames(STAD.BQ2)%in%TP53_mut_sam),11]<-"wt"



colnames(READ.TB)[3]<-"KRAS"
colnames(READ.TB)[4]<-"TP53"
colnames(READ.TB)[5]<-"BRAF"
colnames(READ.TB)[6]<-"PIK3CA"
colnames(READ.TB)[7]<-"EGFR"
colnames(lihc_ri)[16]<-"CDKN2A"
plotmafSummary(maf=maf,rmOutlier = TRUE,addStat = 'median',dashboard = TRUE)

load("pancancer/exp/READ.exp.Rdata")
library(estimate)




single_cox_cell_interact_late[,1]
write.table(exp_READ_fpkm_tumor_aggregate_log2_20000[,single_cox_cell_interact_late[,1]],file="READ_exp.txt",quote=F,sep="\t")
filterCommonGenes(input.f="D:/Users/89800/Desktop/CellRankScore/READ_exp.txt", output.f="READ_exp.gct", id="GeneSymbol")
estimateScore("READ_exp.gct", "READ_exp_estimate_score.gct", platform="illumina")
STAD_score<-read.delim("D:/Users/89800/Desktop/CellRankScore/READ_exp_estimate_score.gct",fill=TRUE)
STAD_score<-t(STAD_score[2:5,])
colnames(STAD_score)<-STAD_score[1,]
STAD_score<-STAD_score[-1,]
rownames(STAD_score)<-STAD_score[,1]
STAD_score<-STAD_score[,-1]
STAD_score<-STAD_score[-1,]
STAD_score<-data.frame(STAD_score)
STAD_score1<-apply(STAD_score,2,as.numeric)
rownames(STAD_score1)<-rownames(STAD_score)
STAD.BQ1<-cbind(STAD.BQ,STAD_score1[rownames(STAD.BQ),])
STAD_score1<-data.frame(STAD_score1)
rownames(STAD_score1)<-substr(rownames(STAD_score1),0,12)
rownames(phe)<-phe[,1]
phe<-phe[,-1]
stad_ri<-single_cox_cell_interact_late[,c(1,9)]

comm<-intersect(rownames(phe),rownames(STAD_score1))
STAD_score1<-cbind(STAD_score1,phe[rownames(STAD_score1),])
#干性得分 
library(xlsx)
TCGA_score_Stemness<-read.xlsx("D://a研究生报告数据//a研究生报告数据//赵希龙-20231014干性得分//stemness_score//TCGA_score.xlsx",sheetIndex  =1)
TCGA_score_Stemness<-TCGA_score_Stemness[which(TCGA_score_Stemness[,2]=="READ"),c(1,4)]
TCGA_score_Stemness[,1]<-substr(TCGA_score_Stemness[,1],1,16)
TCGA_score_Stemness[,1]<-gsub("-",".",TCGA_score_Stemness[,1])
rownames(TCGA_score_Stemness)<-TCGA_score_Stemness[,1]
rownames(TCGA_score_Stemness)<-substr(rownames(TCGA_score_Stemness),0,12)

STAD.BQ2<-cbind(STAD.BQ1,stemness_Score=TCGA_score_Stemness[rownames(STAD.BQ1),2])
stem_cell<-c("Chondrocytes","CLP","CMP","Erythrocytes","GMP","HSCs","MDSCs","MEP","MPP","MSCs","Platelets")
sample<-intersect(TCGA_score_Stemness[,1],single_cox_cell_interact_late[,1])
stem_cell_CTscore<-network_cell_score[stem_cell,sample]
TCGA_score_Stemness1<-TCGA_score_Stemness[which(TCGA_score_Stemness[,1]%in%sample),]
TCGA_score_Stemness1<-TCGA_score_Stemness1[order(TCGA_score_Stemness1[,2],decreasing = T),]

min_value <- -2
max_value <- 2
rownames(TCGA_score_Stemness2)<-TCGA_score_Stemness1[,1]
colnames(TCGA_score_Stemness2)<-"mRNAsi"
TCGA_score_Stemness2<-data.frame(TCGA_score_Stemness1[,-1])
breaks <- seq(min_value, max_value, length.out = 101)
pheatmap(stem_cell_CTscore[,rownames(TCGA_score_Stemness2)],scale="row",breaks=breaks,show_colnames = F,show_rownames = T,cluster_cols = T,annotation_col=TCGA_score_Stemness2)


BRCA.BQ2<-cbind(BRCA.BQ1,TCGA_score_Stemness1[rownames(BRCA.BQ1),2])
colnames(BRCA.BQ2)[14]<-"Stemness Score"
subtype1<-rownames(information)[which(information[,1]=="1")]


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
STAD.BQ3<-STAD.BQ3[,-1]
rownames(survival)<-survival[,1]
survival<-survival[,-1]
BRCA.BQ<-cbind(survival[rownames(BRCA.BQ),3],BRCA.BQ)
colnames(BRCA.BQ)[1]<-"pathologic_stage"
p1<-pheatmap(t(PP),cluster_rows = F,cluster_cols = F,annotation_col =STAD.BQ3,
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
a<-t(exp_READ_fpkm_tumor_aggregate_log2_20000["PDCD1",rownames(READ.TB)])
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
