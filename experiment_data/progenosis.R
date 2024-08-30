#All TCGA processes are identical, only TCGA-SKCM is listed here for reference.
library(CITMIC)
library(parallel)
load(exp_fpkm_tumor_aggregate_log2_20000,file="exp_fpkm_tumor_aggregate_log2_20000.Rdata")
lnScore<-CITMIC(exp_fpkm_tumor_aggregate_log2_20000,cl.cores = 8)
lnScore1<-CITMIC(exp_SKCM_fpkm_tumor_aggregate_log2_20000,cl.cores = 8)
#prognosis
survival<-read.delim("SKCM_survival.txt",header=T,sep = "\t")
survival<-survival[,c("sample","OS","OS.time")]
survival<-survival[which(survival$OS.time!=""),]
survival$sample<-paste0(survival$sample,"A")
survival$sample<-gsub("-",".",survival$sample)

stage<-read.delim("TCGA.SKCM.sampleMap_SKCM_clinicalMatrix",header=T)
stage<-stage[,c("sampleID","pathologic_stage")]
stage<-stage[which(stage$pathologic_stage!=""),]
stage$sampleID<-paste0(stage$sampleID,"A")
stage$sampleID<-gsub("-",".",stage$sampleID)

var_stage_early<-c('Stage IA',"Stage I","Stage IIA","Stage IIC","Stage IIB","I/II NOS","Stage IB","Stage 0","Stage II")
stage_I_II<-stage[which(stage$pathologic_stage%in%var_stage_early),1]

survival<-merge(survival,stage,by.x='sample',by.y='sampleID')
survival<-survival[,c("sample","OS","OS.time")]
survival<-survival[-which(survival$OS.time==0),]


cell_interact_survival<-cbind(survival[which(survival[,1]%in%intersect(survival[,1],colnames(lnScore))),],t(lnScore)[intersect(survival[,1],colnames(lnScore)),])
cell_interact_survival_unI_II<-cell_interact_survival[-which(cell_interact_survival[,1]%in%stage_I_II),]

single_cox_cell_interact<-factor_sing_multi(cell_interact_survival)
single_cox_cell_interact_late<-factor_sing_multi(cell_interact_survival_unI_II)

timeRoc(single_cox_cell_interact,10)
timeRoc(single_cox_cell_interact_late,10)


#vailation set GSE22155 
load("GSE22155_GPL6947_aggregate_20000.Rdata")
load("GSE22155_GPL6102_aggregate_20000.Rdata")
library(parallel)
cl <- makeCluster(8)
clusterExport(cl,varlist = list("GSE22155_GPL6947_aggregate_20000","go_miRNA_score_row","median_inter","matrix_m_m_score","matrix_cell_go_inter","matrix_cell_go_score"))

result_cell_GSE22155_GPL6947_10_350_20000<-parLapply(cl=cl,
                                                     X=colnames(GSE22155_GPL6947_aggregate_20000),function(x){
                                                       Zvalue<-data.frame(cbind(rownames(GSE22155_GPL6947_aggregate_20000),GSE22155_GPL6947_aggregate_20000[,x]))
                                                       Zvalue <- data.frame(Zvalue[!is.infinite(Zvalue[,2]),])
                                                       Zvalue <- data.frame(Zvalue[!is.na(Zvalue[,2]),])
                                                       Zvalue <- data.frame(Zvalue[!is.nan(Zvalue[,2]),])
                                                       Zvalue[,2]<-as.numeric(Zvalue[,2])
                                                       Z<-data.frame(Zvalue[,2])
                                                       rownames(Z)<-Zvalue[,1]
                                                       median_num_SKCM<-median_inter(matrix_cell_go_inter,matrix_cell_go_score,Z)
                                                       score_SKCM<-matrix_m_m_score(median_num_SKCM)
                                                       diag(score_SKCM)<-0
                                                       
                                                       return(score_SKCM)
                                                       
                                                       
                                                     })
names(result_cell_GSE22155_GPL6947_10_350_20000)<-colnames(GSE22155_GPL6947_aggregate_20000)
stopCluster(cl)



Score_rankwalk_GSE22155_GPL6947_20000<-data.frame(row.names=rownames(result_cell_GSE22155_GPL6947_10_350_20000[[1]]))
for(i in names(result_cell_GSE22155_GPL6947_10_350_20000)){
  score<-final_score(result_cell_GSE22155_GPL6947_10_350_20000[[i]])
  Score_rankwalk_GSE22155_GPL6947_20000<-cbind(
    Score_rankwalk_GSE22155_GPL6947_20000,
    score[rownames(Score_rankwalk_GSE22155_GPL6947_20000),])
  print(i)
}
colnames(Score_rankwalk_GSE22155_GPL6947_20000)<-names(result_cell_GSE22155_GPL6947_10_350_20000)


library(parallel)
cl <- makeCluster(8)
clusterExport(cl,varlist = list("GSE22155_GPL6102_aggregate_20000","go_miRNA_score_row","median_inter","matrix_m_m_score","matrix_cell_go_inter","matrix_cell_go_score"))

result_cell_GSE22155_GPL6102_10_350_20000<-parLapply(cl=cl,
                                                     X=colnames(GSE22155_GPL6102_aggregate_20000),function(x){
                                                       Zvalue<-data.frame(cbind(rownames(GSE22155_GPL6102_aggregate_20000),GSE22155_GPL6102_aggregate_20000[,x]))
                                                       Zvalue <- data.frame(Zvalue[!is.infinite(Zvalue[,2]),])
                                                       Zvalue <- data.frame(Zvalue[!is.na(Zvalue[,2]),])
                                                       Zvalue <- data.frame(Zvalue[!is.nan(Zvalue[,2]),])
                                                       Zvalue[,2]<-as.numeric(Zvalue[,2])
                                                       Z<-data.frame(Zvalue[,2])
                                                       rownames(Z)<-Zvalue[,1]
                                                       median_num_SKCM<-median_inter(matrix_cell_go_inter,matrix_cell_go_score,Z)
                                                       score_SKCM<-matrix_m_m_score(median_num_SKCM)
                                                       diag(score_SKCM)<-0
                                                       return(score_SKCM)
                                                       
                                                       
                                                     })
names(result_cell_GSE22155_GPL6102_10_350_20000)<-colnames(GSE22155_GPL6102_aggregate_20000)
stopCluster(cl)
Score_rankwalk_GSE22155_GPL6102_20000<-data.frame(row.names=rownames(result_cell_GSE22155_GPL6102_10_350_20000[[1]]))
for(i in names(result_cell_GSE22155_GPL6102_10_350_20000)){
  score<-final_score(result_cell_GSE22155_GPL6102_10_350_20000[[i]])
  
  Score_rankwalk_GSE22155_GPL6102_20000<-cbind(
    Score_rankwalk_GSE22155_GPL6102_20000,
    score[rownames(Score_rankwalk_GSE22155_GPL6102_20000),])
  #score<-cbind(rownames(score),score[,1])
  #colnames(score)<-c("cell",names(result_RNA_cell_TCGA_SKCM_fpkm)[i])
  #Score_rankwalk_TCGA_SKCM_fpkm<-merge(Score_rankwalk_TCGA_SKCM_fpkm,score,by="cell")
  print(i)
}
colnames(Score_rankwalk_GSE22155_GPL6102_20000)<-names(result_cell_GSE22155_GPL6102_10_350_20000)

Score_rankwalk_GSE22155_20000<-cbind(Score_rankwalk_GSE22155_GPL6102_20000,Score_rankwalk_GSE22155_GPL6947_20000[rownames(Score_rankwalk_GSE22155_GPL6102_20000),])

cell_num<-apply(Score_rankwalk_GSE22155_20000,2,function(x){
  y=log(x,base=10)
  return(y)
})
network_cell_score<-(cell_num-min(cell_num))/(max(cell_num)-min(cell_num))


cell_interact_survival_GSE22155<-cbind(GSE22155_survival[which(GSE22155_survival[,1]%in%intersect(GSE22155_survival[,1],colnames(network_cell_score))),],t(network_cell_score)[intersect(GSE22155_survival[,1],colnames(network_cell_score)),])
colnames(cell_interact_survival_GSE22155)[1:3]<-c("sample","OS","OS.time")
risk<-as.matrix(cell_interact_survival_GSE22155[,rownames(single_cox_cell_interact_late_beta)])%*%as.numeric(single_cox_cell_interact_late_beta[,"beta"])
cell_interact_survival_GSE22155<-cbind(cell_interact_survival_GSE22155[,1:3],RiskScore=risk)

timeRoc(cell_interact_survival_GSE22155,5)



load("GSE19234_GPL570_aggregate_20000.Rdata")
lnScore_GSE19234<-CITMIC(GSE19234_GPL570_aggregate_20000,cl.cores = 8)
GSE19234_survival<-read.csv("vaildation/GSE19234/GSE19234_survival.csv")
cell_interact_survival_GSE19234<-cbind(GSE19234_survival[which(GSE19234_survival[,1]%in%intersect(GSE19234_survival[,1],colnames(lnScore_GSE19234))),],t(lnScore_GSE19234)[intersect(GSE19234_survival[,1],colnames(lnScore_GSE19234)),])
risk<-as.matrix(cell_interact_survival_GSE19234[,rownames(single_cox_cell_interact_late_beta)])%*%as.numeric(single_cox_cell_interact_late_beta[,"beta"])
cell_interact_survival_GSE19234<-cbind(cell_interact_survival_GSE19234[,1:3],RiskScore=risk)


timeRoc(cell_interact_survival_GSE19234,9)

