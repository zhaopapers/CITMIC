library(data.table)
library(tidyr)
library(NbClust)
library(GSVA)
library(kernlab)
library(survival)
library(survminer)
library(clusterProfiler)
library(igraph)
library(pheatmap)
library(timeROC)
library(ggplot2)
library(glmnet)
load("go_cell_inter_10_350.Rdata")


  go_cell_score_row<-function(a,table,del,c){
      score_row<-rep(0,c)

      for(j in 1:length(a)){
        if(a[j]!=""){
          gene<-unlist(strsplit(a[j], split = ","))
          location<-fastmatch::fmatch(gene, table)

          dell<- na.omit(del[location])
          de_score1<-median(as.numeric(dell))
          if (!is.na(de_score1)) {
            score_row[j]<-de_score1
          }

        }

      }
      return(score_row)
    }
    median_inter<-function(matrix_cell_go_inter,matrix_cell_go_jaccard,GEP){

      GEPscore<-cbind(rownames(GEP),GEP[,1])
      table <- GEPscore[,1]
      del <- GEPscore[, 2]
      median_score<-matrix(0,nrow=length(rownames(matrix_cell_go_inter)),ncol=length(colnames(matrix_cell_go_inter)))
      for(k in 1:length(rownames(matrix_cell_go_inter))){

        Genes_vector<-matrix_cell_go_inter[k,]
        row<-go_cell_score_row(Genes_vector,table,del,length(colnames(matrix_cell_go_inter)))
        median_score[k,]<-row

      }
      matrix_median_genes<-median_score*matrix_cell_go_jaccard
      colnames(matrix_median_genes)<-colnames(matrix_cell_go_inter)
      rownames(matrix_median_genes)<-rownames(matrix_cell_go_inter)
      matrix_cell_score<-t(matrix_median_genes)%*%matrix_median_genes
      matrix_cell_score[is.na(matrix_cell_score)]<-0
      diag(matrix_cell_score)<-0
      return(matrix_cell_score)
    }




  random_crosstalk<-function(result_cell,damping=damping){



    adj.final<-as.matrix(result_cell)
    graph = graph_from_adjacency_matrix(adj.final,mode=c("undirected"),weighted=weighted,add.rownames=T)
    temp = page_rank(graph, vids=V(graph), directed=FALSE, damping=damping, weights=NULL)
    rank = temp$vector
    rank1 = as.matrix(rank)


    return(rank1)
  }

timeRoc<-function(risk,year){
  risk$OS.time<-risk$OS.time/365
  tROC<-timeROC(T=risk$OS.time,marker = risk$RiskScore,
                delta=risk$OS,
                cause=1,weighting="marginal",
                times=c(1:year),
                iid=TRUE)
  print(tROC$AUC)
  
}


factor_sing_multi<-function(cell_survival){
  genes<-colnames(cell_survival)[4:length(colnames(cell_survival))]
  outTab= data.frame()
  for(i in genes){
    expr = cell_survival[,i]
    cox = coxph(Surv(OS.time,OS) ~ expr,cell_survival)
    coxsummary = summary(cox)
    if(coxsummary$coefficients[,"Pr(>|z|)"]<=0.05){
      outTab=rbind(
        outTab,cbind(
          cell=i,
          coef=coxsummary$coefficients[,"coef"],
          HR=coxsummary$coefficients[,"exp(coef)"],
          pvalue=coxsummary$coefficients[,"Pr(>|z|)"]
        )
      )
    }
    #print(i)
  }
  cell_survival<-cell_survival[,c("sample","OS","OS.time",outTab[,1])]
  rownames(cell_survival)<-cell_survival[,1]
  cell_survival<-cell_survival[,-1]
  cox = coxph(Surv(OS.time,OS) ~ .,cell_survival)
  coxsummary = summary(cox)
  print(coxsummary)
  cel<-rownames(coxsummary$coefficients)[which(coxsummary$coefficients[,"Pr(>|z|)"]<=0.05)]
  outTab<-cbind(
    HR=coxsummary$coefficients[cel,"exp(coef)"],
    HR.95L=coxsummary$conf.int[cel,"lower .95"],
    HR.95H=coxsummary$conf.int[cel,"upper .95"],
    beta=coxsummary$coefficients[cel,"coef"],
    pvalue=coxsummary$coefficients[cel,"Pr(>|z|)"]
  )
  colnames(outTab)<-c("HR","HR.95L","HR.95H","beta","P-value")
  
  survival_factor_adjusted<-coxsummary$coefficients[which(coxsummary$coefficients[,5]<=0.05),]

  if(length(survival_factor_adjusted)==5){
    singleGene<-rownames(coxsummary$coefficients)[which(coxsummary$coefficients[,5]<0.05)]
    singleGene<-gsub("`","",singleGene)
    G.exp<-as.matrix(cell_survival[,singleGene])
    RiskScore<-G.exp*survival_factor_adjusted[1]
  }else{
    rownames(survival_factor_adjusted)<-gsub("`","",rownames(survival_factor_adjusted))
    G.exp<-as.matrix(cell_survival[,rownames(survival_factor_adjusted)])

    RiskScore<-G.exp%*%survival_factor_adjusted[,1]
  }


  riskresult = cbind(patient = rownames(cell_survival),cell_survival[,1:2],G.exp,RiskScore)
  res.cut <- surv_cutpoint(data.frame(riskresult), time = "OS.time", event = "OS",
                           variables = c("RiskScore"))

  riskresult$group<-ifelse(riskresult$RiskScore>=res.cut[["cutpoint"]][["cutpoint"]],"high","low")
  print(res.cut[["cutpoint"]][["cutpoint"]])
  return(riskresult)
}
#All TCGA processes are identical, only TCGA-SKCM is listed here for reference.
library(CITMIC)
library(parallel)
load(exp_fpkm_tumor_aggregate_log2_20000,file="exp_fpkm_tumor_aggregate_log2_20000.Rdata")
lnScore_SKCM<-CITMIC(exp_fpkm_tumor_aggregate_log2_20000,cl.cores = 8)

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


cell_interact_survival<-cbind(survival[which(survival[,1]%in%intersect(survival[,1],colnames(lnScore_SKCM))),],t(lnScore_SKCM)[intersect(survival[,1],colnames(lnScore_SKCM)),])
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
clusterExport(cl,varlist = list("GSE22155_GPL6947_aggregate_20000","go_cell_score_row","median_inter","matrix_cell_go_inter","matrix_cell_go_score"))


result_cell_GSE22155_GPL6947_10_350_20000<-parLapply(cl=cl,X=colnames(GSE22155_GPL6947_aggregate_20000),function(x){
                           Zvalue<-data.frame(cbind(rownames(GSE22155_GPL6947_aggregate_20000),GSE22155_GPL6947_aggregate_20000[,x]))
                           Zvalue <- data.frame(Zvalue[!is.infinite(Zvalue[,2]),])
                           Zvalue <- data.frame(Zvalue[!is.na(Zvalue[,2]),])
                           Zvalue <- data.frame(Zvalue[!is.nan(Zvalue[,2]),])
                           Zvalue[,2]<-as.numeric(Zvalue[,2])
                           Z<-data.frame(Zvalue[,2])
                           rownames(Z)<-Zvalue[,1]
                           score<-median_inter(matrix_cell_go_inter,matrix_cell_go_jaccard,Z)
                           return(score)
                         })
names(result_cell_GSE22155_GPL6947_10_350_20000)<-colnames(GSE22155_GPL6947_aggregate_20000)
stopCluster(cl)



Score_rankwalk_GSE22155_GPL6947_20000<-data.frame(row.names=rownames(result_cell_GSE22155_GPL6947_10_350_20000[[1]]))
for(i in names(result_cell_GSE22155_GPL6947_10_350_20000)){
  score<-random_crosstalk(result_cell_GSE22155_GPL6947_10_350_20000[[i]])
  Score_rankwalk_GSE22155_GPL6947_20000<-cbind(
    Score_rankwalk_GSE22155_GPL6947_20000,
    score[rownames(Score_rankwalk_GSE22155_GPL6947_20000),])
  print(i)
}
colnames(Score_rankwalk_GSE22155_GPL6947_20000)<-names(result_cell_GSE22155_GPL6947_10_350_20000)


library(parallel)
cl <- makeCluster(8)
clusterExport(cl,varlist = list("GSE22155_GPL6102_aggregate_20000","go_miRNA_score_row","median_inter","matrix_m_m_score","matrix_cell_go_inter","matrix_cell_go_score"))

result_cell_GSE22155_GPL6102_10_350_20000<-parLapply(cl=cl,X=colnames(GSE22155_GPL6102_aggregate_20000),function(x){
                           Zvalue<-data.frame(cbind(rownames(GSE22155_GPL6102_aggregate_20000),GSE22155_GPL6102_aggregate_20000[,x]))
                           Zvalue <- data.frame(Zvalue[!is.infinite(Zvalue[,2]),])
                           Zvalue <- data.frame(Zvalue[!is.na(Zvalue[,2]),])
                           Zvalue <- data.frame(Zvalue[!is.nan(Zvalue[,2]),])
                           Zvalue[,2]<-as.numeric(Zvalue[,2])
                           Z<-data.frame(Zvalue[,2])
                           rownames(Z)<-Zvalue[,1]
                           score<-median_inter(matrix_cell_go_inter,matrix_cell_go_jaccard,Z)
                           return(score)
                         })
names(result_cell_GSE22155_GPL6102_10_350_20000)<-colnames(GSE22155_GPL6102_aggregate_20000)
stopCluster(cl)
Score_rankwalk_GSE22155_GPL6102_20000<-data.frame(row.names=rownames(result_cell_GSE22155_GPL6102_10_350_20000[[1]]))
for(i in names(result_cell_GSE22155_GPL6102_10_350_20000)){
  score<-random_crosstalk(result_cell_GSE22155_GPL6102_10_350_20000[[i]])
  
  Score_rankwalk_GSE22155_GPL6102_20000<-cbind(
    Score_rankwalk_GSE22155_GPL6102_20000,
    score[rownames(Score_rankwalk_GSE22155_GPL6102_20000),])

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

