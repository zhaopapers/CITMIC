#ID transfer function
ID_Tran_Gene_GEO<-function(a,b){
  probe_gene<-cbind(b[,1],b[,2])
  probe_location<-c()
  Probes<-c()
  for(i in 1:dim(probe_gene)[1]){
    if(grepl("///",probe_gene[i,2])){
      symbols<-unlist(strsplit(probe_gene[i,2],"///"))
      probes<-cbind(rep(probe_gene[i,1],length(symbols)),symbols)
      Probes<-rbind(Probes,probes)
      probe_location<-c(probe_location,i)
    }
  }
  Probe_gene<-probe_gene[-probe_location,]
  Probe_gene<-rbind(Probe_gene,Probes)
  Probe_gene<-Probe_gene[-which(Probe_gene[,2]==""),]
  colnames(Probe_gene)<-c("ID","Symbol")
  Probe_gene<-data.frame(Probe_gene)
  colnames(a)[1]<-"ID"
  data_geo_exp<-merge(Probe_gene,a,by.x="ID",by.y="ID")
  
  sy<-unique(data_geo_exp[,2])
  DeleteGene<-c()
  data_MiRNA_gene_exp1<-c()
  for(i in 1:length(sy)){
    matrix<-data_geo_exp[which(data_geo_exp[,2]==sy[i]),]
    if(dim(matrix)[1]>=2){
      
      deleteGene<-which(data_geo_exp[,2]==sy[i])
      genes<-matrix[1,2]
      y= matrix[,3:dim(data_geo_exp)[2]]
      y=apply(y,2,as.numeric)
      sum1<-colMeans(y)
      sum1<-c(genes,sum1)
      
      DeleteGene<-c(DeleteGene,deleteGene)
      data_MiRNA_gene_exp1<-rbind(data_MiRNA_gene_exp1,sum1)
    }
    print(i)
  }
  data_geo_exp<-data_geo_exp[-DeleteGene,]
  colnames(data_geo_exp)[2]<-c("gene")
  colnames(data_MiRNA_gene_exp1)[1]<-c("gene")
  data_geo_exp<-rbind(data_geo_exp[,2:dim(data_geo_exp)[2]],data_MiRNA_gene_exp1)
  rownames(data_geo_exp)<-data_geo_exp[,1]
  data_geo_exp<-data_geo_exp[,-1]
  
  return(data_geo_exp)
}



# cell infiltration
#SDY144 GSE22155
GSE52005<-read.delim("GSE52005_symbol_aggregate_20000.csv")

#SDY67 
SDY67<-read.delim("SDY67_tmp_20000.csv") # Select the SDY67 attempted in MCPcounter

#SDY420 SDY311 from xCell preprocess
load("sdy311.rds")
load("sdy420.rds")

#GSE86363 
library(data.table)
GSE86363<-fread("GSE86363_xp133A.txt",sep = "\t")
GSE86363<-data.frame(GSE86363)
GPL96<-read.csv("GPL96.csv")
GPL96<-GPL96[,c("ID","Gene.Symbol")]
exp_geo_GSE86363<-merge(GPL96,GSE86363,by.x="ID",by.y="V1")
exp_geo_GSE86363<-exp_geo_GSE86363[,-1]
GSE86363_symbol<-ID_Tran_Gene_GEO(GSE86363,GPL96)
GSE86363_symbol_1<-apply(GSE86363_symbol,2,as.numeric)
rownames(GSE86363_symbol_1)<-rownames(GSE86363_symbol)
GSE86363_symbol_log2<-log2(GSE86363_symbol_1+1)
threshold <- 0.5 * ncol(GSE86363_symbol_log2)
GSE86363_symbol_log2_20000<-GSE86363_symbol_log2[rowSums(GSE86363_symbol_log2[, -1] == 0) < threshold, ]

#TCGA data preprocess 
#All TCGA processes are identical, only TCGA-SKCM is listed here for reference.
library(tibble)
library(clusterProfiler)
exp_fpkm<-read.table("TCGA-SKCM.htseq_fpkm.tsv",header=T) # SKCM
mapp<-read.table("gencode.v22.annotation.gene.probeMap",header = T)#map

exp_fpkm<-merge(mapp,exp_fpkm,by.x="id",by.y="Ensembl_ID")
exp_fpkm_tumor<-exp_fpkm[,c("gene",colnames(exp_fpkm)[which(substr(colnames(exp_fpkm),14,14)%in%c("0"))])]
colnames(exp_fpkm_tumor)[1]<-"gene"
exp_fpkm_tumor[,2:ncol(exp_fpkm_tumor)]<-(2^exp_fpkm_tumor[,2:ncol(exp_fpkm_tumor)])-1
exp_fpkm_tumor_aggregate<-aggregate(exp_fpkm_tumor[,-1],by=list(exp_fpkm_tumor$gene),median)

gene <- bitr(exp_fpkm_tumor_aggregate[,1],fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = "org.Hs.eg.db")
exp_fpkm_tumor_aggregate<-exp_fpkm_tumor_aggregate[which(exp_fpkm_tumor_aggregate[,1]%in%gene[,1]),]
rownames(exp_fpkm_tumor_aggregate)<-NULL

exp_fpkm_tumor_aggregate<-column_to_rownames(exp_fpkm_tumor_aggregate,var=colnames(exp_fpkm_tumor_aggregate)[1])
exp_fpkm_tumor_aggregate_log2<-log2(exp_fpkm_tumor_aggregate+1)

threshold <- 0.5 * ncol(exp_fpkm_tumor_aggregate_log2)
exp_fpkm_tumor_aggregate_log2_20000<-exp_fpkm_tumor_aggregate_log2[rowSums(exp_fpkm_tumor_aggregate_log2[, -1] == 0) < threshold, ]

#vailation 
#GSE22155 GSE19234
GSE22155_GPL6102<-read.delim("vaildation/GSE22155/GSE22155-GPL6102_series_matrix.txt")
GSE22155_GPL6947<-read.delim("vaildation/GSE22155/GSE22155-GPL6947_series_matrix.txt")

GPL6102<-na.omit(read.csv("vaildation/GSE22155/gpl6102.csv"))
GPL6947<-na.omit(read.csv("vaildation/GSE22155/GPL6947.csv"))

GSE22155_survival<-read.csv("vaildation/GSE22155/GSE22155_survival.csv")

GPL6102<-GPL6102[which(GPL6102[,2]!=""),]
GPL6947<-GPL6947[which(GPL6947[,2]!=""),]

GSE22155_GPL6102<-merge(GPL6102,GSE22155_GPL6102,by.x="ID",by.y="ID_REF")
GSE22155_GPL6947<-merge(GPL6947,GSE22155_GPL6947,by.x="Probe_Id",by.y="ID_REF")

GSE22155_GPL6102<-GSE22155_GPL6102[,-1]  
GSE22155_GPL6947<-GSE22155_GPL6947[,-1]

GSE22155_GPL6102_aggregate<-aggregate(GSE22155_GPL6102[,-1],by=list(GSE22155_GPL6102$Symbol),median)
GSE22155_GPL6947_aggregate<-aggregate(GSE22155_GPL6947[,-1],by=list(GSE22155_GPL6947$Symbol),median)

rownames(GSE22155_GPL6102_aggregate)<-NULL
GSE22155_GPL6102_aggregate<-column_to_rownames(GSE22155_GPL6102_aggregate,var=colnames(GSE22155_GPL6102_aggregate)[1])

rownames(GSE22155_GPL6947_aggregate)<-NULL
GSE22155_GPL6947_aggregate<-column_to_rownames(GSE22155_GPL6947_aggregate,var=colnames(GSE22155_GPL6947_aggregate)[1])

threshold <- 0.5 * ncol(GSE22155_GPL6102_aggregate)
GSE22155_GPL6102_aggregate_20000<-GSE22155_GPL6102_aggregate[rowSums(GSE22155_GPL6102_aggregate == 0) < threshold, ]

threshold <- 0.5 * ncol(GSE22155_GPL6947_aggregate)
GSE22155_GPL6947_aggregate_20000<-GSE22155_GPL6947_aggregate[rowSums(GSE22155_GPL6947_aggregate == 0) < threshold, ]

#GSE19234
GPL570<-read.csv("vaildation/GSE19234/GPL570.csv")
GSE19234<-read.delim("vaildation/GSE19234/GSE19234_series_matrix.txt")
GSE19234_GPL570<-merge(GPL570,GSE19234,by.x="ID",by.y="ID_REF")
GSE19234_GPL570<-GSE19234_GPL570[,-1]
GSE19234_GPL570_aggregate<-aggregate(GSE19234_GPL570[,-1],by=list(GSE19234_GPL570$Symbol),median)
rownames(GSE19234_GPL570_aggregate)<-NULL
GSE19234_GPL570_aggregate<-column_to_rownames(GSE19234_GPL570_aggregate,var=colnames(GSE19234_GPL570_aggregate)[1])
threshold <- 0.5 * ncol(GSE19234_GPL570_aggregate)
GSE19234_GPL570_aggregate_20000<-GSE19234_GPL570_aggregate[rowSums(GSE19234_GPL570_aggregate == 0) < threshold, ]

