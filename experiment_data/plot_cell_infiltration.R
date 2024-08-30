#SDY144 

#CITMIC
GSE52005<-read.csv("GSE52005_symbol_aggregate_20000.csv",row.names = 1)
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
lm22<-read.delim("D:/Users/89800/Documents/Tencent Files/898003629/FileRecv/LM22.txt",row.names = 1)
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
