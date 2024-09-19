#correlation pheatmap
lnScore_SKCM<-CITMIC(exp_fpkm_tumor_aggregate_log2_20000,cl.cores = 8)
lnScore_SKCM
cell_inter<-data.frame(matrix(0,nrow=71,ncol=71))
colnames(cell_inter)<-rownames(lnScore_SKCM)
rownames(cell_inter)<-rownames(lnScore_SKCM)

for(i in 1:71){
  for(j in 1:71){
    cor<-cor.test(lnScore_SKCM[i,],lnScore_SKCM[j,],method="pearson")
    if(cor[["p.value"]]<0.05){
      cell_inter[i,j]<-cor[["estimate"]][["cor"]]
      
    }
    
  }
  
}


library(pheatmap)
my_colors <- colorRampPalette(c("#5394cd","white", "#c12e34"))(100)
pheatmap(cell_inter, 
         color = my_colors,
         breaks = seq(-1, 1, length.out = 101))  

#View the top 500 bp most relevant to the subgroup
GSE86363_myeloid<-c("pDCs","Monocytes","Eosinophils","DCs",
"Neutrophils","Mast cell","M1 Macrophages","M2 Macrophages",
"Macrophages","iDCs","mDCs")
cell_myeloid_go<-matrix_cell_go_score[,GSE86363_myeloid]
sum_r<-data.frame(apply(cell_myeloid_go,1,sum))
sum_r<-cbind(rownames(sum_r),sum_r[,1])
sum_r<-data.frame(sum_r)
sum_r[,2]<-as.numeric(sum_r[,2])
sum_r<-sum_r[order(sum_r[,2],decreasing = T),]
edge_list_myeloid<-melt(cell_1_go[sum_r[1:52,1],])
edge_list_myeloid<-edge_list_myeloid[order(edge_list_myeloid[,3],decreasing = T),]
edge_list_myeloid<-edge_list_myeloid[1:500,]
write.csv(edge_list_myeloid,"edge_list_myeloid.csv")


cell_1<-c("Th1 cells","CTLs","Cytotoxic cells","Activated CD8+ T cells",
          "nTreg","DCs","Tregs","Tr1","T cells","iTregs","Immature  B cells",
          "Activated CD4+ T cells","Th17 cells","Tfh","MDSCs")
cell_2<-c("Pericyte","mv Endothelial cells","Smooth muscle cell",
          "Plasma cells","MSCs","Preadipocytes","Muscle cell","Platelets",
          "CD4+ memory T cells")

sum_r<-data.frame(apply(cell_1_go,1,sum))
sum_r<-cbind(rownames(sum_r),sum_r[,1])
sum_r<-data.frame(sum_r)
sum_r[,2]<-as.numeric(sum_r[,2])
sum_r<-sum_r[order(sum_r[,2],decreasing = T),]

edge_list_1<-melt(cell_1_go[sum_r[1:46,1],])#The bp with the highest score in the subgroup is the bp that best represents the correlation within the group
edge_list_1<-edge_list_1[order(edge_list_1[,3],decreasing = T),]
edge_list_1<-edge_list_1[1:500,]
write.csv(edge_list_1,"edge_list_1.csv")

cell_2_go<-matrix_cell_go_score[,cell_2]
sum_r<-data.frame(apply(cell_2_go,1,sum))
sum_r<-cbind(rownames(sum_r),sum_r[,1])
sum_r<-data.frame(sum_r)
sum_r[,2]<-as.numeric(sum_r[,2])
sum_r<-sum_r[order(sum_r[,2],decreasing = T),]

edge_list_2<-melt(cell_2_go[sum_r[1:58,1],])
edge_list_2<-edge_list_2[order(edge_list_2[,3],decreasing = T),]
edge_list_2<-edge_list_2[1:500,]
write.csv(edge_list_2,"edge_list_2.csv")
