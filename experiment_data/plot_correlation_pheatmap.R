#correlation pheatmap
lnScore_SKCM<-CITMIC(exp_fpkm_tumor_aggregate_log2_20000,cl.cores = 8)
lnScore_SKCM
cell_inter<-data.frame(matrix(0,nrow=86,ncol=86))
colnames(cell_inter)<-rownames(lnScore_SKCM)
rownames(cell_inter)<-rownames(lnScore_SKCM)

for(i in 1:86){
  for(j in 1:86){
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
