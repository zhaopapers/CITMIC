#correlation pheatmap

cell_inter<-data.frame(matrix(0,nrow=86,ncol=86))
colnames(cell_inter)<-rownames(Inscore_SKCM)
rownames(cell_inter)<-rownames(Inscore_SKCM)


result<-c()
for(i in 1:86){
  for(j in 1:86){
    cor<-cor.test(Inscore_SKCM[i,],Inscore_SKCM[j,],method="pearson")
    if(cor[["p.value"]]<0.05){
      cell_inter[i,j]<-cor[["estimate"]][["cor"]]
      
    }
    
  }
  
}


library(pheatmap)


my_colors <- colorRampPalette(c("#5394cd","white", "#c12e34"))(100)


pheatmap(cell_inter, 
         color = my_colors,
         breaks = seq(-1, 1, length.out = 101))  # breaks 设置颜色映射的范围，使 0 为白色
pheatmap(cell_inter)
