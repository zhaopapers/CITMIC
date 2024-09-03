load("result_score_TCGA_Cancer.Rdata")

var_names<-ls()
pancancer_sample<-c()
for(i in var_names){
  pancancer_type<-cbind(substr(i,16,19),colnames(get(i)))
  pancancer_sample <- rbind(pancancer_sample,pancancer_type)
  print(i)
}

for(i in var_names){
  pancancer<-cbind(pancancer,get(i)[rownames(pancancer),])
  print(i)
}
cell_num<-apply(pancancer,2,function(x){
  y=log(x,base=10)
  return(y)
})

pan_network_cell_score<-(cell_num-min(cell_num))/(max(cell_num)-min(cell_num))
library(Rtsne)
tsne_out1 <- Rtsne(t(pan_network_cell_score),pca=TRUE,theta=0.0,num_theads=0) # Run TSNE
tsne_result2 = as.data.frame(tsne_out1$Y)
colnames(tsne_result2) = c("tSNE1","tSNE2")

pancancer_sample<-data.frame(pancancer_sample)
pancancer_sample[,1]<-gsub("GBM_","GBM",pancancer_sample[,1])
pancancer_sample[,1]<-gsub("ACC_","GBM",pancancer_sample[,1])
pancancer_sample[,1]<-gsub("UVM_","UVM",pancancer_sample[,1])
pancancer_sample[,1]<-gsub("OV_f","OV",pancancer_sample[,1])
pancancer_sample[,1]<-gsub("LGG_","LGG",pancancer_sample[,1])
pancancer_sample[,1]<-gsub("UCS_","UCS",pancancer_sample[,1])

gs<-cbind(tsne_result2,pancancer_sample)
sample_labels <- gs%>% group_by(X1)%>%slice_sample(n = 1)%>%ungroup()
centroids <- aggregate(cbind(tSNE1, tSNE2) ~ X1, gs, mean)
ggplot(gs, aes(x = tSNE1, y = tSNE2, color = X1)) +
  geom_point() +
  geom_text_repel(data = centroids, aes(label = X1), size = 5, color = "black") +
  theme_minimal()
