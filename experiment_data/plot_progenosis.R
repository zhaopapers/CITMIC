#Scatter plot 
rownames(single_cox_cell_ssgsea)<-single_cox_cell_ssgsea[,1]
single_cox_cell_ssgsea<-single_cox_cell_ssgsea[,-1]
apply(single_cox_cell_ssgsea,2,as.numeric)
bubble2<-apply(single_cox_cell_ssgsea,2,as.numeric)
rownames(bubble2)<-rownames(single_cox_cell_ssgsea)
bubble2<-data.frame(bubble2)
bubble2[which(bubble2$pvalue>0.05),2]<-0
bubble2[which(bubble2$pvalue<0.05&bubble2$pvalue>0.01),2]<-1
bubble2[which(bubble2$pvalue<0.01&bubble2$pvalue>0.001),2]<-2
bubble2[which(bubble2$pvalue<0.001),2]<-3
bubble2[which(bubble2$coef>0),1]<-1
bubble2[which(bubble2$coef<0),1]<-0
bubble2<-bubble2[,-3]
bubble2<-cbind(bubble2,"ssgsea")
bubble2<-cbind(rownames(bubble2),bubble2)

bubble2<-cbind(bubble2,cell_type[rownames(bubble2),])
colnames(bubble2)<-c("cell","HR","P","method","type")
 



rownames(single_cox_cell_interact)<-single_cox_cell_interact[,1]
single_cox_cell_interact<-single_cox_cell_interact[,-1]
apply(single_cox_cell_interact,2,as.numeric)
bubble1<-apply(single_cox_cell_interact,2,as.numeric)
rownames(bubble1)<-rownames(single_cox_cell_interact)
bubble1<-data.frame(bubble1)
bubble1[which(bubble1$pvalue>0.05),2]<-0
bubble1[which(bubble1$pvalue<0.05&bubble1$pvalue>0.01),2]<-1
bubble1[which(bubble1$pvalue<0.01&bubble1$pvalue>0.001),2]<-2
bubble1[which(bubble1$pvalue<0.001),2]<-3
bubble1[which(bubble1$coef>0),1]<-1
bubble1[which(bubble1$coef<0),1]<-0
bubble1<-bubble1[,-3]
colnames(bubble1)[2]<-"Pval"
bubble1<-cbind(bubble1,"intetact_score")
bubble1<-cbind(rownames(bubble1),bubble1)

bubble1<-cbind(bubble1,cell_type[rownames(bubble1),])
colnames(bubble1)<-c("cell","HR","P","method","type")
write.csv(bubble1,"bubble1.csv")
write.csv(bubble2,"bubble2.csv")
bubble<-rbind(bubble1,bubble2)
bubble<-read.csv("bubble.csv")

p_cor <-ggplot(bubble[which(bubble$type%in%"others"),], aes(cell,method)) 

p_others<-p_cor + theme(axis.text.x = element_text(vjust = 0.6, hjust = 0.6, angle = 90),
                         axis.title.y = element_blank(),
                         axis.text.y=element_blank(),
                         legend.position="none")+
  coord_fixed(ratio=1)+
  geom_point(aes(size=P,color=factor(HR),shape=factor(P)), stroke = 1)+
  scale_size(range = c(0,3))+ scale_color_manual(values = c("1" = "#c12e34", "0" = "#0098d9"))+
  scale_shape_manual(values=c("0"=4, "1"=16, "2"=16,"3"=16))
library(patchwork)
p_Lymphoids/(p_Myeloids+p_stem+p_stromal)
