#Scatter plot 
rownames(single_cox_cell)<-single_cox_cell[,1]
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


single_cox_cell_interact_late[,3]<-single_cox_cell_interact_late[,3]/30

riskresult<-cell_interact_survival1
res.cut <- surv_cutpoint(data.frame(riskresult), time = "OS.time", event = "OS",
                         variables = c("RiskScore"))
riskresult$OS.time<-riskresult$OS.time/30
riskresult$group<-ifelse(riskresult$RiskScore>=res.cut[["cutpoint"]][["cutpoint"]],"high","low")
kmfit<-survfit(Surv(OS.time,OS)~group,data=riskresult)

ggsurvplot(kmfit, data=riskresult,
           pval =TRUE, 
           xlab="Months",
           title="GSE19234",
           palette=c("#ef6d6d","#5470c6")
           
)
PlotScatter(riskresult,cutoff.point=res.cut[["cutpoint"]][["cutpoint"]])
colnames(single_cox_cell_interact_late)[2:3]<-c("status","time")
PlotScatter<-function(single_cox_cell_interact_late,
                      status.0='Alive',
                      status.1='Dead',
                      TitleYlab_A='Risk Score',
                      TitleYlab_B='Survival Time',
                      TitleXlab='Rank',
                      TitleLegend_A='Risk Group',
                      TitleLegend_B='Status',
                      color.A=c(low='#5470c6',high='#ef6d6d'),
                      color.B=c(status.0='#5470c6',status.1='#ef6d6d'),
                      cutoff.point=NULL
){
  data<-data.frame(riskresult)
  data=data[order(data[,"RiskScore"],decreasing = F),]
  cutoff.point.x<-length(which(data$group=="low"))
  if(is.null(cutoff.point)){
    cutoff.point.y<-median(data[,"RiskScore"])
  }else{
    cutoff.point.y<-cutoff.point
  }
  
  
  data$RiskScore=round(data$RiskScore,1)
  
  `Group` = data$group
  #figure A risk plot
  #rearange colorA
  color.A=c(color.A['low'],color.A['high'])
  names(color.A)=c("low","high")
  fA = ggplot(data = data,
              aes_string(
                x = 1:nrow(data),
                y = data$RiskScore,
                color=factor(`Group`)
              )
  ) +
    geom_point(size = 2) +
    scale_color_manual(name=TitleLegend_A,values = color.A) +
    geom_vline(
      xintercept = cutoff.point.x,
      linetype = 'dotted',
      size = 1
    ) +
    #bg
    theme(
      panel.grid = element_blank(),
      panel.background = element_blank())+
    #x-axis
    theme(
      axis.ticks.x = element_blank(),
      axis.line.x = element_blank(),
      axis.text.x = element_blank(),
      axis.title.x = element_blank()
    ) +
    #y-axis
    theme(
      axis.title.y = element_text(
        size = 14,vjust = 1,angle = 90,family="sans"),
      axis.text.y = element_text(size=11,family = "sans"),
      axis.line.y = element_line(size=0.5,colour = "black"),
      axis.ticks.y = element_line(size = 0.5,colour = "black"))+
    #legend
    theme(legend.title = element_text(size = 13,family = "sans"),
          legend.text = element_text(size=12,family = "sans"))+
    coord_trans()+
    ylab(TitleYlab_A)+
    scale_x_continuous(expand = c(0,3))
  
  cutoff.label=paste0('cutoff: ',round(cutoff.point.y,2))
  
  fA=fA+ annotate("text",
                  x=cutoff.point.x,
                  y=cutoff.point.y,
                  label=cutoff.label,
                  family="sans",
                  size=5,
                  fontface="plain",
                  colour="black")
  
  
  #fB
  color.B=c(color.B['status.0'],color.B['status.1'])
  names(color.B)=c(status.0,status.1)
  fB=ggplot(data = data,
            aes_string(
              x = 1:nrow(data),
              y = data[, "OS.time"],
              color=factor(ifelse(data[,"OS"]==0,status.0,status.1)))
  ) +
    geom_point(size=2)+
    scale_color_manual(name=TitleLegend_B,values = color.B) +
    geom_vline(
      xintercept = cutoff.point.x,
      linetype = 'dotted',
      size = 1
    )  +
    theme(
      panel.grid = element_blank(),
      panel.background = element_blank())+
    #x a-xis
    theme(
      axis.line.x = element_line(size=0.5,colour = "black"),
      axis.ticks.x = element_line(size=0.5,colour = "black"),
      axis.text.x = element_text(size=11,family = "sans"),
      axis.title.x = element_text(size = 14,family="sans")
    ) +
    #y-axis
    theme(
      axis.title.y = element_text(
        size = 14,vjust = 2,angle = 90,family="sans"),
      axis.text.y = element_text(size=11,family = "sans"),
      axis.ticks.y = element_line(size = 0.5),
      axis.line.y = element_line(size=0.5,colour = "black")
    )+
    theme(legend.title = element_text(size = 13,family = "sans"),
          legend.text = element_text(size=12,family = "sans"))+
    ylab(TitleYlab_B)+xlab(TitleXlab)+
    coord_trans()+
    scale_x_continuous(expand = c(0,3))
  
  
  egg::ggarrange(
    fA,
    fB,
    ncol = 1,
    labels = c('A', 'B'),
    label.args = list(gp = grid::gpar(font = 2, cex =1.5,
                                      family="sans"))
  )
  
}

#auc

plotAUCcurve(tROC,conf.int=F,col="blue")
plot(tROC,time=1,col='orange')
plot(tROC,time=3,col='red',add=T)
plot(tROC,time=5,col='green',add=T)
plot(tROC,time=7,col='yellow',add=T)
plot(tROC,time=9,col='blue',add=T)
legend(0.6,0.3,c(paste("AUC of 1 Year =",round(tROC$AUC,3)[1]),
                 paste("AUC of 3 Year =",round(tROC$AUC,3)[2]),#
                 paste("AUC of 5 Year =",round(tROC$AUC,3)[3]),
                 paste("AUC of 7 Year =",round(tROC$AUC,3)[4]),
                 paste("AUC of 9 Year =",round(tROC$AUC,3)[5])),
                 
       x.intersp=0.6, y.intersp=1.0,
       lty= 1 ,lwd= 2,col=c('orange',"red","green","yellow","blue"),
       bty = "n",# bty������???
       seg.len=1,cex=0.9)# 

