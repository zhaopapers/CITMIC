#single cox bubble
factor_s<-function(cell_survival){
  cell<-cell_survival[4:length(colnames(cell_survival))]
  cells<-colnames(cell)
  outTab= data.frame()
  for(i in cells){
    expr = cell_survival[,i]
    cox = coxph(Surv(OS.time,OS) ~ expr,cell_survival)
    coxsummary = summary(cox)
      outTab=rbind(
        outTab,cbind(
          cell=i,
          coef=coxsummary$coefficients[,"coef"],
          HR=coxsummary$coefficients[,"exp(coef)"],
          pvalue=coxsummary$coefficients[,"Pr(>|z|)"]
        )
      )
    print(i)
  }
  G.exp<-as.matrix(cell_survival[,outTab[,1]])
  print(outTab)
  coef<-as.numeric(outTab[,2])
  RiskScore<-G.exp%*%coef
  riskresult = cbind(patient = cell_survival[,1],cell_survival[,2:3],G.exp,RiskScore)
  res.cut <- surv_cutpoint(data.frame(riskresult), time = "OS.time", event = "OS",
                           variables = c("RiskScore"))
  riskresult$group<-ifelse(riskresult$RiskScore>=res.cut[["cutpoint"]][["cutpoint"]],"high","low")
  return(outTab)
  
}


cell_type<-read.csv("cell.type.csv",header=F,row.names = 1)
lnScore_SKCM<-CITMIC(exp_fpkm_tumor_aggregate_log2_20000,cl.cores = 8)

#ssgsea t(scale(t(exp_SKCM_fpkm_tumor_aggregate_log2_20000)))
library(GSVA) 
SKCM_ssGSEA<-gsva(as.matrix(t(scale(t(exp_fpkm_tumor_aggregate_log2_20000)))),TMEcell_list,method="ssgsea",kcdf="Gaussian")

#cibersort  exp_fpkm_tumor_aggregate_20000
SKCM_cibersort<-my_CIBERSORT(exp_fpkm_tumor_aggregate_20000,lm22, perm=100, QN=TRUE, cores = 3)
SKCM_cibersort<-SKCM_cibersort$proportions
SKCM_cibersort<-t(SKCM_cibersort)
SKCM_cibersort<-SKCM_cibersort[1:22,]

#quanTIseq
library(immunedeconv)
SKCM_quantiseq<-deconvolute(as.matrix(exp_fpkm_tumor_aggregate_20000),"quantiseq")
SKCM_quantiseq<-data.frame(SKCM_quantiseq)
rownames(SKCM_quantiseq)<-SKCM_quantiseq[,1]
SKCM_quantiseq<-SKCM_quantiseq[,-1]

#MCPcounter
library(MCPcounter)
SKCM_MCPcounter<-MCPcounter.estimate(exp_fpkm_tumor_aggregate_20000,featuresType=c("affy133P2_probesets","HUGO_symbols","ENTREZ_ID","ENSEMBL_ID")[2],
                                    probesets=read.table("probesets.txt",sep="\t",stringsAsFactors=FALSE,colClasses="character"),
                                    genes=read.table("genes.txt",sep="\t",stringsAsFactors=FALSE,header=TRUE,colClasses="character",check.names=FALSE)
                                    
)

#EPIC exp_fpkm_tumor_aggregate_20000
library(EPIC)
SKCM_EPIC<-deconvolute(exp_fpkm_tumor_aggregate_20000,"epic")
SKCM_EPIC<-data.frame(SKCM_EPIC)
rownames(SKCM_EPIC)<-SKCM_EPIC[,1]
SKCM_EPIC<-SKCM_EPIC[,-1]

#xCell   exp_SKCM_fpkm_tumor_aggregate_log2_20000
library(xCell)
SKCM_xcell <-xCellAnalysis(exp_fpkm_tumor_aggregate_log2_20000)
SKCM_xcell<-SKCM_xcell[1:64,]

cell_interact_survival<-cbind(survival[which(survival[,1]%in%intersect(survival[,1],colnames(lnScore_SKCM))),],t(lnScore_SKCM)[intersect(survival[,1],colnames(lnScore_SKCM)),])
single_cox_cell_CITMIC<-factor_s(cell_interact_survival)
write.csv(single_cox_cell_CITMIC,"single_cox_cell_CITMIC.csv")

cell_interact_survival<-cbind(survival[which(survival[,1]%in%intersect(survival[,1],colnames(SKCM_ssGSEA))),],t(SKCM_ssGSEA)[intersect(survival[,1],colnames(SKCM_ssGSEA)),])
single_cox_cell_ssGSEA<-factor_s(cell_interact_survival)
write.csv(single_cox_cell_ssGSEA,"single_cox_cell_ssGSEA.csv")

cell_interact_survival<-cbind(survival[which(survival[,1]%in%intersect(survival[,1],colnames(SKCM_xcell))),],t(SKCM_xcell)[intersect(survival[,1],colnames(SKCM_xcell)),])
single_cox_cell_xcell<-factor_s(cell_interact_survival)
write.csv(single_cox_cell_xcell,"single_cox_cell_xcell.csv")

cell_interact_survival<-cbind(survival[which(survival[,1]%in%intersect(survival[,1],colnames(tcga_epic))),],t(tcga_epic)[intersect(survival[,1],colnames(tcga_epic)),])
single_cox_cell_EPIC<-factor_s(cell_interact_survival)
write.csv(single_cox_cell_EPIC,"single_cox_cell_EPIC.csv")

cell_interact_survival<-cbind(survival[which(survival[,1]%in%intersect(survival[,1],colnames(SKCM_MCPcounter))),],t(SKCM_MCPcounter)[intersect(survival[,1],colnames(SKCM_MCPcounter)),])
single_cox_cell_MCPcounter<-factor_s(cell_interact_survival)
write.csv(single_cox_cell_MCPcounter,"single_cox_cell_MCPcounter.csv")

cell_interact_survival<-cbind(survival[which(survival[,1]%in%intersect(survival[,1],colnames(SKCM_cibersort))),],t(SKCM_cibersort)[intersect(survival[,1],colnames(SKCM_cibersort)),])
single_cox_cell_cibersort<-factor_s(cell_interact_survival)
write.csv(single_cox_cell_cibersort,"single_cox_cell_cibersort.csv")

cell_interact_survival<-cbind(survival[which(survival[,1]%in%intersect(survival[,1],colnames(SKCM_quantiseq))),],t(SKCM_quantiseq)[intersect(survival[,1],colnames(SKCM_quantiseq)),])
single_cox_cell_quantiseq<-factor_s(cell_interact_survival)
write.csv(single_cox_cell_quantiseq,"single_cox_cell_quantiseq.csv")


rownames(single_cox_cell_ssGSEA)<-single_cox_cell_ssGSEA[,1]
single_cox_cell_ssGSEA<-single_cox_cell_ssGSEA[,-1]
bubble<-apply(single_cox_cell_ssGSEA,2,as.numeric)
rownames(bubble)<-rownames(single_cox_cell_ssGSEA)
bubble<-data.frame(bubble)
bubble[which(bubble$pvalue>0.05),2]<-0
bubble[which(bubble$pvalue<0.05&bubble$pvalue>0.01),2]<-1
bubble[which(bubble$pvalue<0.01&bubble$pvalue>0.001),2]<-2
bubble[which(bubble$pvalue<0.001),2]<-3
bubble[which(bubble$coef>0),1]<-1
bubble[which(bubble$coef<0),1]<-0
bubble<-bubble[,-3]
bubble<-cbind(bubble,"ssGSEA")
bubble<-cbind(rownames(bubble),bubble)

bubble<-cbind(bubble,cell_type[rownames(bubble),])
colnames(bubble)<-c("cell","HR","P","method","type")

write.csv(bubble,"bubble_CITMIC.csv")

#Integration different methods for the same cell type
dotpot<-read.csv("cell_interact_Lymphoids_cell.csv")
dotpot[,3]<-as.factor(dotpot[,3])
dotpot[,2]<-as.factor(dotpot[,2])


p_Lymphoids_cell<-ggplot(Dotpot, aes(cell,method))  + scale_fill_gradient2(low = "#FFFFFF",mid = "#FFFFFF",high = "#FFFFFF",na.value ="#FFFFFF", midpoint = 6 ) +
  theme(axis.text.x = element_text(vjust = 0.6, hjust = 0.6, angle = 90))+
  coord_fixed(ratio=1)+
  geom_point(aes(size=P,shape=P,color=HR))+
  scale_color_manual(values = c("1" = "#d0648a", "0" = "#4ea397"))+
  scale_shape_manual(values = c("0" = 4, "1" = 16,"2" = 16, "3" = 16))+
  scale_size_manual(values = c("0" = 6,"1" = 4, "2" = 6, "3" = 8))+  
  theme(legend.position="none")  

p_Lymphoids_cell/(p_Myeloids_cell+p_Stem_cell+p_Stromal_cell+ plot_layout(widths = c(14,11,10)))/(p_Others_cell+plot_spacer()+ plot_layout(widths = c(9,26)))





#survival curve
riskresult<-cell_interact_survival
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

#scatter graph
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

#time AUC

plotAUCcurve(tROC,conf.int=F,col="blue")
plot(tROC,time=1,col='orange')
plot(tROC,time=3,col='red',add=T)
plot(tROC,time=5,col='green',add=T)
plot(tROC,time=7,col='yellow',add=T)
plot(tROC,time=9,col='blue',add=T)
legend(0.6,0.3,c(paste("AUC of 1 Year =",round(tROC$AUC,3)[1]),
                 paste("AUC of 3 Year =",round(tROC$AUC,3)[3]),
                 paste("AUC of 5 Year =",round(tROC$AUC,3)[5]),
                 paste("AUC of 7 Year =",round(tROC$AUC,3)[7]),
                 paste("AUC of 9 Year =",round(tROC$AUC,3)[9])),
                 
       x.intersp=0.6, y.intersp=1.0,
       lty= 1 ,lwd= 2,col=c('orange',"red","green","yellow","blue"),
       bty = "n",
       seg.len=1,cex=0.9)# 

#single cell survival curve
cells<-c("Platelets","Erythrocytes","MSCs","Endothelial cells","Fibroblasts","LECs","Muscle cell",
  "mv Endothelial cells","Preadipocytes","Smooth muscle cell")
result<-list()

for(i in cells){
  cell_sin<-cell_interact_survival[,c("OS","OS.time",i)]
  group<-ifelse(cell_sin[,i]>=median(cell_sin[,i]),"high","low")
  cell_sin_sur<-cbind(cell_sin,group)
  colnames(cell_sin_sur)<-c("OS","OS.time","score","group")
  cell_sin_sur$OS.time=cell_sin_sur$OS.time/30
  cell_sin_sur<-list(cell_sin_sur)
  result<-c(result,cell_sin_sur)
}
names(result)<-cells

splots <- lapply(names(result), function(g){
  i = which(names(result) == g)
  test<-as.data.frame(result[i])
  colnames(test)<-c("OS","OS.time","score","group")
  sfit<-survfit(Surv(OS.time,OS)~group,data=test)
  p = survminer::ggsurvplot(sfit, pval = TRUE, palette = c("red","darkblue"),
                            data = test, legend = c(0.8, 0.8), title = names(result)[[i]])
  p2 = p$plot + theme(plot.title = element_text(hjust = 0.5))
  return(p2)
})
patchwork::wrap_plots(splots)+patchwork::plot_layout(guides = "collect",nrow=2)



