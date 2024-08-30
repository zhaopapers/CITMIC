#SDY144 
GSE52005<-read.delim("GSE52005_symbol_aggregate_20000.csv")
library(immunedeconv)
tcga_quantiseq<-deconvolute(as.matrix(GSE52005),"quantiseq")
tcga_mcp<-deconvolute(GSE52005,"mcp_counter")
tcga_epic<-deconvolute(GSE52005,"epic")
deco_xcell<-deconvolute(GSE52005,"xcell")
tcga_abis<-deconvolute(GSE52005,"abis")
lnscore_GSE52005<-lnscore(GSE52005)

load()
