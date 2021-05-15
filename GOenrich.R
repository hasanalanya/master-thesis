rm(list=ls())
rm(list=ls())
rm(list=ls())

library(clusterProfiler)
library(dplyr)
library(tidyr)
library(DOSE)
library(GO.db)
library(org.Hs.eg.db)
library(GSEABase)
num_gene<-10 
setwd(" ")
source('deg_GoTerm_clusterProfiler.R')
setwd(" ")
for (i in 1:6479){
  deg=moduleList[i][[1]]
  index<-length(deg)
  if (index>num_gene){
    pathway=enrichGO(gene=deg,OrgDb='org.Hs.eg.db',keyType = "ENSEMBL",ont="BP", pAdjustMethod = "BH",qvalueCutoff=0.05)
    path_table<-deg_GoTerm_clusterProfiler(pathway)
    write.table(path_table,file=paste0("path_table_module_",i,".txt"),sep="\t",row.names=F,col.names=T,quote=F) 
    
  }
}

################################################################3############3####################################
#for (i in 1:200){
  #deg=moduleList[i][[1]]
  #path_table<-deg
  #write.table(path_table,file=paste0("path_table_module_",i,".txt"),sep="\t",row.names=F,col.names=T,quote=F) 
#}
##################################################################################################################