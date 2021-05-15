#coexpression network, no weitht network, random walk find high connectivity module
rm(list=ls())
rm(list=ls())
rm(list=ls())

require(matrixStats)
require(data.table)
require(igraph)
require(Hmisc)
require(compiler)
require(reshape2)

quantile_cutoff=0.99
num_nodes=5
connectivity=0.5

setwd("...")
source('networkson_change.R')

path_raw<-"..."
setwd(path_raw)
exp<-as.matrix(read.csv(".../...txt",header=T,row.names = 1,sep="\t"))

mean_vector = as.matrix(rowMeans(exp))
dropindex = which(mean_vector<1)
exp= exp[-dropindex, ]

#exp<-t(exp)

corMatrix = makeCorTable(exp, cutoff= quantile_cutoff, mode="spearman", self=F, debug= F)#Creat co-expression matrix (corMatrix)
corNet = makeCorNet(corMatrix)#Create network graph based on the co-expression matrix
moduleList= makeModuleList(corNet, debug = F)#Use random walk method to extract interactive module
save(file="module_list.Rdata",moduleList)
cytoscapematerial = annotateModulesByCC(corNet, moduleList, cutCluster = num_nodes, cutCC = connectivity, debug = F)#Extract the module list with number of nodes (cutCluster) >5; If connectivity (cutCC) > 0.5, label as "HighCC", otherwise "LowCC".
write.node.cytoscape(cytoscapematerial$nodeTable, 'Output.txt')#Export

