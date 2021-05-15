#build coexpression network, we only extract top 1% sub-network
#Then identify the module and generate module-realted list
rm(list=ls())
rm(list=ls())
rm(list=ls())

require(matrixStats)
require(data.table)
require(igraph)
require(Hmisc)
require(compiler)
require(reshape2)

setwd(" ")
source('networkson_change.R')

top_quantile=0.99 #extract the coexpressed gene pair with r value at top 1%
num_nodes=5
transitivity=0.5# a cutoff to identify module in which the genes have high Connectivity
exp_th=1

setwd(" ")
#load("your_tpm_exp.Rdata")#load symbol_exp, thus I use symbol_exp afterwards
symbol_exp <-as.matrix(read.csv(".../....txt",header=T,row.names = 1,sep="\t"))
mean_value<-as.matrix(rowMeans(symbol_exp))
index<-which(mean_value>exp_th)
symbol_exp<-symbol_exp[index,]#only keep the genes average TPM >1

path_out<-"..."
setwd(path_out)

corMatrix = makeCorTable(symbol_exp, cutoff = top_quantile, mode = "spearman", self=F, debug =F)##Creat co-expression matrix (corMatrix); take longest time
save(file="coexp_network_top1.Rdata",corMatrix)

corNet = makeCorNet(corMatrix)#Creat the network graph based on the co-expression matrix
moduleList= makeModuleList(corNet, debug = F) #Use random walk method to extract interactive module; take some time
save(file="module_list.Rdata",moduleList)

cytoscapematerial = annotateModulesByCC(corNet, moduleList, cutCluster = num_nodes, cutCC = transitivity, debug = F)#Extract the modue list with number of nodes (cutCluster) >5; If connectivity (cutCC) > 0.5, label as "HighCC", otherwise "LowCC".

write.node.cytoscape(cytoscapematerial$nodeTable, 'cytoNodefile.txt')#output the module number and its transitivity/Connectivity

############################################################33
#overlap between the genes in each module and marker genes 
#overlap significance was determined by Hypergeometric distribution
#(in your case, you marker genes should be the DEGs between AD and normal control,DEGs cutoff: fdr<0.05 or fdr<0.01)
rm(list=ls())
rm(list=ls())
rm(list=ls())

setwd("...")
source('two_gene_list_over_overlap.R')

setwd("...")
#load("your_tpm_exp.Rdata")#load symbol_exp
symbol_exp <-as.matrix(read.csv(".../...txt",header=T,row.names = 1,sep="\t"))
mean_value<-as.matrix(rowMeans(symbol_exp))
index<-which(mean_value>exp_th)
symbol_exp<-symbol_exp[index,]
bg_num<-dim(symbol_exp)[1]
#rm(clinical_exp,symbol_exp,ensem_of_symbol_exp)

setwd("...")#change to your markers,in your case, the DEGs between AD and normal;
KIRC_marker<-as.matrix(read.csv("....txt",header=T,sep="\t"))
up_marker<-as.matrix(KIRC_marker[,])# a list of gene symbol of your marker genes

path_raw<-"..."
setwd(path_raw)
load("module_list.Rdata")
module_info<-as.matrix(read.csv("cytoNodefile.txt",header=T,sep="\t"))
result_end<-NULL
for (i in 1:dim(module_info)[1]){
  num_module<-as.numeric(module_info[i,1])
  gene_list_1<-as.matrix(moduleList[num_module][[1]])
  result_each<-two_gene_list_overlap(gene_list_1,up_marker,bg_num)
  result_end<-rbind(result_end,result_each)
}
result_end<-cbind(as.matrix(module_info[,1]),result_end)
colnames(result_end)[1]<-"module_num"

path_out<-paste0(path_raw,"each_module_vs_JAM_C3_marker")#--------change
dir.create(path_out)
setwd(path_out)
write.table(result_end,file="overlap_module_genes_vs_JAM_C3_markers.txt",sep="\t",row.names=F,col.names=T,quote=F)
############################################
#calculate the degree, betweeness and closeness of each node for a subnetwork
rm(list=ls())
rm(list=ls())
rm(list=ls())

setwd("code_path\\")
source('igraph_network_features.R')

path_raw<-"network_path\\"
setwd(path_raw)
load("coexp_network_top1.Rdata")#
load("module_list.Rdata")

int_gene<-as.matrix(moduleList[253][[1]])#------------change 235 to the your module number
int_corMatrix<-corMatrix[int_gene,int_gene]

corNet<-coExpressNetwork(int_corMatrix)
features<-networkNodeAnno(corNet)

write.table(features,file="features_subnetwork.txt",sep="\t",row.names=F,col.names=T,quote=F)