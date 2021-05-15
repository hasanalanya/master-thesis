rm(list=ls())
rm(list=ls())
rm(list=ls())
library(readxl)
library(tidyverse)
#Based on the name of the cohort.
country = paste("Singapore" , "C")
DiffExpGeneSing <- read_excel('...xlsx', sheet = 1)
DiffExpGeneSing_f <- DiffExpGeneSing %>% filter(gene_biotype == 'protein_coding') %>% mutate(label = case_when(
  (adjp_sarc < 0.1 & coef_sarc < 0) ~ "2", (adjp_sarc < 0.1 & coef_sarc > 0) ~ "1", TRUE~ "3" )) %>% select(Gene, label, everything())

pos_jama <- read.csv('.../Singapore_list.txt', sep = "\t")
# Summary
num2cluster <- as.data.frame(table(pos_jama$modulist)) 
colnames(num2cluster) <-c("Cluster", "Number of Genes")

# x = Number of genes overlapped in cluster 1 and cluster 2
# m = Number of genes in cluster 2
# n = Number of genes Total - m
# k = Number of genes in cluster 1
# Total = 19983 in this case

res = matrix(1:2*length(unique(pos_jama$modulist)), nrow = 2, ncol = length(unique(pos_jama$modulist)))
res_p = matrix(1:2*length(unique(pos_jama$modulist)), nrow = 2, ncol = length(unique(pos_jama$modulist)))


up=list()
down = list()

for(i in 1:2){
  # Cluster in GTEX
  print(paste0(i, ' in project DiffExpSing'))
  for(j in 1:length(unique(pos_jama$modulist))){
    print(paste0(j, ' cluster in project jama'))
    N = intersect(DiffExpGeneSing_f$Gene, pos_jama$Gene)
    cluster2 = as.character(pos_jama %>% filter(modulist == unique(pos_jama$modulist)[j]) %>% select(Gene) %>% pull(Gene))
    cluster1 = DiffExpGeneSing_f %>% filter(label == i) %>% select(Gene) %>% pull(Gene)
    x = intersect(cluster1, cluster2)
    m = intersect(N, cluster1)
    n = intersect(N, cluster2)
    
    p.value <-  phyper(length(x) -1, length(m), length(N)-length(m), length(n), lower.tail=FALSE)
    res[i,j] = p.value
    if(p.value < 0.05){
      res_p[i,j] = "*"
    }else{
      res_p[i,j] = ""
    }
    
    # Extracting the number of genes which are upregulated or downregulated in each cluster
    if(i ==1){
      up = append(up, length(x))
    }
    
    if(i ==2){
      down = append(down, length(x))
    }
  }
}

colnames = list()
for(o in unique(pos_jama$modulist)){
  colnames = append(colnames,paste0(country,o))
}


num2cluster$Up = up
num2cluster$Down = down



rownames(res) = c("UP","Down")
colnames(res) = colnames


my_palette <- colorRampPalette(c("cadetblue4", "cadetblue3", "cadetblue2"))(n = 100)
res_test <- heatmap.2(res, cellnote = res_p,notecol="black", notecex=1,density.info="none", trace="none",
                      margins = c(6,6), col = my_palette, cexRow = 0.8
                      , cexCol = 0.8,srtCol=45,  adjCol = c(1,0.6))

save.image('Res_ht_diff_Vs_Sing.Rdata')

