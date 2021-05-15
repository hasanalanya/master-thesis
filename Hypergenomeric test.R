rm(list=ls())
rm(list=ls())
rm(list=ls())

pos_GTEX <- read.csv('.../.txt', sep = "\t")
pos_Jamaica <- read.csv('.../txt', sep = "\t")

# x = Number of genes overlapped in cluster 1 and cluster 2
# m = Number of genes in cluster 2
# n = Number of genes Total - m
# k = Number of genes in cluster 1
# Total = 19983 in this case

res = matrix(1:length(unique(pos_GTEX$modulist))*length(unique(pos_Jamaica$modulist)), nrow = length(unique(pos_GTEX$modulist)), ncol = length(unique(pos_Jamaica$modulist)))
res_p = matrix(1:length(unique(pos_GTEX$modulist))*length(unique(pos_Jamaica$modulist)), nrow = length(unique(pos_GTEX$modulist)), ncol = length(unique(pos_Jamaica$modulist)))


for(i in 1:length(unique(pos_GTEX$modulist))){
  # Cluster in GTEX
  print(paste0(i, ' cluster in project GTEX'))
  for(j in 1:length(unique(pos_Jamaica$modulist))){
    print(paste0(j, ' cluster in project Jamaica'))
    N = intersect(pos_GTEX$Gene, pos_Jamaica$Gene)
    cluster2 = pos_Jamaica %>% filter(modulist == unique(pos_Jamaica$modulist)[j]) %>% select(Gene) %>% pull(Gene)
    cluster1 = pos_GTEX %>% filter(modulist == unique(pos_GTEX$modulist)[i]) %>% select(Gene) %>% pull(Gene)
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
  }
}

rownames = list()
for(i in unique(pos_GTEX$modulist)){
  rownames = append(rownames,paste0("GTEX_C",i))
}

colnames = list()
for(i in unique(pos_Jamaica$modulist)){
  colnames = append(colnames,paste0("Jamaica_C",i))
}

rownames(res) = rownames
colnames(res) = colnames


my_palette <- colorRampPalette(c("cadetblue4", "cadetblue3", "cadetblue2"))(n = 100)
res_test <- heatmap.2(res, cellnote = res_p,notecol="black", notecex=1,density.info="none", trace="none",
                      margins = c(10,10), col = my_palette, cexRow = 0.8
                      , cexCol = 0.8,srtCol=45,  adjCol = c(1,0.6))
save.image('Res_ht_gtexvsCohortJamaica_update.Rdata')

