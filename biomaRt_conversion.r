
list<- as.matrix(read.csv(" ",header=T,sep="\t"))

setwd(" ")

map_inform<-as.matrix(read.csv("mart_export.txt",header=T,sep="\t"))
rownames(map_inform)<-map_inform[,1]
over_ensg<-as.matrix(intersect(list,map_inform[,1]))
result<-as.matrix(map_inform[over_ensg,2])

write.table(result, file = " .txt", sep = "\t",row.names = F, col.names = F)
