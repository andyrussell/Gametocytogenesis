setwd("~/OneDrive\ -\ University\ of\ Glasgow/Bioinformatics/OverexpressionTC/MauscriptFig")

library(kohonen)

library(viridis)

#####============== DATASET DOWNLOAD ===================================

FoldChangeTable<-read.table(file="CombinedFoldChangesSOM.txt",sep="\t",header=TRUE, row.names = 1)

Genes<-read.table(file="PBGeneTable.txt", sep="\t",header=TRUE,quote="\"'")



#####============== CREATION OF SOM MODEL ===================================

#create a grid,this one will give 64 clusters 8x8

som_grid <- somgrid(xdim = 8, ydim=8, topo="hexagonal") 

#prepare the fold change data as matrix object

som_data<-as.matrix(na.omit(FoldChangeTable[,c(1:12)]))

#run som clustering

set.seed(2560)
som_model <- som(som_data, 
                 grid=som_grid, 
                 rlen=200,
                 alpha=c(0.05,0.01),
                 n.hood = "circular",
                 keep.data = TRUE )

Output<-as.data.frame(cbind("GeneId"=rownames(som_model$data),som_model$data, "Cluster"=som_model$unit.classif))

#gene IDs added

Output <-merge(Output, Genes[,c(1,3)], by="GeneId")

write.csv(moze, file="SOMClusteringTable.csv")



###============ CLUSTER VISUALISATION =====================================

pdf(file="SOMtimeline.pdf", height =10, width =7)

par(mfrow =c(3,2))

for (i in 1: length(colnames(som_model$data))){
  file_name<- colnames(som_model$data)[i]
  plot(som_model, type = "property", property = som_model$codes[,i],palette.name=magma, main =file_name)
}

dev.off()

#overlaying male/female/asexual genes on the dataset
male_genes<-as.data.frame(som_data[row.names(FoldChangeTable[FoldChangeTable$Type=="male",]),])
female_genes<-as.data.frame(som_data[row.names(FoldChangeTable[FoldChangeTable$Type=="female",]),])
asexual_genes<-as.data.frame(som_data[row.names(FoldChangeTable[FoldChangeTable$Type=="asexual",]),])
gametocyte_genes<-as.data.frame(som_data[row.names(FoldChangeTable[FoldChangeTable$Type=="gametocyte",]),])

male_genes$colo<-"#016c00"
female_genes$colo<-"#a52b1e"
asexual_genes$colo<-"#0052c5"
gametocyte_genes$colo<-"#8B008B"

all_annotations<-rbind(male_genes,female_genes,asexual_genes,gametocyte_genes)

pdf(file="MFA_genes_SOM.pdf", height =8, width =8)

par(mfrow =c(1,1))

mappp<-map(som_model,as.matrix(all_annotations[,c(1:12)]))

plot(som_model, type = "mapping", main="Overlays", pch = 16, cex=0.8,classif = mappp, col=all_annotations$colo)

dev.off()

pdf(file="MFA_genes_SOM.pdf", height =8, width =8)

par(mfrow =c(1,1))

mappp<-map(som_model,as.matrix(all_annotations[,c(1:12)]))

plot(som_model, type = "mapping", main="Overlays", pch = 1, cex=0.3,classif = mappp, col=all_annotations$colo)

legend(x=1,y=0.2,
       legend=c("male", "female", "gametocyte", "asexual"), 
       col=c("#016c00","#a52b1e","#8B008B","#0052c5"), pch=16, ncol=4, bty="n")

dev.off()

pdf(file="Profiles_SOM.pdf", height =8, width =8)

plot(som_model, type = "codes",codeRendering="lines",palette.name=coolBlueHotRed)

dev.off()

