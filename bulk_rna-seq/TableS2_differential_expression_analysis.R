
###===================================DOWNLOAD THE DATA=======================================

setwd("~/OneDrive - University of Glasgow/Bioinformatics/OverexpressionTC/MauscriptFig/")

library("DESeq2")
library("kohonen")

#### download sample plan and metadata 
plan<-read.table(file="sample_plan.txt",sep="\t",header=TRUE,row.names=6)

### Download all relevant read counts
data<-read.table(file=paste("Counts/",as.character(plan[1,2]),sep=""), sep="\t",header=FALSE, col.names=c("GeneId",as.character(row.names(plan)[1])))
for (i in 2:nrow(plan)) {
  xy<-read.table(file=paste("Counts/",as.character(plan[i,2]),sep=""), sep="\t",header=FALSE,col.names=c("GeneId",as.character(row.names(plan)[i])))
  data<-merge(data,xy, by="GeneId")
}

### And table with gene annotations

Genes<-read.table(file="PBGeneTable.txt", sep="\t",header=TRUE,quote="\"'")


###=================GENERATION OF DIFFERENTIAL EXPRESSION TABLES FOR EACH TIME POINT======================

#differential expression function for a time point

DEtime<-function(x){
  sub_plan<-plan[plan$time==x&plan$Strain=="AP2Goe",]
  sub_data<-data[,c("GeneId",row.names(sub_plan))]
  dbp<-DESeqDataSetFromMatrix(countData= sub_data, colData=sub_plan[,-1], design=formula(~rapamycin), tidy = TRUE)
  dbp<-DESeq(dbp)
  res<-results(dbp, tidy=TRUE)
  names(res)[1]<-"GeneId"
  final<-merge(	res,Genes[,c(1,2)], by="GeneId")
  return(final)
}

#Generation of object with a DE table for each timepoint
DElist<-list()

for (i in c("0h","1h","2h","4h","6h","8h","12h","18h","24h","30h","44h")){
  DElist[[i]]<-DEtime(i)}

#Extract the fold changes only into a separate table
allFC<-DElist[[1]][,c(1,8)]

for (i in 1:length(DElist)){
  allFC<-merge(allFC,DElist[[i]][,c(1,3)], by= "GeneId",suffixes=c(names(DElist)[i-1], names(DElist[i])))}

#removal of NA's
allFC<-na.omit(allFC)

###=================GENERATION OF DIFFERENTIAL EXPRESSION TABLES FOR MALE/FEMALE/ASEXUAL COMPARAISONS ======================


#male to female differntial expression:
sub_plan<-plan[plan$Strain=="WT",]
sub_data<-data[,c("GeneId",row.names(sub_plan))]
dbp<-DESeqDataSetFromMatrix(countData= sub_data, colData=sub_plan[,-1], design=formula(~time), tidy = TRUE)
dbp<-DESeq(dbp)
DElist$MvsF<- results(dbp, c("time", "male", "female"), tidy=TRUE)
names(DElist$MvsF)[1]<-"GeneId"

DElist$MvsA<- results(dbp, c("time", "male", "asexual"),tidy=TRUE)
names(DElist$MvsA)[1]<-"GeneId"

DElist$FvsA<- results(dbp, c("time", "female", "asexual"),tidy=TRUE)
names(DElist$FvsA)[1]<-"GeneId"

##add male/female exoression differences to the rapamycin induced expression table:

FullDataset<-merge(allFC,DElist$MvsF[,c(1,3)], by="GeneId",suffixes=c("44h","MvsF"))




###=========CLASSIFICATION OF GENES INTO MALE, FAMALE AND ASEXUAL ENRICHED============

AllChanges<-data.frame("GeneId"=DElist$MvsA$GeneId,"MSfc"=DElist$MvsA$log2FoldChange,
                       "MSpa"=DElist$MvsA$padj,"FSfc"=DElist$FvsA$log2FoldChange,
                       "FSpa"=DElist$FvsA$padj,"MFfc"=DElist$MvsF$log2FoldChange,
                       "MFpa"=DElist$MvsF$padj)

#add extra column indicating the gene enrochment
AllChanges$Type<-"ambigious"
AllChanges$Type[(AllChanges$MSfc>2&AllChanges$MSpa<0.01)&(AllChanges$FSfc>2&AllChanges$MSpa<0.01)]<-"gametocyte"
AllChanges$Type[(AllChanges$MSfc>2&AllChanges$MSpa<0.01)&(AllChanges$MFfc>2&AllChanges$MFpa<0.01)]<-"male"
AllChanges$Type[(AllChanges$FSfc>2&AllChanges$MSpa<0.01)&(AllChanges$MFfc<2&AllChanges$MFpa<0.01)]<-"female"
AllChanges$Type[(AllChanges$MSfc<2&AllChanges$MSpa<0.01)&(AllChanges$FSfc<2&AllChanges$MSpa<0.01)]<-"asexual"

#merge wiith the fold change data
FullDataset<-merge(FullDataset,AllChanges[,c("GeneId", "Type")], by="GeneId")


###==============EXPORT ALL DE TABLES FOR FURTHER ANALYSIS ==========================

##All Fold change tables with statistics

for (i in 1:length(DElist)){
  write.csv(DElist[[i]], file=paste(names(DElist)[i],"fold_change_table.csv",sep="_"),row.names = FALSE)
}

### and combined fold changes table for downstream analysis:

write.csv(FullDataset,file="CombinedFoldChanges.csv",row.names=FALSE)
