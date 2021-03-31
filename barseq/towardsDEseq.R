library(tidyverse)
counts<-read_tsv("/Users/theo/gametocytes/gametocyteBarseq/ap2gtimecourse/counts.txt",col_names = T)

phenodata<-read_tsv("/Users/theo/gametocytes/gametocyteBarseq/ap2gtimecourse/phenodata.txt",col_names = T)
phenodata$X1=gsub("X820", "820", phenodata$X1)

narrow<-counts %>% gather("sample","count",-X1)
colnames(narrow)=c("gene","X1","count")
narrow<-narrow %>% filter(count>0)

both<-inner_join(phenodata,narrow)
both<- both %>% group_by(X1) %>% mutate(proportion=count/sum(count))




library(DESeq)
cds = newCountDataSet( counts, condition )
counts2<-as.matrix(select(counts, -X1))
rownames(counts2)=counts$X1


timepoints=unique(phenodata$time)
i=9
cols=colnames(counts2)

performDE <- function(filter){
abc=tibble()

  filter2 <- filter %>% filter(Strain == "AP2Goe")
  print(unique(filter$time))
  if (length(unique(filter2$rapamycin))>1)
  {
    counttable=counts2[,base::match(filter2$X1,cols)]
    condition=filter2$rapamycin
    cds = newCountDataSet( counttable, condition )
    cds = estimateSizeFactors( cds )
    cds = estimateDispersions( cds )
    res = nbinomTest( cds, "YES", "NO" )
    abc<-as.tibble(res) %>% arrange(pval)
    abc$timepoint=unique(filter$time)
  }
  return(abc);
}




new<- phenodata %>% group_by(time) %>% filter(time!="3h") %>% do(performDE(.))


oneversion<-joint
asexual<-read_csv("/Users/theo/gametocytes/gametocyteBarseq/otherdata/asexualscreendata.csv",col_names = T)
merge<-right_join(oneversion,asexual)
new$id=gsub(".1$","",new$id)
new$time2=as.numeric(gsub("h","",new$time))
merge2<-inner_join(new,merge,by=c("id"="current_version_ID"))
filt <- merge2 %>% group_by(timepoint) %>% arrange(pval) %>% filter(row_number()<50) %>% mutate(n=row_number())
ggplot(filt ,aes(y=-n,x=time2,fill=class))+geom_tile()+geom_line(data=filter(filt,!is.na(class)),aes(group=id,color=class,alpha=0.5))


filt <- merge2 %>% group_by(timepoint) %>% arrange(log2FoldChange) %>% filter(row_number()<25) %>% mutate(n=row_number())
pal<-c("white","red","blue","green","gray","gray")
ggplot(filt ,aes(y=-n,x=time2,fill=class))+geom_tile()+geom_line(data=filter(filt,!is.na(class),class!="No change"),aes(group=id,color=class,alpha=0.5))+labs(y="Genes ranked at each timepoint by log2-fold change",x="Timepoint (hours)")+scale_fill_manual(values=pal)+scale_color_manual(values=pal)

pal<-c("red","blue","green","gray","gray")
filt <- merge2 %>% group_by(timepoint) %>% arrange(pval) %>% filter(row_number()<2500) %>% mutate(n=row_number())
library(viridis)
ggplot(filt ,aes(y=-n,x=time2,fill=minmax))+geom_tile()+geom_line(data=filter(filt,!is.na(class),class!="No change"),aes(group=id,color=class,alpha=0.5))+labs(y="Genes ranked at each timepoint by p-value",x="Timepoint (hours)")+scale_fill_viridis()



ggplot(filt)

cds = newCountDataSet( counts2, phenodata$ )

gois=c("PBANKA_0519800.1")
both<-inner_join(both,newmerge,by=c("gene"="transcript")) 
both <- mutate(both,label=paste0(gene.y,class))
bothavg<-both %>% group_by(timew,rapamycin,gene,Strain,label) %>% summarise(proportion=mean(proportion))

ggplot(both %>% filter(gene %in% gois,Strain=="AP2Goe",timew!="schizont",timew!="gams"),aes(color=rapamycin,y=proportion,x=as.numeric(timew)))+geom_point(alpha=0.2)+facet_wrap(~label,scales="free_y")+geom_line(aes(group=rapamycin),data=bothavg%>% filter(gene%in% gois,Strain=="AP2Goe",timew!="schizont",timew!="gams"))



newmerge<-filter(merge,class %in% c("lossOfMales","lossOfFemales","lossOfBoth")) %>% mutate(transcript=paste0(current_version_ID,".1"))



gois=newmerge$current_version_ID
sextranscriptomes<-read_tsv("/Users/theo/gametocytes/gametocyteBarseq/otherdata/sextranscriptomes.tsv")
sextranscriptomes<- sextranscriptomes %>% gather("typerep","count",-GeneID) %>% separate(typerep,c('type','rep'),sep=-2) %>% filter(count>0)
sextranscriptomes <- sextranscriptomes %>% group_by(type,rep)%>% mutate(proportion=count/sum(count))

ggplot(sextranscriptomes %>% filter(GeneID %in% gois),aes(x=type,y=proportion,color=type))+geom_point(alpha=0.5)+facet_wrap(~GeneID,scales="free_y")+ expand_limits(y = 0)
