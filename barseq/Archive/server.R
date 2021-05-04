library(shiny)
library(tidyverse)
library(reshape2)
library(zoo)
library(here)
here()
setwd("/Users/theo/Dropbox/Sanger/GamGithub/gametocyteBarseq")

loadAndCollapse<- function(file){
  counts<-read.csv(file, header=T,stringsAsFactors=FALSE)
  counts<- counts[counts$barcode!="no_match",]
  fcounts<- as.matrix (counts [grep("\\.1$",names(counts))] ) #forward
  rcounts<- as.matrix (counts [grep("\\.2$",names(counts))] ) #reverse
  alltcounts<-fcounts+rcounts
  colnames(alltcounts)<-gsub("\\.1", "", colnames(alltcounts) )
  rownames(alltcounts)=counts$gene
  alltcounts<-as.data.frame(t(alltcounts))
  alltcounts$Index=1:nrow(alltcounts)
  alltcounts$file=file
  row.names(alltcounts)<-NULL
  return(alltcounts)
}
files<-c("PbSTM101.csv","replicateA.csv","replicateB.csv", "counts_24160.csv")

dfs<-lapply(files,loadAndCollapse)
bigdf<-plyr::rbind.fill(dfs)
bigdf$Run=gsub(".csv","",bigdf$file)
metadata<-read.csv("SampleAssignments.csv")
####
#DETECTIVE
'
new<-bigdf %>% filter(file=="counts_24160.csv")

indexes<-read_csv("indexesmessedup.csv") %>% filter(Type!="")
join<-inner_join(new,indexes)
narrow<-melt(join,id.vars=c("Run","Type","Id","Index","Mouse"),variable.name="gene")
narrow$value<-as.numeric(as.character(narrow$value))
narrow$value[is.na(narrow$value)]=0
narrow<-narrow %>% group_by(Id,Index) %>% mutate(prop=value/sum(value))
ggplot(filter(narrow,gene=="PBANKA_051490"),aes(y=prop,x=Type,color=Id))+geom_point()
summary<-narrow %>% group_by(Type,gene) %>% summarise(meanprop=mean(prop))
top<-summary %>% filter(Type=="neg") %>% arrange(-meanprop)
narrow$Type=as.character(narrow$Type)
narrow[narrow$Type=="GFPsort",]$Type="G"
narrow[narrow$Type=="RFPsort",]$Type="R"
narrow[narrow$Type=="neg",]$Type="N"
narrow$Id=as.character(narrow$Id)
#narrow[narrow$Id=="gfp4",]$Id="A"
#narrow[narrow$Id=="rfp14",]$Id="A"
narrow$gene=factor(as.character(narrow$gene),levels=as.character(top$gene))

ggplot(filter(narrow,gene %in% top$gene[1:120]),aes(y=prop,x=Type,color=Id))+geom_point()+facet_wrap(~gene,scales="free")+expand_limits(y=0)+scale_color_brewer(palette="Set1")
ggsave("hello.pdf",width=20,height=14)

thenext<- narrow %>% group_by(Type,Id,gene,Mouse) %>% summarise(prop=mean(prop)) %>% filter(!(Id %in% c("neg1","neg2","neg5")))
ggplot(filter(thenext,gene %in% top$gene[1:180]),aes(y=prop,x=Type,color=Mouse,group=Mouse))+geom_line(alpha=0.8)+geom_point(alpha=0.8)+facet_wrap(~gene,scales="free")+expand_limits(y=0)+scale_color_brewer(palette="Set1")
ggsave("hello.pdf",width=25,height=20)


normal<- filter(thenext,gene %in% c("PBANKA_051500","PBANKA_071830")) %>% group_by(Type,Mouse) %>% summarise(prop=mean(prop))
thenext2<-inner_join(thenext,normal,by=c("Mouse","Type")) %>% ungroup() %>% mutate(prop=prop.x/prop.y)
ggplot(filter(thenext2,gene %in% top$gene[1:180]),aes(y=prop,x=Type,color=Mouse,group=Mouse))+geom_line(alpha=0.8)+geom_point(alpha=0.8)+facet_wrap(~gene,scales="free")+expand_limits(y=0)+scale_color_brewer(palette="Set1")
ggsave("hello.pdf",width=20,height=14)


tops<-thenext %>% group_by(Type,gene) %>% arrange(prop)  %>% filter (gene %in% top$gene[1:80]) %>% mutate(pord=row_number())
new<-tops %>% ungroup() %>% select(gene,pord,Id) %>%spread(Id,pord) %>% select(-gene)
heatmap(as.matrix(new))
#GFP3, RFP15
#GFP13,RFP14
#GFP4,RFP6
#neg 1,2,5 and neg 10,11,12 seem different
join=inner_join(tops,tops,by="gene") %>% filter(Type.x != Type.y)
table(join$Id.x,join$Id.y)

filt<-newmat[sums>300,]
colsu=colSums(filt)
div<-sweep(filt,2,colsu,`/`)
heatmap(div)

colnames(filt)=paste0("V",1:28)
normalgrowth=filt
normalgrowth=filt[rownames(filt) %in% c("PBANKA_143520","PBANKA_140920","PBANKA_041340"),
                  1:28 %in% c(1,2,3,4,5,6,10,11,12,13,14,15)]
colsu=colSums(normalgrowth)
div<-sweep(normalgrowth,2,colsu,`/`)
heatmap(div)

sub<-div[, !(colnames(div) %in% c("V9","V16","V7","V8"))]
abc<-prcomp(t(sub))
df2<-as.tibble(abc$x)
df2$name=colnames(sub)
ggplot(df2,aes(x=PC1,y=PC2,label=name))+geom_point()+geom_text(color="red",nudge_y=0.002)
####
'


bigdf<-merge(bigdf,metadata,by=c("Index","Run"),all.x=T)





bigdf <- bigdf  %>% filter(Condition!="DoNotAnalyse")
bigdf$Run
bigdf$Index<-NULL
bigdf$file<-NULL
bigdf$Description<-NULL
bigdf$Mouse=as.factor(bigdf$Mouse)
narrow<-melt(bigdf,id.vars=c("Pool","Condition","Mouse","Replicate","Run"),variable.name="gene")
narrow<-narrow %>% filter(!is.na(Pool))
narrow$value<-as.numeric(as.character(narrow$value))
narrow$value[is.na(narrow$value)]=0


narrow<-narrow %>% group_by(Pool,Condition,Mouse,Replicate) %>% mutate(proportion=value/sum(value))
wellrepresentedinpassage<-narrow %>% filter(Condition=="control" || Condition =="Negsort")%>%group_by(gene,Pool) %>%filter(max(proportion)>0.0001) %>% summarise(toinclude=T)
narrow<-narrow %>% left_join(wellrepresentedinpassage) %>% filter(toinclude==T)
narrow <- narrow %>% mutate(value=value+0.5)
narrow<-narrow %>% group_by(Pool,Condition,Mouse,Replicate) %>% mutate(proportion=(value)/sum(value))
narrow<- narrow %>% mutate(logprop=log2(proportion))
totals<-narrow%>% group_by(Pool,Condition,Mouse,Replicate,Run)%>% summarise(sum=sum(value))
ggplot(totals,aes(x=paste(Pool,Condition,Mouse,Replicate),fill=Pool,y=sum))+geom_bar(stat="identity")+facet_wrap(~Run,scales="free",ncol=3) +theme_bw()+ theme(axis.text.x = element_text(angle = 90, hjust = 1,size=2))


ggsave("terrible.pdf",width=8,height=5)


#
#InputGrid<-read.csv("InputGrid.csv")
#gatheredInput<-InputGrid %>% gather("Pool","MeantToBeIn",-gene) %>% filter(MeantToBeIn>0)


#Now model variance
#narrow<-merge(narrow,gatheredInput)
variance<-narrow %>% group_by(Pool,Condition,Mouse,gene) %>% summarise(sd=sd(logprop),val1=logprop[1],val2=logprop[2],meanlogprop=mean(logprop))
cors<-variance %>% group_by(Condition,Mouse,Pool) %>% summarise( cor=cor(val1,val2),n=n())
ggplot(variance,aes(x=val1,y=val2,color=Condition))+
  geom_rect(data=cors,xmin= -50,xmax=50,ymin= -50,ymax=50,aes(fill=-cor,x=0,y=0))+
  geom_point(alpha=0.4,size=0.1)+facet_grid(Pool~Condition+Mouse)+scale_fill_distiller(palette="Greys")+theme( panel.background = element_rect(fill="white")) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+scale_color_brewer(palette="Set1")+
  theme(axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank())

#ggplot(variance,aes(x=val1,y=val2,color=Condition))+
#  geom_point(alpha=0.4,size=0.1)+facet_grid(Pool~Condition+Mouse)+theme( panel.background = element_rect(fill="white")) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+scale_color_brewer(palette="Set1")+
#  theme(axis.title=element_blank(),
 #       axis.text=element_blank(),
  #      axis.ticks=element_blank())



#ggplot(variance,aes(x=meanlogprop,y=abs(sd),color=Condition,linetype=Mouse))+geom_smooth()+facet_grid(Pool~Mouse)



#variance<-filter(variance, !(Pool=="H" & Mouse==3))
variance$Condition=as.character(variance$Condition)
variance<-variance %>% ungroup() %>% mutate(Condition=ifelse(Pool=="SP3" & Condition=="control","Negsort",Condition))

variance<-variance %>% group_by(Pool,Condition,Mouse) %>% arrange(-meanlogprop) %>% mutate(medsd=rollmean(x=sd,k=11, fill="extend"))
variance<- variance %>% mutate(monotonicmedsd=cummax(medsd))

#ggplot(variance,aes(x=meanlogprop,y=monotonicmedsd,color=Condition,linetype=Mouse))+geom_point()+facet_wrap(~Pool)


#ggplot(variance,aes(x=meanlogprop,y=monotonicmedsd,color=Condition))+
#  geom_line()+facet_grid(Pool~Condition+Mouse)+theme( panel.background = element_rect(fill="white")) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+scale_color_brewer(palette="Set1")+
 # theme(axis.title=element_blank(),
  #      axis.text=element_blank(),
   #     axis.ticks=element_blank())

noncontrol<-variance %>% filter(Condition!="Negsort")
control<-variance %>% filter(Condition=="Negsort")

comparison<-inner_join(ungroup(control),ungroup(noncontrol),by=c("Mouse","gene","Pool"))

comparison<- comparison %>% mutate(difference=meanlogprop.y-meanlogprop.x)
comparison<- comparison %>% mutate(differencesd=sqrt(monotonicmedsd.y^2+monotonicmedsd.x^2))

comparison<- comparison %>% mutate(Condition=Condition.y)


gaussianMeanAndVarianceSplit<- function(vals, variances,return){
  df<-data.frame(value=vals,variance=variances)
  df<-df[complete.cases(df),]
  
  vals=df$value 
  variances=df$variance
  if(length(vals)==1){
    var<-variances[1]
    mean<-vals[1]
  }
  else{
    
    precs=1/variances
    
    mean=sum(vals*precs)/sum(precs)
    var1=(1/sum(precs))*(1/(length(vals)-1))*sum((vals-mean)**2/variances)
    var2=1/sum(precs) 
    if(is.na(var1)){var1<-0}
    var<-max(var1,var2)
  }
  if (return=="mean"){
    return(mean)
  }
  if (return=="variance"){
    return(var)
  }
}



#comparisonmerge<- comparison %>% group_by(Condition,gene,Pool) %>% summarise(difference=gaussianMeanAndVarianceSplit(difference,differencesd^2,"mean"),differencesd=sqrt(gaussianMeanAndVarianceSplit(difference,differencesd^2,"variance")))
comparisonmerge<- comparison %>% group_by(Condition,gene) %>% summarise(difference=gaussianMeanAndVarianceSplit(difference,differencesd^2,"mean"),differencesd=sqrt(gaussianMeanAndVarianceSplit(difference,differencesd^2,"variance")))

comparisonmerge<- comparisonmerge %>% mutate(p=1-pnorm(0,mean=difference,sd=differencesd))
comparisonmerge<- comparisonmerge %>% mutate(diffmax=difference+2*differencesd)
comparisonmerge<- comparisonmerge %>% mutate(diffmin=difference-2*differencesd)
comparisonmerge<- comparisonmerge %>% mutate(power=ifelse(diffmin>-1,"notreduced",ifelse(diffmax<(-1),"reduced","nopower")))
GFP<-filter(comparisonmerge,Condition=="GFPsort")
GFP$gene=as.character(GFP$gene)
RFP<-filter(comparisonmerge,Condition=="RFPsort")
RFP$gene=as.character(RFP$gene)
#joint=full_join(GFP,RFP,by=c("gene","Pool"))
joint=full_join(GFP,RFP,by=c("gene"),suffix=c(".gfp",".rfp"))
joint$minmax=pmin(joint$diffmax.gfp,joint$diffmax.rfp)
joint <- mutate(joint, 
                class=case_when(
                  diffmax.gfp < -2 &  diffmax.rfp < -2  ~ "lossOfBoth",
                  diffmax.gfp < -2   ~ "lossOfMales",
                  diffmax.rfp < -2   ~ "lossOfFemales",
                  diffmin.rfp >1 &  diffmin.gfp >1  ~ "gainOfBoth",
                  diffmin.rfp >1   ~ "gainOfFemales",
                  diffmin.gfp >1   ~ "gainOfMales",
                  diffmax.gfp > -1 & diffmax.rfp > -1 ~ "No change",
                  TRUE ~ "Unsure"
                ))
write.csv(joint,"joint.csv")
candidates<-c("PBANKA_123760","PBANKA_130270","PBANKA_082800","PBANKA_090240","PBANKA_141810","PBANKA_143520","PBANKA_071230","PBANKA_145480","PBANKA_010240","PBANKA_041340","PBANKA_071650","PBANKA_144790","PBANKA_142770","PBANKA_143750")





comparison<- comparison %>% mutate(diffmax=difference+2*differencesd)
comparison<- comparison%>% mutate(diffmin=difference-2*differencesd)

comparisonmerge$logp=log10(comparisonmerge$p+0.0000001)


phenolevels=c("Insufficient data","Essential","Slow","Dispensable","Fast","Unselected")  

phenolevelscolor=c("black", "#f90f00", "#0f3791","#007e41","#f83cd9", "darkgray") #standard phenotype colours, colour-blind safe
othercolors=c("#66c2a5","#fc8d62","#8da0cb", "white") #used when another palette is needed
names(phenolevelscolor)=phenolevels

# Define server logic required to draw a histogram
shinyServer(function(input, output) {
   
  
  
  selected<-reactive({
    selpoints<-filter(thetable(),  selected==T)
    if(nrow(selpoints)==1){
      return( selpoints[1, ])
    }
    else{
      return(  nearPoints(thetable(), input$topplot_click)[1, ])
    }
  })
  
  thetable<-reactive({
    if(input$var=="diffmin"){
    spread<-comparisonmerge %>% select(gene,Condition,difference)%>% group_by(gene) %>% spread(key=Condition,value=difference)
    }
    else{
      spread<-comparisonmerge %>% select(gene,Condition,diffmax)%>% group_by(gene) %>% spread(key=Condition,value=diffmax)
      
    }
    spread$incand=as.character(spread$gene) %in% as.character(candidates)
    ggplot(spread,aes(x=GFPsort,y=RFPsort,color=incand))+geom_point()
    
    geneinfo<-read.csv("geneinfo.csv")
    spreaded<-spread
    spreaded<-merge(spreaded, geneinfo,by.y="Old.Gene.ID",by.x="gene")
    bloodstage<-read.csv("bloodstage.csv")
    joint<-merge(joint,bloodstage,by="gene")
    
    
    spreaded<-merge(spreaded,select(bloodstage,gene,phenotype,Relative.Growth.Rate))
    
    newspreaded<-spreaded
    if(input$slow ==F){newspreaded <- filter(newspreaded,phenotype!="Slow" & phenotype!="Essential")}
    if(trimws(input$gene)!="")
    {
    newspreaded$selected=grepl(trimws(input$gene),newspreaded$gene,ignore.case=T) | grepl(trimws(input$gene),newspreaded$current_version_ID,ignore.case=T) | grepl(trimws(input$gene),newspreaded$gene_product,ignore.case=T) | grepl(trimws(input$gene),newspreaded$gene_name,ignore.case=T)
    abc<-strsplit(input$gene, ",",fixed=T)
    if(length(abc)>1){
    newspreaded$selected=as.character(newspreaded$gene) %in% abc
    }
   # newspreaded$selected=newspreaded$gene %in% c("PBANKA_041340","PBANKA_145480","PBANKA_090230","PBANKA_090240","PBANKA_143520","PBANKA_142770","PBANKA_141770","PBANKA_141810","PBANKA_082800","PBANKA_010240","PBANKA_140300","PBANKA_071650","PBANKA_121260","PBANKA_134040","PBANKA_135990","PBANKA_082460","PBANKA_092930","PBANKA_010570","PBANKA_144960","PBANKA_140890","PBANKA_082660","PBANKA_144280","PBANKA_100120","PBANKA_130180","PBANKA_092920","PBANKA_140700","PBANKA_136040","PBANKA_144790","PBANKA_101870","PBANKA_145500","PBANKA_092570","PBANKA_113420","PBANKA_060330","PBANKA_143430","PBANKA_113320","PBANKA_092160","PBANKA_081940","PBANKA_061600","PBANKA_083370","PBANKA_123580","PBANKA_080840","PBANKA_070790","PBANKA_090950","PBANKA_051300","PBANKA_131830","PBANKA_082190","PBANKA_143000","PBANKA_112130","PBANKA_141330","PBANKA_050790","PBANKA_121520","PBANKA_132890","PBANKA_081610","PBANKA_040520","PBANKA_083090","PBANKA_112490","PBANKA_082630","PBANKA_142535","berg07_28S","PBANKA_134010","PBANKA_120210","PBANKA_114320","PBANKA_144410","PBANKA_142220","PBANKA_061610","PBANKA_010380","PBANKA_140590","PBANKA_112300","PBANKA_112290","PBANKA_124660","PBANKA_124400","PBANKA_133680","PBANKA_082650","PBANKA_123760","PBANKA_140840","PBANKA_071350","PBANKA_140880","PBANKA_101980","PBANKA_031270","PBANKA_040650","PBANKA_121400","PBANKA_093540","PBANKA_061340","PBANKA_131980","PBANKA_091990","PBANKA_142720","PBANKA_123570","PBANKA_094220","PBANKA_102340","PBANKA_010130","PBANKA_051500","PBANKA_021400","PBANKA_082790","PBANKA_071190","PBANKA_131840","PBANKA_142900","PBANKA_093180","PBANKA_132640","PBANKA_030260","PBANKA_120910","PBANKA_132290","PBANKA_101130","PBANKA_092810","PBANKA_083230","PBANKA_051550","PBANKA_021580","PBANKA_010370","PBANKA_010230","PBANKA_100990","PBANKA_092630","PBANKA_071090","PBANKA_101000","PBANKA_120480","PBANKA_142550")
    }
    else{
      newspreaded$selected=F;
      }
    
    newspreaded
    
    
    })
  output$tbl = DT::renderDataTable(
    {
      
      thetable()
    
    },selection = 'single', rownames= FALSE)
  
 # output$gene_name = renderText({
  #  selected()$gene
#  })
  output$basicplot = renderPlot({

    g=selected()$gene
    req(g)
    number_ticks <- function(n) {function(limits) pretty(limits, n)}
    ggplot(ungroup(filter(variance,gene==g,Condition%in%c("control","Negsort","GFPsort","RFPsort"))),aes(x=Condition,color=Mouse,y=2^(meanlogprop-0.5),group=Mouse))+geom_point(position=position_dodge(width=0.2))+geom_line()+ggtitle(g)+scale_y_log10(breaks=number_ticks(4))+labs(y="Proportion",x="Population")+geom_errorbar(aes(ymin=2^(meanlogprop-2*monotonicmedsd-.5),ymax=2^(meanlogprop+2*monotonicmedsd-.5)),width=0.2,position="dodge",alpha=0.5)+facet_wrap(~Pool, scales="free")+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
  },res=120)
  output$topplot = renderPlot({
    
    p= ggplot(thetable(),aes(x=GFPsort,y=RFPsort,color=phenotype))+geom_point(size=0.7)+scale_color_manual(values = phenolevelscolor)
    if(input$highlight){
      p= p+geom_point(data=filter(thetable(),gene%in% candidates),shape=1,size=3,color="black")
    }
     p= p+geom_point(data=filter(thetable(),  selected==T)  ,shape=1,size=8,color="red",stroke=2)
    
    p +      geom_vline(xintercept =0)+      geom_hline(yintercept =0)+labs(x="GFP (males)",y="RFP (females)")
    },res=80)
  output$mergedeffectplot = renderPlot({
    
    g=selected()$gene
    req(g)
    number_ticks <- function(n) {function(limits) pretty(limits, n)}
    ggplot(ungroup(filter(comparisonmerge,gene==g)),aes(x=Condition,y=difference,ymin=diffmin,ymax=diffmax,fill=Condition))+geom_bar(stat="identity")+geom_errorbar(width=0.2,alpha=0.5)+scale_fill_brewer(palette="Accent")+labs(x="Population",y="Effect",alpha=0.5)+geom_hline(yintercept=0)+ theme(axis.text.x = element_text(angle = 90, hjust = 1))+ggtitle(g)
    
    },res=120)
  output$gene=renderText({as.character(selected()$gene)})
  output$geneproduct=renderText({paste(as.character(selected()$gene_name),": ",as.character(selected()$gene_product))})
})


'
SP3=filter(comparisonmerge,Pool=="SP3")
Targeted=filter(comparisonmerge,Pool=="Targeted")
#Targeted$Condition=ifelse(Targeted$Condition=="GFPsort","RFPsort","GFPsort")
scatter<-inner_join(SP3,Targeted,by=c("Condition","gene"))
ggplot(scatter,aes(x=difference.x,y=difference.y,color=Condition))+geom_point()+geom_smooth(method="lm")

SP3=filter(comparisonmerge,Pool=="SP3")
Targeted=filter(comparisonmerge,Pool=="PbSTM101")
#Targeted$Condition=ifelse(Targeted$Condition=="GFPsort","RFPsort","GFPsort")
scatter<-inner_join(SP3,Targeted,by=c("Condition","gene"))
ggplot(scatter,aes(x=difference.x,y=difference.y,color=Condition))+geom_point()+geom_smooth(method="lm")



ggplot(joint,aes(x=diffmax.x,y=diffmax.y))+geom_point(alpha=0.1)+facet_wrap(~Pool)+geom_point(color="red",data=filter(joint,gene=="PBANKA_102620"))
'