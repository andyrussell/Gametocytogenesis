library(shiny)
library(tidyverse)
library(reshape2)
library(zoo)
#setwd("C:/Users/Theo/Dropbox/Sanger/CombinedGametocytes")

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
files<-c("PbSTM101.csv","replicateA.csv","replicateB.csv")

dfs<-lapply(files,loadAndCollapse)
bigdf<-plyr::rbind.fill(dfs)
bigdf$Run=gsub(".csv","",bigdf$file)
metadata<-read.csv("SampleAssignments.csv")
bigdf<-merge(bigdf,metadata,by=c("Index","Run"),all.x=T)
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
wellrepresentedinpassage<-narrow %>% filter(Condition=="control")%>%group_by(gene,Pool) %>%filter(max(proportion)>0.0001) %>% summarise(toinclude=T)
narrow<-narrow %>% left_join(wellrepresentedinpassage) %>% filter(toinclude==T)
narrow <- narrow %>% mutate(value=value+0.5)
narrow<-narrow %>% group_by(Pool,Condition,Mouse,Replicate) %>% mutate(proportion=(value)/sum(value))
narrow<- narrow %>% mutate(logprop=log2(proportion))
totals<-narrow%>% group_by(Pool,Condition,Mouse,Replicate,Run)%>% summarise(sum=sum(value))
ggplot(totals,aes(x=paste(Pool,Condition,Mouse,Replicate),fill=Pool,y=sum))+geom_bar(stat="identity")+facet_wrap(~Run,scales="free",ncol=3) +theme_bw()+ theme(axis.text.x = element_text(angle = 90, hjust = 1,size=2))
ggsave("terrible.pdf",width=8,height=5)


forcorrelating<-filter(ungroup(narrow),Pool=="PbSTM101") %>% mutate(thing=paste(Condition,Mouse,Replicate)) %>% select(value,thing,gene) %>% spread(key=thing,value=value) %>% select(-1)

mat<-as.matrix(forcorrelating)
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


variance<-variance %>% group_by(Pool,Condition,Mouse) %>% arrange(-meanlogprop) %>% mutate(medsd=rollmean(x=sd,k=11, fill="extend"))
variance<- variance %>% mutate(monotonicmedsd=cummax(medsd))

#ggplot(variance,aes(x=meanlogprop,y=monotonicmedsd,color=Condition,linetype=Mouse))+geom_point()+facet_wrap(~Pool)


#ggplot(variance,aes(x=meanlogprop,y=monotonicmedsd,color=Condition))+
#  geom_line()+facet_grid(Pool~Condition+Mouse)+theme( panel.background = element_rect(fill="white")) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+scale_color_brewer(palette="Set1")+
 # theme(axis.title=element_blank(),
  #      axis.text=element_blank(),
   #     axis.ticks=element_blank())

noncontrol<-variance %>% filter(Condition!="control")
control<-variance %>% filter(Condition=="control")

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



comparisonmerge<- comparison %>% group_by(Condition,gene) %>% summarise(difference=gaussianMeanAndVarianceSplit(difference,differencesd^2,"mean"),differencesd=sqrt(gaussianMeanAndVarianceSplit(difference,differencesd^2,"variance")))

comparisonmerge<- comparisonmerge %>% mutate(p=1-pnorm(0,mean=difference,sd=differencesd))
comparisonmerge<- comparisonmerge %>% mutate(diffmax=difference+2*differencesd)
comparisonmerge<- comparisonmerge %>% mutate(diffmin=difference-2*differencesd)
comparisonmerge<- comparisonmerge %>% mutate(power=ifelse(diffmin>-1,"notreduced",ifelse(diffmax<(-1),"reduced","nopower")))
GFP<-filter(comparisonmerge,Condition=="GFPsort")
GFP$gene=as.character(GFP$gene)
RFP<-filter(comparisonmerge,Condition=="RFPsort")
RFP$gene=as.character(RFP$gene)
joint=full_join(GFP,RFP,by="gene")
joint$minmax=pmin(joint$diffmax.x,joint$diffmax.y)
write.csv(joint,"joint.csv")
candidates<-c("PBANKA_123760","PBANKA_130270","PBANKA_082800","PBANKA_090240","PBANKA_141810","PBANKA_143520","PBANKA_071230","PBANKA_145480","PBANKA_010240","PBANKA_041340","PBANKA_071650","PBANKA_144790","PBANKA_142770","PBANKA_143750")





comparison<- comparison %>% mutate(diffmax=difference+2*differencesd)
comparison<- comparison%>% mutate(diffmin=difference-2*differencesd)

comparisonmerge$logp=log10(comparisonmerge$p+0.0000001)
spread<-comparisonmerge %>% select(gene,Condition,diffmax)%>% group_by(gene) %>% spread(key=Condition,value=diffmax)

spread$incand=as.character(spread$gene) %in% as.character(candidates)
ggplot(spread,aes(x=GFPsort,y=RFPsort,color=incand))+geom_point()

geneinfo<-read.csv("geneinfo.csv")
spreaded<-spread
spreaded<-merge(spreaded, geneinfo,by.y="Old.Gene.ID",by.x="gene")
bloodstage<-read.csv("bloodstage.csv")
joint<-merge(joint,bloodstage,by="gene")


spreaded<-merge(spreaded,select(bloodstage,gene,phenotype,Relative.Growth.Rate))


phenolevels=c("Insufficient data","Essential","Slow","Dispensable","Fast","Unselected")  

phenolevelscolor=c("black", "#f90f00", "#0f3791","#007e41","#f83cd9", "darkgray") #standard phenotype colours, colour-blind safe
othercolors=c("#66c2a5","#fc8d62","#8da0cb", "white") #used when another palette is needed
names(phenolevelscolor)=phenolevels

# Define server logic required to draw a histogram
shinyServer(function(input, output) {
   
  
  
  selected<-reactive({
    
    nearPoints(thetable(), input$topplot_click)[1, ]
  })
  
  thetable<-reactive({
    
    newspreaded<-spreaded
    if(input$slow ==F){newspreaded <- filter(newspreaded,phenotype!="Slow" & phenotype!="Essential")}
    if(trimws(input$gene)!="")
    {
    newspreaded$selected=grepl(trimws(input$gene),newspreaded$gene,ignore.case=T) | grepl(trimws(input$gene),newspreaded$current_version_ID,ignore.case=T) | grepl(trimws(input$gene),newspreaded$gene_product,ignore.case=T) | grepl(trimws(input$gene),newspreaded$gene_name,ignore.case=T)
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
    ggplot(ungroup(filter(variance,gene==g,Condition%in%c("control","GFPsort","RFPsort"))),aes(x=Condition,color=Mouse,y=2^(meanlogprop-0.5),group=Mouse))+geom_point(position=position_dodge(width=0.2))+geom_line()+ggtitle(g)+scale_y_log10(breaks=number_ticks(4))+labs(y="Proportion",x="Population")+geom_errorbar(aes(ymin=2^(meanlogprop-2*monotonicmedsd-.5),ymax=2^(meanlogprop+2*monotonicmedsd-.5)),width=0.2,position="dodge",alpha=0.5)+facet_wrap(~Pool, scales="free")+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
  },res=120)
  output$topplot = renderPlot({
    
    p= ggplot(thetable(),aes(x=GFPsort,y=RFPsort,color=phenotype))+geom_point()+scale_color_manual(values = phenolevelscolor)
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
