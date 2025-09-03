suppressWarnings(suppressMessages(library(tidyverse)))
suppressWarnings(suppressMessages(library(data.table)))
suppressWarnings(suppressMessages(library(ape)))
suppressWarnings(suppressMessages(library(phytools)))
suppressWarnings(suppressMessages(library(combinat)))
suppressWarnings(suppressMessages(library(gtools)))

setwd("/home/cc2923/project/PNEUMO_CAPSULE/CAPSULE_SWITCH")

TREE1<-root(ladderize(read.tree("GPS.fasttree.Dec2023.tre")),outgroup="Smitis_CP067992")
TREE1<-drop.tip(root(ladderize(read.tree("GPS.fasttree.Dec2023.tre")),outgroup="Smitis_CP067992"),
               tip=c("Reference","Smitis_CP067992","19084_7#27","24775_1#345"))


MMF<-as_tibble(fread("GPS.capsule.switch.selected.tsv",header=TRUE,sep="\t")) %>%
  dplyr::filter(Lane_id %in% TREE1$tip.label) %>%
  mutate(In_silico_serotype=ifelse(In_silico_serotype=="39X","10D",In_silico_serotype)) %>%
  mutate(In_silico_serotype=ifelse(In_silico_serotype=="alternative_aliB_NT","alternative_aliB_NT",In_silico_serotype)) %>%
  mutate(In_silico_serotype=ifelse(In_silico_serotype=="ALTERNATIVE_ALIB_NT","alternative_aliB_NT",In_silico_serotype)) %>%
  mutate(In_silico_serotype=ifelse(In_silico_serotype=="Swiss_NT","Swiss_NT",In_silico_serotype)) %>%
  mutate(In_silico_serotype=ifelse(In_silico_serotype=="SWISS_NT","Swiss_NT",In_silico_serotype)) %>%
  mutate(In_silico_serotype=ifelse(In_silico_serotype=="15B","15B/15C",In_silico_serotype)) %>%
  mutate(In_silico_serotype=ifelse(In_silico_serotype=="15C","15B/15C",In_silico_serotype)) %>%
  dplyr::filter(In_silico_serotype!="SEROGROUP 24") %>%
  dplyr::filter(In_silico_serotype!="32") 

args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("No arguments provided", call.=FALSE)
} else if (length(args)==1) {
  file.name<-args[1]
  serotype.names<-as_tibble(fread(file.name,header=FALSE,sep=" ")) %>% rename(SER1=V1,SER2=V2)
}


if(FALSE){
  SER.PAIRS<-data.frame(SER1=NA,SER2=NA)
  pos<-1
  for(i in unique(MMF$GPSC)){
    tmp.DF<-MMF %>% dplyr::filter(In_silico_serotype!="SEROGROUP 24") %>% 
      dplyr::filter(GPSC==i) 
    Serotypes.DF<-unique(tmp.DF$In_silico_serotype)
    
    if(length(Serotypes.DF)>1){
      XX<-permutations(v=Serotypes.DF,r=2,n=length(Serotypes.DF)) 
      colnames(XX)<-c("SER1","SER2")
      XX<-data.frame(XX) %>% 
        mutate(SER1_=ifelse(SER1<SER2,SER1,SER2),
               SER2_=ifelse(SER1<SER2,SER2,SER1)) %>%
        dplyr::select(SER1_,SER2_) %>% unique() %>% 
        rename(SER1=SER1_,SER2=SER2_)
      
      if(i==unique(MMF$GPSC)[1]){
        SER.PAIRS<-XX
      }else{
        SER.PAIRS<-SER.PAIRS %>% bind_rows(XX) 
      }
    }
    pos<-pos+1
  }
  
  SER.PAIRS.FINAL<-as_tibble(SER.PAIRS) %>% distinct() %>% mutate(SAME.GPSC=1) %>%
    bind_rows(
      data.frame(permutations(v=unique(MMF$In_silico_serotype),r=2,
                              n=length(unique(MMF$In_silico_serotype)))) %>%
        rename(SER1=X1,SER2=X2) %>% 
        mutate(SER1_=ifelse(SER1<SER2,SER1,SER2),
               SER2_=ifelse(SER1<SER2,SER2,SER1)) %>%
        dplyr::select(SER1_,SER2_) %>% unique() %>% 
        rename(SER1=SER1_,SER2=SER2_) %>% 
        dplyr::filter(SER1!="SEROGROUP 24") %>%
        dplyr::filter(SER2!="SEROGROUP 24") %>% mutate(SAME.GPSC=0)
    ) %>% group_by(SER1,SER2) %>%
    arrange(desc(SAME.GPSC)) %>% dplyr::filter(row_number()==1) %>%
    dplyr::filter(!is.na(SER1))
  
  for(i in 1:dim(SER.PAIRS.FINAL[SER.PAIRS.FINAL$SAME.GPSC==1,])[1] ){
    cat(i,",")
    SER1<-SER.PAIRS.FINAL[i,1:2]$SER1
    SER2<-SER.PAIRS.FINAL[i,1:2]$SER2
    
    write.table(SER.PAIRS.FINAL[i,1:2],
                file=paste0("SER.PAIRS.",gsub("/","",SER1),".",gsub("/","",SER2),".CAPSULE.SWITCH.tsv"),
                row.names=FALSE,sep="\t",col.names=FALSE)
    
  }
  
}


#######################################################


{
  set.seed(1)
  nsim<-100
  
  TREE<-TREE1
  #TREE$edge.length<-(TREE$edge.length)^(1/2)
  TREE$edge.length[TREE$edge.length==0]<-max(nodeHeights(TREE))*1e-6
  TREE<-multi2di(TREE)
  
  ser.gain.loss.count<-data.frame(SER1=NA,
                                  SER2=NA,
                                  SER3=NA,
                                  SER1.SER2=NA,     
                                  SER2.SER1=NA,     
                                  SER1.SER3=NA,     
                                  SER3.SER1=NA,     
                                  SER2.SER3=NA,     
                                  SER3.SER2=NA)
  
  loop.pos<-1
  sim.obj.df<-list()
  
  for(rs in 1:dim(serotype.names)[1] ){ 
    
    cat("--->>>",rs,"\n")
    
    MMF1<-MMF %>% mutate(In_silico_serotype=ifelse(In_silico_serotype==serotype.names[rs,]$SER1,
                                                   0,
                                                   ifelse(In_silico_serotype==serotype.names[rs,]$SER2,
                                                          1,2))) 
    
    mdat<-MMF1[,c("Lane_id","In_silico_serotype")]; 
    .rowNamesDF(mdat, make.names=FALSE)<-MMF1$Lane_id
    
    X<-mdat %>% filter(Lane_id %in% TREE$tip.label)
    source.ff<-setNames(X$In_silico_serotype,X$Lane_id)
    source.ff<-as.factor(source.ff)
    
    {
      mtree<-make.simmap(TREE,source.ff,model="ARD",nsim=nsim)
      num.tr<-describe.simmap(mtree)
      ser.gain.loss.count[loop.pos,]<-c(serotype.names[rs,]$SER1,
                                        serotype.names[rs,]$SER2,
                                        "Other",
                                        ifelse("0,1" %in% colnames(num.tr$count), round(mean(num.tr$count[,"0,1"])),0),
                                        ifelse("1,0" %in% colnames(num.tr$count), round(mean(num.tr$count[,"1,0"])),0),
                                        ifelse("0,2" %in% colnames(num.tr$count), round(mean(num.tr$count[,"0,2"])),0),
                                        ifelse("2,0" %in% colnames(num.tr$count), round(mean(num.tr$count[,"2,0"])),0),
                                        ifelse("1,2" %in% colnames(num.tr$count), round(mean(num.tr$count[,"1,2"])),0),
                                        ifelse("2,1" %in% colnames(num.tr$count), round(mean(num.tr$count[,"2,1"])),0) )

    }
    
    if(rs==1){
      sim.obj.df[[rs]]<-mtree
    }
    
    loop.pos<-loop.pos+1
  }
  
  write.table(ser.gain.loss.count,file=paste0(file.name,".tsv"),sep="\t",row.names=FALSE)
}

