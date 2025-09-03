suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ape))
suppressPackageStartupMessages(library(phytools))
suppressPackageStartupMessages(library(parallel))

setwd("/home/cc2923/project/PNEUMO_CAPSULE/CAPSULE_SWITCH")


TREE1<-root(ladderize(read.tree("GPS.fasttree.Dec2023.tre")),outgroup="Smitis_CP067992")
TREE1<-drop.tip(root(ladderize(read.tree("GPS.fasttree.Dec2023.tre")),outgroup="Smitis_CP067992"),
               tip=c("Reference","Smitis_CP067992","19084_7#27","24775_1#345"))
#write.tree(TREE1,"GPS.fasttree.Dec2023.R.tre")

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

##for(i in MMF$In_silico_serotype){
##    write.table(i,file=paste0("SER.",gsub("/","",i),".GAIN.LOSS.tsv"),row.names=FALSE,sep="\t",col.names=FALSE)
##}

args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("No arguments provided", call.=FALSE)
} else if (length(args)==1) {
  file.name<-args[1]
  serotype.names<-as_tibble(fread(file.name,header=FALSE,sep="\t")) %>% rename(SER=V1)
}


{
  set.seed(1)
  nsim<-100
  
  TREE<-TREE1
  #TREE$edge.length<-(TREE$edge.length)^(1/2)
  TREE$edge.length[TREE$edge.length==0]<-max(nodeHeights(TREE))*1e-6
  TREE<-multi2di(TREE)

  write.tree(TREE,file="GPS.SPN.R.TREE.tre")
  write.tree(TREE1,file="GPS.SPN.R.TREE1.tre")

  ##sim.obj.df<-list()
  ser.gain.loss.count<-data.frame(serotype=0,acquisition=0,loss=0)
  loop.pos<-1
  
  for(rs in serotype.names$SER ){ 
    
    cat("--->>>",rs,"\n")
    
    MMF1<-MMF %>% mutate(In_silico_serotype=ifelse(In_silico_serotype==rs,1,0))
    mdat<-MMF1[,c("Lane_id","In_silico_serotype")]; 
    .rowNamesDF(mdat, make.names=FALSE)<-MMF1$Lane_id
    
    X<-mdat %>% filter(Lane_id %in% TREE$tip.label)
    ##source.ff<-setNames(eval(parse(text=paste0("X$",rs))),X$lane_id)
    source.ff<-setNames(X$In_silico_serotype,X$Lane_id)
    source.ff<-as.factor(source.ff)
    source.ff<-relevel(source.ff,ref="1")

    fitARD.ALL<-ace(x=source.ff,phy=TREE,model="ARD",type="discrete",CI=TRUE)
    
    {
      cols1<-setNames(c(rgb(1,0,0,alpha=0.85),rgb(0,0,1,alpha=0.05)),levels(source.ff)) 
      cols2<-setNames(c(rgb(1,0,0,alpha=0.05),rgb(0,0,1,alpha=0.05)),levels(source.ff)) 
      
      pdf(paste0("SER_",gsub("/","",rs),"_ACE.pdf")) #,width=9,height=9
      par(fg="transparent")
      
      plot(TREE,lwd=1,
           show.tip.label=FALSE,
           type="fan",edge.color="#000000",edge.width=1,
           align.tip.label=TRUE,open.angle=0)
      #nodelabels(pie=fitARD.ALL$lik.anc,cex=0.35,piecol=cols2,lwd=0,
      #           col=rgb(1,1,1,alpha=0),bg=rgb(1,1,1,alpha=0),
      #           pch=16)
      tiplabels(pie=to.matrix(source.ff[TREE$tip.label],levels(source.ff)),piecol=cols1,
                cex=0.35,col=rgb(1,1,1,alpha=0),bg=rgb(1,1,1,alpha=0),
                pch=16,lwd=0)
      #add.simmap.legend(colors=cols,prompt=FALSE,x=0.9*par()$usr[1],y=0.8*par()$usr[3],fsize=0.8)
      #add.scale.bar(length=20)
      dev.off()
    }
    
    {
      mtree<-make.simmap(TREE,source.ff,model="ARD",nsim=nsim)
      obj<-densityMap(mtree,states=levels(source.ff)[2:1],plot=FALSE)
      
      pdf(paste0("SER_",gsub("/","",rs),"_SIMMAP.pdf"),width=9,height=9);
      cols<-setNames(c(rgb(1,0,0,alpha=0.65),rgb(0,0,1,alpha=0.05)),levels(source.ff)) 
      plot(obj,fsize=c(0.00000006,1),ftype="i",type="fan",colors=cols1,lwd=1.5)
      ###nodelabels(pie=fitARD.ALL$lik.anc,cex=0.35,piecol=cols1,pch = 16)
      tiplabels(col=ifelse(source.ff[TREE$tip.label]==1,cols[1],cols[2]),cex=2.0,
                #border="transparent",
                pch=19) #pch=16 ,bg="#FFFFFF"
      #tiplabels(pie=to.matrix(source.ff[TREE$tip.label],levels(source.ff)),piecol=cols1,
      #          cex=0.35,col=rgb(1,1,1,alpha=0),
      #          bg=rgb(1,1,1,alpha=0),
      #          pch=19,lwd=0)
      dev.off()
      #plot(mtree,cols,type="fan",fsize=0.001,ftype="i",lwd=1.5.5)
      #add.simmap.legend(colors=cols,prompt=FALSE,x=0.9*par()$usr[1],
      #                  y=0.8*par()$usr[3],fsize=0.008)
      #nodelabels(pie=fitARD.ALL$lik.anc,cex=0.35,piecol=cols1,pch=16)
      #tiplabels(pie=to.matrix(source.ff[TREE$tip.label],
      #                        levels(source.ff)),piecol=cols,cex=0.35,pch=16)
      
      num.tr<-describe.simmap(mtree)
      #ser.gain.loss.count<-rbind(ser.gain.loss.count,c(rs,num.tr$Tr[2,1],num.tr$Tr[1,2]))
      if(!is.matrix(num.tr$count)){
        if("0,1" %in% colnames(num.tr$count) ){
          ser.gain.loss.count[loop.pos,]<-c(rs,0,round(mean(num.tr$count[,"1,0"])))
        }else{
          ser.gain.loss.count[loop.pos,]<-(c(rs,round(mean(num.tr$count[,"0,1"])),0))
        }
      }else{
        ser.gain.loss.count[loop.pos,]<-c(rs,round(mean(num.tr$count[,"0,1"])),
                                          round(mean(num.tr$count[,"1,0"])))
      }
    }
    
    ##sim.obj.df[[rs]]<-mtree
    
    loop.pos<-loop.pos+1
  }
}

##saveRDS(sim.obj.df,file="sim.obj.df.rds")
write.table(ser.gain.loss.count,file=paste0("ser.gain.loss.count.",file.name,".tsv"),sep="\t",row.names=FALSE)




