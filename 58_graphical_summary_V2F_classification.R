
suppressMessages(library("plyr", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("data.table", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("crayon", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("withr", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("ggplot2", lib.loc = "/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("farver", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("labeling", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("optparse", lib.loc = "/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("dplyr", lib.loc = "/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("withr", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("backports", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("broom", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("rstudioapi", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("cli", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("tzdb", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("svglite", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("ggeasy", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("sandwich", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("digest", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("tidyverse", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("RColorBrewer", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("svglite", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("cowplot", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("org.Hs.eg.db", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("ggforce", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("tidyr", lib.loc = "/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))

opt = NULL

options(warn = 1)

graphical_summary_peaks = function(option_list)
{
  opt_in = option_list
  opt <<- option_list
  
  cat("All options:\n")
  printList(opt)
  
  
  #### READ and transform type ----
  
  type = opt$type
  
  cat("TYPE_\n")
  cat(sprintf(as.character(type)))
  cat("\n")
  
  #### READ and transform out ----
  
  out = opt$out
  
  cat("out_\n")
  cat(sprintf(as.character(out)))
  cat("\n")
  
  
  #### Read ATAC_CHIP_seq----
  
  ATAC_CHIP_seq<-readRDS(file=opt$ATAC_CHIP_seq_BINNED)
  
  ATAC_CHIP_seq$cell_type<-factor(ATAC_CHIP_seq$cell_type,
                                  levels=c('HSC','MEP','GMP'),
                                  ordered=T)
  
  cat("ATAC_CHIP_seq_0\n")
  cat(str(ATAC_CHIP_seq))
  cat("\n")
  
  ATAC_CHIP_seq_subset<-unique(ATAC_CHIP_seq[,-c(which(colnames(ATAC_CHIP_seq) == 'HSC_HiChIP_loops'),which(colnames(ATAC_CHIP_seq) == 'CMP_HiChIP_loops'),
                                                           which(colnames(ATAC_CHIP_seq) == 'GMP_HiChIP_loops'),which(colnames(ATAC_CHIP_seq) == 'MEP_HiChIP_loops'),
                                                           which(colnames(ATAC_CHIP_seq) == 'chr_VAR'),which(colnames(ATAC_CHIP_seq) == 'start_VAR'),which(colnames(ATAC_CHIP_seq) == 'end_VAR'))])
  
  cat("ATAC_CHIP_seq_subset_0\n")
  cat(str(ATAC_CHIP_seq_subset))
  cat("\n")
  
  ATAC_CHIP_seq_subset.m<-melt(ATAC_CHIP_seq_subset, id.vars=c("Peak_name","chr","start","end","cell_type"), variable.name='HiChIP_loops', value.name='BIN')
  
  ATAC_CHIP_seq_subset.m$HiChIP_loops<-gsub("BIN_","",ATAC_CHIP_seq_subset.m$HiChIP_loops)
  
  
  
  
  ATAC_CHIP_seq_subset.m$HiChIP_loops<-factor(ATAC_CHIP_seq_subset.m$HiChIP_loops,
                                                   levels=c('HSC_HiChIP_loops','MEP_HiChIP_loops','GMP_HiChIP_loops','CMP_HiChIP_loops'),
                                                   ordered=T)
  
  ATAC_CHIP_seq_subset.m$BIN<-factor(ATAC_CHIP_seq_subset.m$BIN,
                                          levels=c('[0,1)','[1,2)','[2,5)','[5,100)'),
                                          ordered=T)
  
  ATAC_CHIP_seq_subset.m$BIN<-revalue(ATAC_CHIP_seq_subset.m$BIN,
                                           c('[0,1)' = '0',
                                             '[1,2)' = '1',
                                             '[2,5)' = '2-4',
                                             '[5,100)' = '>5'))
  
  ATAC_CHIP_seq_subset.m$PRE_CLASS<-'NOT_FOUND'
  
  cat("ATAC_CHIP_seq_subset.m_0\n")
  cat(str(ATAC_CHIP_seq_subset.m))
  cat("\n")
  cat(sprintf(as.character(names(summary(ATAC_CHIP_seq_subset.m$HiChIP_loops)))))
  cat("\n")
  cat(sprintf(as.character(summary(ATAC_CHIP_seq_subset.m$HiChIP_loops))))
  cat("\n")
  cat(sprintf(as.character(names(summary(ATAC_CHIP_seq_subset.m$BIN)))))
  cat("\n")
  cat(sprintf(as.character(summary(ATAC_CHIP_seq_subset.m$BIN))))
  cat("\n")
  cat(sprintf(as.character(names(summary(ATAC_CHIP_seq_subset.m$cell_type)))))
  cat("\n")
  cat(sprintf(as.character(summary(ATAC_CHIP_seq_subset.m$cell_type))))
  cat("\n")
  
  #### Read HITS----
  
  HITS<-readRDS(file=opt$HITS)
  
  cat("HITS_0\n")
  cat(str(HITS))
  cat("\n")
  
  #### merge ATAC_CHIP_seq_subset + HITS ----
  
  ATAC_CHIP_seq_subset.m_HITS<-merge(ATAC_CHIP_seq_subset.m,
                            HITS,
                            by=c('chr','start','end','cell_type','Peak_name'),
                            all.y=T)
  
  cat("ATAC_CHIP_seq_subset.m_HITS_0\n")
  cat(str(ATAC_CHIP_seq_subset.m_HITS))
  cat("\n")
  
  
 indx.dep<-c(which(colnames(ATAC_CHIP_seq_subset.m_HITS) =='VAR_38'),which(colnames(ATAC_CHIP_seq_subset.m_HITS) =='chr_VAR'),
             which(colnames(ATAC_CHIP_seq_subset.m_HITS) =='start_VAR'),which(colnames(ATAC_CHIP_seq_subset.m_HITS) =='end_VAR'),which(colnames(ATAC_CHIP_seq_subset.m_HITS) =='PRE_CLASS'))
 
 cat("indx.dep_0\n")
 cat(str(indx.dep))
 cat("\n")
 
 ATAC_CHIP_seq_subset.m_HITS<-unique(ATAC_CHIP_seq_subset.m_HITS[,-indx.dep])
 
 ATAC_CHIP_seq_subset.m_HITS$CLASS<-'HIT'
 
 
 cat("ATAC_CHIP_seq_subset.m_HITS_1\n")
 cat(str(ATAC_CHIP_seq_subset.m_HITS))
 cat("\n")
 
 #### merge ATAC_CHIP_seq_subset.m + ATAC_CHIP_seq_subset.m_HITS ----
 
 
 ATAC_CHIP_seq_subset.m<-merge(ATAC_CHIP_seq_subset.m_HITS,
                               ATAC_CHIP_seq_subset.m,
                               by=c('Peak_name','chr','start','end','cell_type','HiChIP_loops','BIN'),
                               all=T)
 
 ATAC_CHIP_seq_subset.m$PRE_CLASS[which(ATAC_CHIP_seq_subset.m$CLASS == 'HIT')]<-'HIT'
 

 cat("ATAC_CHIP_seq_subset.m_1\n")
 cat(str(ATAC_CHIP_seq_subset.m))
 cat("\n")
 cat(sprintf(as.character(names(summary(ATAC_CHIP_seq_subset.m$HiChIP_loops)))))
 cat("\n")
 cat(sprintf(as.character(summary(ATAC_CHIP_seq_subset.m$HiChIP_loops))))
 cat("\n")
 cat(sprintf(as.character(names(summary(ATAC_CHIP_seq_subset.m$BIN)))))
 cat("\n")
 cat(sprintf(as.character(summary(ATAC_CHIP_seq_subset.m$BIN))))
 cat("\n")
 cat(sprintf(as.character(names(summary(ATAC_CHIP_seq_subset.m$cell_type)))))
 cat("\n")
 cat(sprintf(as.character(summary(ATAC_CHIP_seq_subset.m$cell_type))))
 cat("\n")
 cat(sprintf(as.character(names(summary(as.factor(ATAC_CHIP_seq_subset.m$CLASS))))))
 cat("\n")
 cat(sprintf(as.character(summary(as.factor(ATAC_CHIP_seq_subset.m$CLASS)))))
 cat("\n")
 cat(sprintf(as.character(names(summary(as.factor(ATAC_CHIP_seq_subset.m$PRE_CLASS))))))
 cat("\n")
 cat(sprintf(as.character(summary(as.factor(ATAC_CHIP_seq_subset.m$PRE_CLASS)))))
 cat("\n")
 
 
 
 Accepted_chr<-c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11",
   "chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21",
   "chr22","chr23","chrX","chrY")
 
 
 ATAC_CHIP_seq_subset.m$CLASS[which(ATAC_CHIP_seq_subset.m$chr%in%Accepted_chr &
                                      ATAC_CHIP_seq_subset.m$PRE_CLASS =='NOT_FOUND')]<-'NO_HITs'
 
 ATAC_CHIP_seq_subset.m$CLASS[is.na(ATAC_CHIP_seq_subset.m$CLASS)]<-'ATAC_chromosome_not_in_GWAS'
 
 ATAC_CHIP_seq_subset.m$CLASS<-factor(ATAC_CHIP_seq_subset.m$CLASS,
                                      levels=c('ATAC_chromosome_not_in_GWAS','NO_HITs','HIT'),
                                      ordered=T)
 
 
 


 cat("ATAC_CHIP_seq_subset.m_2\n")
 cat(str(ATAC_CHIP_seq_subset.m))
 cat("\n")
 cat(sprintf(as.character(names(summary(ATAC_CHIP_seq_subset.m$HiChIP_loops)))))
 cat("\n")
 cat(sprintf(as.character(summary(ATAC_CHIP_seq_subset.m$HiChIP_loops))))
 cat("\n")
 cat(sprintf(as.character(names(summary(ATAC_CHIP_seq_subset.m$BIN)))))
 cat("\n")
 cat(sprintf(as.character(summary(ATAC_CHIP_seq_subset.m$BIN))))
 cat("\n")
 cat(sprintf(as.character(names(summary(ATAC_CHIP_seq_subset.m$cell_type)))))
 cat("\n")
 cat(sprintf(as.character(summary(ATAC_CHIP_seq_subset.m$cell_type))))
 cat("\n")
 cat(sprintf(as.character(names(summary(ATAC_CHIP_seq_subset.m$CLASS)))))
 cat("\n")
 cat(sprintf(as.character(summary(ATAC_CHIP_seq_subset.m$CLASS))))
 cat("\n")
 cat(sprintf(as.character(names(summary(as.factor(ATAC_CHIP_seq_subset.m$PRE_CLASS))))))
 cat("\n")
 cat(sprintf(as.character(summary(as.factor(ATAC_CHIP_seq_subset.m$PRE_CLASS)))))
 cat("\n")
 
 check<-droplevels(ATAC_CHIP_seq_subset.m[which(ATAC_CHIP_seq_subset.m$CLASS == 'ATAC_chromosome_not_in_GWAS'),])
 
 cat("check_0\n")
 cat(str(check))
 cat("\n")
 cat(sprintf(as.character(names(summary(check$HiChIP_loops)))))
 cat("\n")
 cat(sprintf(as.character(summary(check$HiChIP_loops))))
 cat("\n")
 cat(sprintf(as.character(names(summary(check$BIN)))))
 cat("\n")
 cat(sprintf(as.character(summary(check$BIN))))
 cat("\n")
 cat(sprintf(as.character(names(summary(check$cell_type)))))
 cat("\n")
 cat(sprintf(as.character(summary(check$cell_type))))
 cat("\n")
 cat(sprintf(as.character(names(summary(check$CLASS)))))
 cat("\n")
 cat(sprintf(as.character(summary(check$CLASS))))
 cat("\n")
 cat(sprintf(as.character(names(summary(as.factor(check$PRE_CLASS))))))
 cat("\n")
 cat(sprintf(as.character(summary(as.factor(check$PRE_CLASS)))))
 cat("\n")
 cat(sprintf(as.character(names(summary(as.factor(check$chr))))))
 cat("\n")
 cat(sprintf(as.character(summary(as.factor(check$chr)))))
 cat("\n")
 
 
  ### Print ----
  
  
  
  path_graphs<-paste(out,'ATAC_Peaks_characterization','/',sep='')
  
  if (file.exists(path_graphs)){
    
  }else{
    
    dir.create(file.path(path_graphs))
    
  }#path_graphs
  
  
  DEBUG<-1
  
  
 
  
  ATAC_CHIP_seq_subset.m.dt<-data.table(ATAC_CHIP_seq_subset.m, key=c("cell_type","HiChIP_loops","BIN","CLASS"))
  
  
  
  Freq.table<-as.data.frame(ATAC_CHIP_seq_subset.m.dt[,.(instances =.N), by=key(ATAC_CHIP_seq_subset.m.dt)], stringsAsFactors=F)
  
  if(DEBUG == 1)
  {
    cat("Freq.table_0\n")
    cat(str(Freq.table))
    cat("\n")
  }
  
  ATAC_CHIP_seq_subset.m.dt<-data.table(ATAC_CHIP_seq_subset.m, key=c("cell_type","HiChIP_loops","BIN"))
  
  
  
  Freq.TOTAL<-as.data.frame(ATAC_CHIP_seq_subset.m.dt[,.(TOTAL =.N), by=key(ATAC_CHIP_seq_subset.m.dt)], stringsAsFactors=F)
  
  if(DEBUG == 1)
  {
    cat("Freq.TOTAL_0\n")
    cat(str(Freq.TOTAL))
    cat("\n")
  }
  
  Freq.table<-merge(Freq.table,
                    Freq.TOTAL,
                    by=c("cell_type","HiChIP_loops","BIN"))
  
  Freq.table$Perc<-round(100*(Freq.table$instances/Freq.table$TOTAL),1)
  
  if(DEBUG == 1)
  {
    cat("Freq.table_1\n")
    cat(str(Freq.table))
    cat("\n")
  }
  
  
  
  # quit(status = 1)
  
  #### Graph
  
  breaks.Rank<-(seq(0,100,by=10))
  labels.Rank<-as.character(breaks.Rank)
  
  
  cat(sprintf(as.character(labels.Rank)))
  cat("\n")
  
  vector_colors<-brewer.pal(8,"Set1")
  
  cat("vector_colors_0\n")
  cat(str(vector_colors))
  cat("\n")
  
  
  stacked_barplot<-ggplot(data=Freq.table,
                          aes(x=BIN, y=Perc, fill=CLASS)) +
    geom_bar(stat="identity",colour='black')+
    theme_classic()+
    theme(axis.title.y=element_text(size=18, family="sans"),
          axis.text.y=element_text(angle=0,size=14, color="black", family="sans"))+
    scale_y_continuous(name=paste("Percentage of ATAC Peaks",sep=" "),breaks=breaks.Rank,labels=labels.Rank,
                       limits=c(breaks.Rank[1],breaks.Rank[length(breaks.Rank)]+1))+
    scale_x_discrete(name=NULL, drop=F)+
    scale_fill_manual(values=vector_colors,drop=F)+
    ggeasy::easy_center_title()
  
  
  
  stacked_barplot<-stacked_barplot+
    facet_grid(cell_type ~ HiChIP_loops)+
    theme_cowplot(font_size = 14)+
    theme( strip.background = element_blank(),
           strip.placement = "inside",
           strip.text = element_text(size=14),
           panel.spacing = unit(0.2, "lines"),
           panel.background=element_rect(fill="white"),
           panel.border=element_rect(colour="gray",size=1),
           panel.grid.major = element_blank(),
           panel.grid.minor = element_blank())+
    theme(axis.title.y=element_text(size=18, family="sans"),
          axis.title.x=element_blank(),
          axis.text.y=element_text(size=14, family="sans"),
          axis.text.x=element_text(size=14, family="sans"))+
    theme(legend.title = element_blank(), #change legend title font size
          legend.text = element_text(size=12))+ #change legend text font size
    theme(legend.position = "bottom")+
    guides(fill=guide_legend(nrow=2,byrow=TRUE))
  
  setwd(path_graphs)
  
  svgname<-paste(paste("Stacked_barplot",'Peak_characterization', sep='_'),".svg",sep='')
  
  
  svglite(svgname, width = 10, height = 12)
  print(stacked_barplot)
  dev.off()
  
  ########################################### GRAPH 2 ##########################################
  
  A<-round(summary(Freq.TOTAL$TOTAL[!is.na(Freq.TOTAL$TOTAL)]),0)
  
  step<-abs(A[6]-A[1])/5
  
  if(step == 0)
  {
    
    step<-1
  }
  
  if(DEBUG ==1)
  {
    cat("Summary_TOTAL\n")
    cat(sprintf(as.character(names(A))))
    cat("\n")
    cat(sprintf(as.character(A)))
    cat("\n")
    cat(sprintf(as.character(step)))
    cat("\n")
  }
  
  
  breaks.Rank<-unique(sort(c(0,50,100,seq(from= A[1], to=A[6]+step,by=step))))
  labels.Rank<-as.character(round(breaks.Rank,0))
  
  if(DEBUG ==1)
  {
    cat("labels.Rank\n")
    cat(sprintf(as.character(labels.Rank)))
    cat("\n")
  }
  
  ALL_HITS<-ggplot(data=Freq.TOTAL,
                   aes(y=TOTAL,
                       x=BIN)) +
    geom_bar(stat="identity",colour='black',fill="steelblue")+
    theme_classic()+
    theme(axis.title.y=element_text(size=18, family="sans"),
          axis.text.y=element_text(angle=0,size=14, color="black", family="sans"))+
    scale_y_continuous(name=paste("# ATAC Peaks",sep=" "),breaks=breaks.Rank,labels=labels.Rank,
                       limits=c(breaks.Rank[1],breaks.Rank[length(breaks.Rank)]+1))+
    scale_x_discrete(name=NULL, drop=F)+
    ggeasy::easy_center_title()
  
  ALL_HITS<-ALL_HITS+
    facet_grid(cell_type ~ HiChIP_loops)+
    theme_cowplot(font_size = 14)+
    theme( strip.background = element_blank(),
           strip.placement = "inside",
           strip.text = element_text(size=14),
           panel.spacing = unit(0.2, "lines"),
           panel.background=element_rect(fill="white"),
           panel.border=element_rect(colour="gray",size=1),
           panel.grid.major = element_blank(),
           panel.grid.minor = element_blank())+
    theme(axis.title.y=element_text(size=18, family="sans"),
          axis.title.x=element_blank(),
          axis.text.y=element_text(size=14, family="sans"),
          axis.text.x=element_text(size=14, family="sans"))+
    theme(legend.key.size = unit(1.5, 'cm'), #change legend key size
          legend.key.height = unit(1.5, 'cm'), #change legend key height
          legend.key.width = unit(1, 'cm'), #change legend key width
          legend.title = element_blank(), #change legend title font size
          legend.text = element_text(size=10))+ #change legend text font size
    theme(legend.position = "none")
  
  
  cat("ALL_HITS DONE\n")
  
  setwd(path_graphs)
  
  svgname<-paste(paste("ALL_HITS",'Peak_characterization', sep='_'),".svg",sep='')
  
  
  svglite(svgname, width = 10, height = 12)
  print(ALL_HITS)
  dev.off()
  
  ##### Graph 3 length of the peaks -----
  
  ATAC_CHIP_seq_subset.m$length<-ATAC_CHIP_seq_subset.m$end - ATAC_CHIP_seq_subset.m$start
  
  cat("ATAC_CHIP_seq_subset.m_3\n")
  cat(str(ATAC_CHIP_seq_subset.m))
  cat("\n")
 
  
  A<-round(summary(ATAC_CHIP_seq_subset.m$length[!is.na(ATAC_CHIP_seq_subset.m$length)]),0)
  
  step<-abs(A[6]-A[1])/5
  
  if(step == 0)
  {
    
    step<-1
  }
  
  if(DEBUG ==1)
  {
    cat("Summary_length\n")
    cat(sprintf(as.character(names(A))))
    cat("\n")
    cat(sprintf(as.character(A)))
    cat("\n")
    cat(sprintf(as.character(step)))
    cat("\n")
  }
  
  
  breaks.Rank<-unique(sort(c(0,1000,seq(from= A[1], to=A[6]+step,by=step))))
  labels.Rank<-as.character(round(breaks.Rank,0))
  
  if(DEBUG ==1)
  {
    cat("labels.Rank\n")
    cat(sprintf(as.character(labels.Rank)))
    cat("\n")
  }
  
  
  length_plot<-ggplot(data=ATAC_CHIP_seq_subset.m,
                   aes(y=length,
                       x=BIN,
                       fill=CLASS)) +
    geom_violin(scale = "width", adjust = 1, trim = TRUE)+
    scale_fill_manual(values=vector_colors,drop=F)+
    theme(axis.title.y=element_text(size=18, family="sans"),
          axis.text.y=element_text(angle=0,size=14, color="black", family="sans"))+
    scale_y_continuous(name=paste("Length of ATAC Peaks",sep=" "),breaks=breaks.Rank,labels=labels.Rank,
                       limits=c(breaks.Rank[1],breaks.Rank[length(breaks.Rank)]+1))+
    theme_classic()+
    ggeasy::easy_center_title()
  
  length_plot<-length_plot+
    facet_grid(cell_type ~ HiChIP_loops)+
    theme_cowplot(font_size = 14)+
    scale_x_discrete(name=NULL)+
    theme( strip.background = element_blank(),
           strip.placement = "inside",
           strip.text = element_text(size=14),
           panel.spacing = unit(0.2, "lines"),
           panel.background=element_rect(fill="white"),
           panel.border=element_rect(colour="gray",size=1),
           panel.grid.major = element_blank(),
           panel.grid.minor = element_blank())+
    theme(axis.title.y=element_text(size=18, family="sans"),
          axis.title.x=element_blank(),
          axis.text.y=element_text(size=14, family="sans"),
          axis.text.x=element_text(size=14, family="sans"))+
    theme(legend.title = element_blank(), #change legend title font size
          legend.text = element_text(size=12))+ #change legend text font size
    theme(legend.position = "bottom")+
    guides(fill=guide_legend(nrow=1,byrow=TRUE))
  
  
  cat("length_plot DONE\n")
  
  setwd(path_graphs)
  
  svgname<-paste(paste("length_plot",'Peak_characterization', sep='_'),".svg",sep='')
  
  
  svglite(svgname, width = 10, height = 12)
  print(length_plot)
  dev.off()
  
  
}


graphical_summary = function(option_list)
{
 
  
  opt_in = option_list
  opt <<- option_list
  
  cat("All options:\n")
  printList(opt)
  
  
  #### READ and transform type ----
  
  type = opt$type
  
  cat("TYPE_\n")
  cat(sprintf(as.character(type)))
  cat("\n")
  
  #### READ and transform out ----
  
  out = opt$out
  
  cat("out_\n")
  cat(sprintf(as.character(out)))
  cat("\n")
  
  #### Read annotation files----
  
  Annotation_file<-readRDS(file=opt$Annotation_file_LiftOver)
  
  cat("Annotation_file_0\n")
  cat(str(Annotation_file))
  cat("\n")
  
  Annotation_file_subset<-unique(Annotation_file[,c(which(colnames(Annotation_file) == 'VAR_38'),
                                                    which(colnames(Annotation_file) == 'rs'),
                                                    which(colnames(Annotation_file) == 'V2F_CLASSIF'))])
  
  cat("Annotation_file_subset_0\n")
  cat(str(Annotation_file_subset))
  cat("\n")
  
   #### Read ATAC_CHIP_seq----
  
  ATAC_CHIP_seq<-readRDS(file=opt$ATAC_CHIP_seq_BINNED)
  
  ATAC_CHIP_seq$cell_type<-factor(ATAC_CHIP_seq$cell_type,
                                          levels=c('HSC','MEP','GMP'),
                                          ordered=T)
  
  cat("ATAC_CHIP_seq_0\n")
  cat(str(ATAC_CHIP_seq))
  cat("\n")
 
  
  #### Read HITS----
  
  HITS<-readRDS(file=opt$HITS)
  
  cat("HITS_0\n")
  cat(str(HITS))
  cat("\n")
 
  
  #### merge ATAC_CHIP_seq + HITS ----
  
  ATAC_CHIP_seq_HITS<-merge(ATAC_CHIP_seq,
                       HITS,
                       by=c('chr','start','end','cell_type','Peak_name'),
                       all.y=T)
  
  cat("ATAC_CHIP_seq_HITS_0\n")
  cat(str(ATAC_CHIP_seq_HITS))
  cat("\n")
  
  
  ATAC_CHIP_seq_HITS_subset<-unique(ATAC_CHIP_seq_HITS[,-c(which(colnames(ATAC_CHIP_seq_HITS) == 'HSC_HiChIP_loops'),which(colnames(ATAC_CHIP_seq_HITS) == 'CMP_HiChIP_loops'),
                                                    which(colnames(ATAC_CHIP_seq_HITS) == 'GMP_HiChIP_loops'),which(colnames(ATAC_CHIP_seq_HITS) == 'MEP_HiChIP_loops'),
                                                    which(colnames(ATAC_CHIP_seq_HITS) == 'chr_VAR'),which(colnames(ATAC_CHIP_seq_HITS) == 'start_VAR'),which(colnames(ATAC_CHIP_seq_HITS) == 'end_VAR'))])
  
  cat("ATAC_CHIP_seq_HITS_subset_0\n")
  cat(str(ATAC_CHIP_seq_HITS_subset))
  cat("\n")
  
  ATAC_CHIP_seq_HITS_subset.m<-melt(ATAC_CHIP_seq_HITS_subset, id.vars=c("Peak_name","chr","start","end","cell_type",'VAR_38','rs'), variable.name='HiChIP_loops', value.name='BIN')
  
  ATAC_CHIP_seq_HITS_subset.m$HiChIP_loops<-gsub("BIN_","",ATAC_CHIP_seq_HITS_subset.m$HiChIP_loops)
  
  
  
  
  ATAC_CHIP_seq_HITS_subset.m$HiChIP_loops<-factor(ATAC_CHIP_seq_HITS_subset.m$HiChIP_loops,
                                          levels=c('HSC_HiChIP_loops','MEP_HiChIP_loops','GMP_HiChIP_loops','CMP_HiChIP_loops'),
                                          ordered=T)
  
  ATAC_CHIP_seq_HITS_subset.m$BIN<-factor(ATAC_CHIP_seq_HITS_subset.m$BIN,
                                          levels=c('[0,1)','[1,2)','[2,5)','[5,100)'),
                                          ordered=T)
  
  ATAC_CHIP_seq_HITS_subset.m$BIN<-revalue(ATAC_CHIP_seq_HITS_subset.m$BIN,
                                                    c('[0,1)' = '0',
                                                   '[1,2)' = '1',
                                                   '[2,5)' = '2-4',
                                                   '[5,100)' = '>5'))
  
  cat("ATAC_CHIP_seq_HITS_subset.m_0\n")
  cat(str(ATAC_CHIP_seq_HITS_subset.m))
  cat("\n")
  cat(sprintf(as.character(names(summary(ATAC_CHIP_seq_HITS_subset.m$HiChIP_loops)))))
  cat("\n")
  cat(sprintf(as.character(summary(ATAC_CHIP_seq_HITS_subset.m$HiChIP_loops))))
  cat("\n")
  cat(sprintf(as.character(names(summary(ATAC_CHIP_seq_HITS_subset.m$BIN)))))
  cat("\n")
  cat(sprintf(as.character(summary(ATAC_CHIP_seq_HITS_subset.m$BIN))))
  cat("\n")
  cat(sprintf(as.character(names(summary(ATAC_CHIP_seq_HITS_subset.m$cell_type)))))
  cat("\n")
  cat(sprintf(as.character(summary(ATAC_CHIP_seq_HITS_subset.m$cell_type))))
  cat("\n")
  
  ### Print ----
  
  
  
  Type_of_variables<-gsub(out,"",as.character(opt$Annotation_file))
  
  cat("Type_of_variables_0\n")
  str(Type_of_variables)
  cat("\n")
  
  Type_of_variables<-gsub("\\.rds","",Type_of_variables)
  Type_of_variables<-gsub("_liftOver","",Type_of_variables)
  
  cat("Type_of_variables_1\n")
  str(Type_of_variables)
  cat("\n")
  
  
  path_graphs<-paste(out,Type_of_variables,'/',sep='')
  
  if (file.exists(path_graphs)){
    
  }else{
    
    dir.create(file.path(path_graphs))
    
  }#path_graphs
  
  
  DEBUG<-1
  
  
  Annotation_file_subset_HITS<-Annotation_file_subset[which(Annotation_file_subset$VAR_38%in%ATAC_CHIP_seq_HITS_subset.m$VAR_38),]
  
  if(DEBUG ==1)
  {
    cat("Annotation_file_subset_HITS_0\n")
    cat(str(Annotation_file_subset_HITS))
    cat("\n")
  }
  Annotation_file_subset_HITS<-merge(Annotation_file_subset_HITS,
                                         ATAC_CHIP_seq_HITS_subset.m,
                                         by=c("VAR_38","rs"))
  
  if(DEBUG ==1)
  {
    cat("Annotation_file_subset_HITS_1\n")
    cat(str(Annotation_file_subset_HITS))
    cat("\n")
  }
  
  setwd(path_graphs)
  
  filename<-paste(paste("Annotated_HITS",features_sel, sep='_'),".rds",sep='')
  
  saveRDS(Annotation_file_subset_sel_HITS, file=filename)
  
  filename<-paste(paste("Annotated_HITS",features_sel, sep='_'),".tsv",sep='')
  
  write.table(Annotation_file_subset_sel_HITS, file=filename, sep="\t",quote=F, row.names = F)
  
  
  Annotation_file_subset_HITS.dt<-data.table(Annotation_file_subset_HITS, key=c("cell_type","HiChIP_loops","BIN","V2F_CLASSIF"))
  
  
  
  Freq.table<-as.data.frame(Annotation_file_subset_HITS.dt[,.(instances =.N), by=key(Annotation_file_subset_HITS.dt)], stringsAsFactors=F)
  
  if(DEBUG == 1)
  {
    cat("Freq.table_0\n")
    cat(str(Freq.table))
    cat("\n")
  }
  
  Annotation_file_subset_HITS.dt<-data.table(Annotation_file_subset_HITS, key=c("cell_type","HiChIP_loops","BIN"))
  
  
  
  Freq.TOTAL<-as.data.frame(Annotation_file_subset_HITS.dt[,.(TOTAL =.N), by=key(Annotation_file_subset_HITS.dt)], stringsAsFactors=F)
  
  if(DEBUG == 1)
  {
    cat("Freq.TOTAL_0\n")
    cat(str(Freq.TOTAL))
    cat("\n")
  }
  
  Freq.table<-merge(Freq.table,
                    Freq.TOTAL,
                    by=c("cell_type","HiChIP_loops","BIN"))
  
  Freq.table$Perc<-round(100*(Freq.table$instances/Freq.table$TOTAL),1)
  
  if(DEBUG == 1)
  {
    cat("Freq.table_1\n")
    cat(str(Freq.table))
    cat("\n")
  }
  
  
  
  # quit(status = 1)
  
  #### Graph
  
  breaks.Rank<-(seq(0,100,by=10))
  labels.Rank<-as.character(breaks.Rank)
  
  
  cat(sprintf(as.character(labels.Rank)))
  cat("\n")
  
  vector_colors<-brewer.pal(8,"Set1")
  
  cat("vector_colors_0\n")
  cat(str(vector_colors))
  cat("\n")
  
  
  stacked_barplot<-ggplot(data=Freq.table,
                          aes(x=BIN, y=Perc, fill=V2F_CLASSIF)) +
    geom_bar(stat="identity",colour='black')+
    theme_classic()+
    theme(axis.title.y=element_text(size=18, family="sans"),
          axis.text.y=element_text(angle=0,size=14, color="black", family="sans"))+
    scale_y_continuous(name=paste("Percentage of variants",sep=" "),breaks=breaks.Rank,labels=labels.Rank,
                       limits=c(breaks.Rank[1],breaks.Rank[length(breaks.Rank)]+1))+
    scale_x_discrete(name=NULL, drop=F)+
    scale_fill_manual(values=vector_colors,drop=F)+
    ggeasy::easy_center_title()
  
  
  
  stacked_barplot<-stacked_barplot+
    facet_grid(cell_type ~ HiChIP_loops)+
    theme_cowplot(font_size = 14)+
    theme( strip.background = element_blank(),
           strip.placement = "inside",
           strip.text = element_text(size=14),
           panel.spacing = unit(0.2, "lines"),
           panel.background=element_rect(fill="white"),
           panel.border=element_rect(colour="gray",size=1),
           panel.grid.major = element_blank(),
           panel.grid.minor = element_blank())+
    theme(axis.title.y=element_text(size=18, family="sans"),
          axis.title.x=element_blank(),
          axis.text.y=element_text(size=14, family="sans"),
          axis.text.x=element_text(size=14, family="sans"))+
    theme(legend.title = element_blank(), #change legend title font size
          legend.text = element_text(size=12))+ #change legend text font size
    theme(legend.position = "bottom")+
    guides(fill=guide_legend(nrow=2,byrow=TRUE))
  
  setwd(path_graphs)
  
  svgname<-paste(paste("Stacked_barplot",Type_of_variables, sep='_'),".svg",sep='')
  
  
  svglite(svgname, width = 10, height = 12)
  print(stacked_barplot)
  dev.off()
  
  ########################################### GRAPH 2 ##########################################
  
  A<-round(summary(Freq.TOTAL$TOTAL[!is.na(Freq.TOTAL$TOTAL)]),0)
  
  step<-abs(A[6]-A[1])/5
  
  if(step == 0)
  {
    
    step<-1
  }
  
  if(DEBUG ==1)
  {
    cat("Summary_TOTAL\n")
    cat(sprintf(as.character(names(A))))
    cat("\n")
    cat(sprintf(as.character(A)))
    cat("\n")
    cat(sprintf(as.character(step)))
    cat("\n")
  }
  
  
  breaks.Rank<-unique(sort(c(0,seq(from= A[1], to=A[6]+step,by=step))))
  labels.Rank<-as.character(round(breaks.Rank,0))
  
  if(DEBUG ==1)
  {
    cat("labels.Rank\n")
    cat(sprintf(as.character(labels.Rank)))
    cat("\n")
  }
  
  ALL_HITS<-ggplot(data=Freq.TOTAL,
                   aes(y=TOTAL,
                       x=BIN)) +
    geom_bar(stat="identity",colour='black',fill="steelblue")+
    theme_classic()+
    theme(axis.title.y=element_text(size=18, family="sans"),
          axis.text.y=element_text(angle=0,size=14, color="black", family="sans"))+
    scale_y_continuous(name=paste("# variants",sep=" "),breaks=breaks.Rank,labels=labels.Rank,
                       limits=c(breaks.Rank[1],breaks.Rank[length(breaks.Rank)]+1))+
    scale_x_discrete(name=NULL, drop=F)+
    ggeasy::easy_center_title()
  
  ALL_HITS<-ALL_HITS+
    facet_grid(cell_type ~ HiChIP_loops)+
    theme_cowplot(font_size = 14)+
    theme( strip.background = element_blank(),
           strip.placement = "inside",
           strip.text = element_text(size=14),
           panel.spacing = unit(0.2, "lines"),
           panel.background=element_rect(fill="white"),
           panel.border=element_rect(colour="gray",size=1),
           panel.grid.major = element_blank(),
           panel.grid.minor = element_blank())+
    theme(axis.title.y=element_text(size=18, family="sans"),
          axis.title.x=element_blank(),
          axis.text.y=element_text(size=14, family="sans"),
          axis.text.x=element_text(size=14, family="sans"))+
    theme(legend.key.size = unit(1.5, 'cm'), #change legend key size
          legend.key.height = unit(1.5, 'cm'), #change legend key height
          legend.key.width = unit(1, 'cm'), #change legend key width
          legend.title = element_blank(), #change legend title font size
          legend.text = element_text(size=10))+ #change legend text font size
    theme(legend.position = "none")
  
  
  cat("ALL_HITS DONE\n")
  
  setwd(path_graphs)
  
  svgname<-paste(paste("ALL_HITS",Type_of_variables, sep='_'),".svg",sep='')
  
  
  svglite(svgname, width = 10, height = 12)
  print(ALL_HITS)
  dev.off()
  
  
 
 
}



printList = function(l, prefix = "    ") {
  list.df = data.frame(val_name = names(l), value = as.character(l))
  list_strs = apply(list.df, MARGIN = 1, FUN = function(x) { paste(x, collapse = " = ")})
  cat(paste(paste(paste0(prefix, list_strs), collapse = "\n"), "\n"))
}


#### main script ----

main = function() {
  cmd_line = commandArgs()
  cat("Command line:\n")
  cat(paste(gsub("--file=", "", cmd_line[4], fixed=T),
            paste(cmd_line[6:length(cmd_line)], collapse = " "),
            "\n\n"))
  option_list <- list(
    make_option(c("--ATAC_CHIP_seq_BINNED"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--Annotation_file_LiftOver"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--HITS"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--type"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--out"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required.")
  )
  parser = OptionParser(usage = "140__Rscript_v106.R
                        --subset type
                        --TranscriptEXP FILE.txt
                        --cadd FILE.txt
                        --ncboost FILE.txt
                        --type type
                        --out filename",
                        option_list = option_list)
  opt <<- parse_args(parser)
  
  graphical_summary_peaks(opt)
  graphical_summary(opt)

  
}


###########################################################################

system.time( main() )