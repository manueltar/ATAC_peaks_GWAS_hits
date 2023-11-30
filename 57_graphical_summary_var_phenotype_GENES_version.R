
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
  
  Annotation_file_subset<-unique(Annotation_file[,c(which(colnames(Annotation_file) == 'VAR_38'),which(colnames(Annotation_file) == 'rs'),
                                                    which(colnames(Annotation_file) == 'ensembl_gene_id'),which(colnames(Annotation_file) == 'HGNC'),
                                                    which(colnames(Annotation_file) == 'phenotype'),which(colnames(Annotation_file) == 'lineage'),
                                                    which(colnames(Annotation_file) == 'variable'),
                                             which(colnames(Annotation_file) == 'value'),which(colnames(Annotation_file) == 'value_Z_score'))])
  
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
  
  ### LOOP print ----
  
  array_features<-levels(Annotation_file_subset$variable)
  
  cat("array_features_0\n")
  cat(str(array_features))
  cat("\n")
  
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
  
  
  DEBUG<-0
  
  for(i in 1:length(array_features))
  {
    features_sel<-array_features[i]
    
    cat("------------------------------------------------------------------------------------------------------------------------->\t")
    cat(sprintf(as.character(features_sel)))
    cat("\n")
    
    Annotation_file_subset_sel<-droplevels(Annotation_file_subset[which(Annotation_file_subset$variable%in%features_sel),])
    
    if(DEBUG ==1)
    {
      cat("Annotation_file_subset_sel_0\n")
      cat(str(Annotation_file_subset_sel))
      cat("\n")
    }
    
    Annotation_file_subset_sel_HITS<-Annotation_file_subset_sel[which(Annotation_file_subset_sel$VAR_38%in%ATAC_CHIP_seq_HITS_subset.m$VAR_38),]
    
    if(DEBUG ==1)
    {
      cat("Annotation_file_subset_sel_HITS_0\n")
      cat(str(Annotation_file_subset_sel_HITS))
      cat("\n")
    }
    
   
    
    
    array_lineages<-levels(Annotation_file_subset_sel_HITS$lineage)
    
    if(DEBUG ==1)
    {
      cat("array_lineages_0\n")
      cat(str(array_lineages))
      cat("\n")
    }
    
    for(k in 1:length(array_lineages))
    {
      lineage_sel<-array_lineages[k]
      
      cat("---------------------------------->\t")
      cat(sprintf(as.character(lineage_sel)))
      cat("\n")
      
      Annotation_file_subset_sel_HITS_lineage<-droplevels(Annotation_file_subset_sel_HITS[which(Annotation_file_subset_sel_HITS$lineage == lineage_sel),])
      
      if(DEBUG ==1)
      {
        cat("Annotation_file_subset_sel_HITS_lineage_0\n")
        cat(str(Annotation_file_subset_sel_HITS_lineage))
        cat("\n")
      }
      
   
      
      Annotation_file_subset_sel_HITS_lineage<-merge(Annotation_file_subset_sel_HITS_lineage,
                                                     ATAC_CHIP_seq_HITS_subset.m,
                                                     by=c("VAR_38","rs"))
      
      if(DEBUG ==1)
      {
        cat("Annotation_file_subset_sel_HITS_lineage_1\n")
        cat(str(Annotation_file_subset_sel_HITS_lineage))
        cat("\n")
      }
      
      setwd(path_graphs)
      
      filename<-paste(paste("Annotated_HITS",lineage_sel,features_sel, sep='_'),".rds",sep='')
      
      saveRDS(Annotation_file_subset_sel_HITS_lineage, file=filename)
      
      filename<-paste(paste("Annotated_HITS",lineage_sel,features_sel, sep='_'),"tsv",sep='')
      
      write.table(Annotation_file_subset_sel_HITS_lineage, file=filename, sep="\t",quote=F, row.names = F)
      
      
      ########################################### GRAPH ##########################################
      
      A<-round(summary(Annotation_file_subset_sel_HITS_lineage$value_Z_score[!is.na(Annotation_file_subset_sel_HITS_lineage$value_Z_score)]),2)
      
      step<-abs(A[6]-A[1])/5
      
      if(step == 0)
      {
        
        step<-1
      }
      
      if(DEBUG ==1)
      {
        cat("Summary_value_Z_score\n")
        cat(sprintf(as.character(names(A))))
        cat("\n")
        cat(sprintf(as.character(A)))
        cat("\n")
        cat(sprintf(as.character(step)))
        cat("\n")
      }
      
      
      breaks.Rank<-unique(sort(c(0,seq(from= A[1], to=A[6]+step,by=step))))
      labels.Rank<-as.character(round(breaks.Rank,2))
      
      if(DEBUG ==1)
      {
        cat("labels.Rank\n")
        cat(sprintf(as.character(labels.Rank)))
        cat("\n")
      }
      
      #### sina _plot ----
      
      graph_variable<-ggplot()+
        geom_sina(data=Annotation_file_subset_sel_HITS_lineage,
                  aes(x=BIN, 
                      y=value_Z_score), size=3)+
        geom_boxplot(data=Annotation_file_subset_sel_HITS_lineage,
                     aes(x=BIN, 
                         y=value_Z_score),
                     notch = TRUE,
                     notchwidth = 0.5)+
        theme_bw()+
        scale_x_discrete(name=NULL, drop=F)+
        scale_y_continuous(name=paste(features_sel,"Z score", sep='   '),
                           breaks=breaks.Rank,labels=labels.Rank, 
                           limits=c(breaks.Rank[1]-0.01,breaks.Rank[length(breaks.Rank)]+0.01))
      
      if(DEBUG == 1)
      {
        cat("Graph_Part_I:-----------------------------------------------------------------------------------------------\n")
        
      }
      
      
      graph_variable<-graph_variable+
        facet_grid(Annotation_file_subset_sel_HITS_lineage$cell_type ~ Annotation_file_subset_sel_HITS_lineage$HiChIP_loops, scales='free_x', space='free_x') +
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
        theme(legend.position = "hidden")+
        ggeasy::easy_center_title()
      
      if(DEBUG == 1)
      {
        cat("Graph_Part_II:-----------------------------------------------------------------------------------------------\n")
        
      }
      
      setwd(path_graphs)
      
      svgname<-paste(paste(lineage_sel,features_sel,"Z_score", sep='_'),".svg",sep='')
      # svgname<-'test.svg'
      makesvg = TRUE
      
      if (makesvg == TRUE)
      {
        ggsave(svgname, plot= graph_variable,
               device="svg",
               height=10, width=12)
      }
      
      if(DEBUG == 1)
      {
        cat("Graph_Part_END:-----------------------------------------------------------------------------------------------\n")
       
      }
      
    }#k in 1:length(array_lineages)
  }# i in 1:length(array_features)
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
  
  graphical_summary(opt)

  
}


###########################################################################

system.time( main() )