
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
suppressMessages(library("BiocGenerics", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("S4Vectors", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("IRanges", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("GenomeInfoDb", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("GenomicRanges", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("Biobase", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("AnnotationDbi", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("GO.db", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("org.Hs.eg.db", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("TxDb.Hsapiens.UCSC.hg19.knownGene", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("rtracklayer", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("cowplot", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("liftOver", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("rCNV", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))

opt = NULL

options(warn = 1)

data_wrangling_MPRA = function(option_list)
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
  
  
  #### READ Annotation_file ----
  
  Annotation_file<-readRDS(file=opt$Annotation_file)
  
  
  cat("Annotation_file_0\n")
  cat(str(Annotation_file))
  cat("\n")
  cat(str(unique(Annotation_file$VAR)))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(Annotation_file$Prioritisation_bins))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(Annotation_file$Prioritisation_bins)))))
  cat("\n")
  
  
  
  Annotation_file$chr<-gsub("_.+$","",Annotation_file$VAR)
  Annotation_file$pos<-gsub("^[^_]+_","",Annotation_file$VAR)
  Annotation_file$pos<-as.integer(gsub("_.+$","",Annotation_file$pos))
  Annotation_file$ref<-gsub("^[^_]+_[^_]+_","",Annotation_file$VAR)
  Annotation_file$ref<-gsub("_.+$","",Annotation_file$ref)
  Annotation_file$alt<-gsub("^[^_]+_[^_]+_[^_]+_","",Annotation_file$VAR)
  
  cat("Annotation_file_1\n")
  cat(str(Annotation_file))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(Annotation_file$chr))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(Annotation_file$chr)))))
  cat("\n")
  
  Annotation_file$chr<-factor(Annotation_file$chr,
                       levels=c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11",
                                "chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21",
                                "chr22","chr23","chrX","chrY"), ordered=T)
  
  
  Annotation_file<-Annotation_file[order(Annotation_file$chr,Annotation_file$pos),]
  
  cat("Annotation_file_2\n")
  cat(str(Annotation_file))
  cat("\n")
  cat(str(unique(Annotation_file$chr)))
  cat("\n")
  cat(str(unique(Annotation_file$pos)))
  cat("\n")
  cat(str(unique(Annotation_file$ref)))
  cat("\n")
  cat(str(unique(Annotation_file$alt)))
  cat("\n")
  
  indx.int<-c(which(colnames(Annotation_file) == 'VAR'),
              which(colnames(Annotation_file) == 'chr'),which(colnames(Annotation_file) == 'pos'),
              which(colnames(Annotation_file) == 'ref'),which(colnames(Annotation_file) == 'alt'))
  
  Annotation_file_subset<-unique(Annotation_file[,indx.int])
  
  cat("Annotation_file_subset_0\n")
  cat(str(Annotation_file_subset))
  cat("\n")
  
  #### LiftOver 37 -> 38 ----
  
  gr_VARS <- GRanges(
    seqnames = as.character(gsub("chr","",Annotation_file_subset$chr)),
    ranges=IRanges(
      start=as.numeric(Annotation_file_subset$pos),
      end=as.numeric(Annotation_file_subset$pos),
      name=Annotation_file_subset$VAR))
  
  # cat("gr_VARS\n")
  # str(gr_VARS)
  # cat("\n")
  
  VAR_df<-data.frame(chr=as.character(paste('chr',seqnames(gr_VARS), sep='')),
                     pos=start(gr_VARS),
                     ref=Annotation_file_subset$ref,
                     alt=Annotation_file_subset$alt,
                     VAR=Annotation_file_subset$VAR,
                     stringsAsFactors = F)
  
  cat("VAR_df_\n")
  str(VAR_df)
  cat("\n")
  
  #path = system.file(package="liftOver", "extdata", "Hg19Tohg38.over.chain")
  ch = import.chain("/home/manuel.tardaguila/reference_files/hg19ToHg38.over.chain")
  
  seqlevelsStyle(gr_VARS) = "UCSC"  # necessary
  gr_VARS38 = liftOver(gr_VARS, ch)
  gr_VARS38 = unlist(gr_VARS38)
  genome(gr_VARS38) = "hg38"
  
  if(length(gr_VARS38) >0)
  {
    
    chr_38<-as.character(seqnames(gr_VARS38))
    names_38<-as.character(names(gr_VARS38))
    
    ref_VAR38<-gsub("^chr[^_]+_[0-9]+_","",names_38)
    ref_VAR38<-gsub("_.+$","",ref_VAR38)
    
    
    # cat("ref_VAR38\n")
    # cat(sprintf(as.character(ref_VAR38)))
    # cat("\n")
    
    alt_VAR38<-gsub("^chr[^_]+_[0-9]+_[^_]+_","",names_38)
    # alt_VAR38<-gsub("_.+$","",alt_VAR38)
    
    
    # cat("alt_VAR38\n")
    # cat(sprintf(as.character(alt_VAR38)))
    # cat("\n")
    
    
    
    
    VAR_38_df<-data.frame(chr=as.character(seqnames(gr_VARS38)),
                          pos_38=start(gr_VARS38),
                          ref=ref_VAR38,
                          alt=alt_VAR38,
                          VAR=names(gr_VARS38),
                          stringsAsFactors = F)
    
    VAR_38_df$VAR_38<-paste(VAR_38_df$chr,VAR_38_df$pos_38,VAR_38_df$ref,VAR_38_df$alt,sep='_')
    
    
    cat("VAR_38_df_1\n")
    str(VAR_38_df)
    cat("\n")
    
    
    VAR_DEF_df<-unique(merge(VAR_df,
                             VAR_38_df,
                             by=c("chr","ref","alt","VAR"),
                             all=T))
    
    VAR_DEF_df$VAR[is.na(VAR_DEF_df$VAR)]<-"ABSENT"
    
    cat("VAR_DEF_df_2\n")
    str(VAR_DEF_df)
    cat("\n")
    # 
    
    
    
    df<-merge(Annotation_file,
                     VAR_DEF_df,
                     by=c("VAR","chr","pos","ref","alt"),
                     all=T)
    
    cat("df_2\n")
    str(df)
    cat("\n")
    
    #### SAVE ----
    
    filename<-gsub(out,"",as.character(opt$Annotation_file))
    
    cat("filename_0\n")
    str(filename)
    cat("\n")
    
    filename<-gsub("\\.rds","",filename)
    
    cat("filename_1\n")
    str(filename)
    cat("\n")
    
    
    setwd(out)
    
    saveRDS(df, file=paste(filename,"_liftOver",".rds",sep=''))
    
    
  }else{
    
    
    stop("NO_LIFT_OVER\n")
    
  }# length(gr_VARS38) >0
  
  
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
    make_option(c("--Annotation_file"), type="numeric", default=NULL, 
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
  
  data_wrangling_MPRA(opt)
  
  
}


###########################################################################

system.time( main() )