
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
suppressMessages(library("org.Hs.eg.db", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))


opt = NULL

options(warn = 1)

intersect_with_Annotations = function(option_list)
{
  suppressMessages(library("BiocGenerics", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
  suppressMessages(library("S4Vectors", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
  suppressMessages(library("IRanges", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
  suppressMessages(library("GenomeInfoDb", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
  suppressMessages(library("GenomicRanges", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
  suppressMessages(library("Biobase", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
  suppressMessages(library("AnnotationDbi", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
  suppressMessages(library("GO.db", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
  suppressMessages(library("TxDb.Hsapiens.UCSC.hg19.knownGene", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
  suppressMessages(library("rtracklayer", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
  
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
  
  Annotation_file_subset<-unique(Annotation_file[,c(which(colnames(Annotation_file)%in%c('chr','pos_38','ref','alt','VAR_38','rs')))])
  
  Annotation_file_subset$chr<-as.character(Annotation_file_subset$chr)
  
  cat("Annotation_file_subset_0\n")
  cat(str(Annotation_file_subset))
  cat("\n")
  cat(str(unique(Annotation_file_subset$VAR_38)))
  cat("\n")
  
  Annotation_file_subset<-Annotation_file_subset[!is.na(Annotation_file_subset$pos_38),]
  
  cat("Annotation_file_subset_1\n")
  cat(str(Annotation_file_subset))
  cat("\n")
  cat(str(unique(Annotation_file_subset$VAR_38)))
  cat("\n")
  
  # 
  gr_VARS <- GRanges(
    seqnames = as.character(Annotation_file_subset$chr),
    name2=Annotation_file_subset$rs,
    strand="*",
    ranges=IRanges(
      start=Annotation_file_subset$pos_38,
      end=Annotation_file_subset$pos_38,
      names=Annotation_file_subset$VAR_38))
  
  cat("gr_VARS_0\n")
  cat(str(gr_VARS))
  cat("\n")
  
  #### Read annotation files----
  
  ATAC_CHIP_seq<-readRDS(file=opt$ATAC_CHIP_seq_BINNED)
  
  cat("ATAC_CHIP_seq_0\n")
  cat(str(ATAC_CHIP_seq))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(ATAC_CHIP_seq$chr))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(ATAC_CHIP_seq$chr)))))
  cat("\n")
  
  
  gr_ATAC_CHIP_seq <- GRanges(
    seqnames = as.character(ATAC_CHIP_seq$chr),
    name2=ATAC_CHIP_seq$cell_type,
    name3=ATAC_CHIP_seq$Peak_name,
    strand="*",
    ranges=IRanges(
      start=ATAC_CHIP_seq$start,
      end=ATAC_CHIP_seq$end,
      names=paste(ATAC_CHIP_seq$chr,ATAC_CHIP_seq$start,ATAC_CHIP_seq$end,ATAC_CHIP_seq$cell_type, sep='_')))
  
  # cat("gr_ATAC_CHIP_seq\n")
  # cat(str(gr_ATAC_CHIP_seq))
  # cat("\n")
 
  
 
  
  #### find overlaps ----
  
  #### find overlap with SNP ---
  
  m <- findOverlaps(gr_VARS,gr_ATAC_CHIP_seq)
  
  cat("m\n")
  cat(str(m))
  cat("\n")
  
  subjectHits_ATAC<-subjectHits(m)
  
  cat("subjectHits_ATAC\n")
  cat(str(subjectHits_ATAC))
  cat("\n")
  
  queryHits_VARS<-queryHits(m)
  
  cat("queryHits_VARS\n")
  cat(str(queryHits_VARS))
  cat("\n")
  
  VARS_2 <- data.frame(chr_VAR=as.character(seqnames(gr_VARS)),
                       rs=as.character(gr_VARS$name2),
                                start_VAR=as.integer(start(gr_VARS)),
                                end_VAR=as.integer(end(gr_VARS)),
                                VAR_38=names(gr_VARS), stringsAsFactors = F)
  
  cat("VARS_2_0\n")
  cat(str(VARS_2))
  cat("\n")
  
  VARS_2_hits<-VARS_2[queryHits_VARS,]
  
  cat("VARS_2_0\n")
  cat(str(VARS_2_hits))
  cat("\n")
  
  # setwd(out)
  # 
  # write.table(VARS_2_hits,file="VARS_2_hits.tsv",sep="\t", quote=F, row.names = F)
  

  ATAC_CHIP_seq_2 <- data.frame(Peak_name=as.character(gr_ATAC_CHIP_seq$name3),
                                chr=as.character(seqnames(gr_ATAC_CHIP_seq)),
                   start=as.integer(start(gr_ATAC_CHIP_seq)),
                   end=as.integer(end(gr_ATAC_CHIP_seq)),
                   cell_type=as.character(gr_ATAC_CHIP_seq$name2),
                    stringsAsFactors = F)
  
  cat("ATAC_CHIP_seq_2_0\n")
  cat(str(ATAC_CHIP_seq_2))
  cat("\n")
  
  ATAC_CHIP_seq_2_hits<-ATAC_CHIP_seq_2[subjectHits_ATAC,]
  
  cat("ATAC_CHIP_seq_2_0\n")
  cat(str(ATAC_CHIP_seq_2_hits))
  cat("\n")
  
  # setwd(out)
  # 
  # write.table(ATAC_CHIP_seq_2_hits,file="ATAC_CHIP_seq_2_hits.tsv",sep="\t", quote=F, row.names = F)
  
  
  #### The table is prepared to cbind both data.frames ----
  
  Overlap_table<-cbind(ATAC_CHIP_seq_2_hits,VARS_2_hits)
  
  cat("Overlap_table_0\n")
  cat(str(Overlap_table))
  cat("\n")
  
  
  
  #### SAVE ----
  
  filename<-gsub(out,"",as.character(opt$Annotation_file))
  
  cat("filename_0\n")
  str(filename)
  cat("\n")
  
  filename<-gsub("\\.rds","",filename)
  filename<-gsub("_liftOver","",filename)
  
  cat("filename_1\n")
  str(filename)
  cat("\n")
  
  setwd(out)
  
  saveRDS(Overlap_table, file=paste("HITS_",filename,".rds",sep=''))
  write.table(Overlap_table, file=paste("HITS_",filename,".tsv",sep=''), sep="\t",quote=F,row.names = F)
  
  
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
  
  intersect_with_Annotations(opt)

  
}


###########################################################################

system.time( main() )