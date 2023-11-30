
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
suppressMessages(library("tidyverse", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1//"))
suppressMessages(library("svglite", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("ggeasy", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("sandwich", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("digest", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("ggforce", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))


opt = NULL

options(warn = 1)

Data_wrangling = function(option_list)
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
  
  cat("OUT_\n")
  cat(sprintf(as.character(out)))
  cat("\n")
  
  #### READ and transform BIN_1 ----
  
  BIN_1 = opt$BIN_1
  
  cat("BIN_1_\n")
  cat(sprintf(as.character(BIN_1)))
  cat("\n")
  
  #### READ and transform BIN_2 ----
  
  BIN_2 = opt$BIN_2
  
  cat("BIN_2_\n")
  cat(sprintf(as.character(BIN_2)))
  cat("\n")
  
  #### READ and transform BIN_3 ----
  
  BIN_3 = opt$BIN_3
  
  cat("BIN_3_\n")
  cat(sprintf(as.character(BIN_3)))
  cat("\n")
  
  #### READ and transform BIN_4 ----
  
  BIN_4 = opt$BIN_4
  
  cat("BIN_4_\n")
  cat(sprintf(as.character(BIN_4)))
  cat("\n")
  
  #### Read ATAC_CHIP_seq----
  
  ATAC_CHIP_seq<-as.data.frame(fread(file=opt$ATAC_CHIP_seq, sep="\t", header=T), stringsAsFactors=F)
  
  colnames(ATAC_CHIP_seq)[6]<-'HSC_HiChIP_loops'
  colnames(ATAC_CHIP_seq)[7]<-'CMP_HiChIP_loops'
  colnames(ATAC_CHIP_seq)[8]<-'GMP_HiChIP_loops'
  colnames(ATAC_CHIP_seq)[9]<-'MEP_HiChIP_loops'
  
  
  cat("ATAC_CHIP_seq_0\n")
  cat(str(ATAC_CHIP_seq))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(ATAC_CHIP_seq$cell_type))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(ATAC_CHIP_seq$cell_type)))))
  cat("\n")
  cat(sprintf(as.character(names(summary(ATAC_CHIP_seq[,6])))))
  cat("\n")
  cat(sprintf(as.character(summary(ATAC_CHIP_seq[,6]))))
  cat("\n")
  cat(sprintf(as.character(names(summary(ATAC_CHIP_seq[,7])))))
  cat("\n")
  cat(sprintf(as.character(summary(ATAC_CHIP_seq[,7]))))
  cat("\n")
  cat(sprintf(as.character(names(summary(ATAC_CHIP_seq[,8])))))
  cat("\n")
  cat(sprintf(as.character(summary(ATAC_CHIP_seq[,8]))))
  cat("\n")
  cat(sprintf(as.character(names(summary(ATAC_CHIP_seq[,9])))))
  cat("\n")
  cat(sprintf(as.character(summary(ATAC_CHIP_seq[,9]))))
  cat("\n")
  
  ### Binarize per cell type ----
  
  ATAC_CHIP_seq.dt<-data.table(ATAC_CHIP_seq, key='cell_type')
  
  
  ATAC_CHIP_seq_breaks_per_cell_type<-as.data.frame(ATAC_CHIP_seq.dt[,.(min_HSC_HiChIP_loops=as.numeric(summary(HSC_HiChIP_loops)[1]),
                                                                     Q1_HSC_HiChIP_loops=as.numeric(summary(HSC_HiChIP_loops)[2]),
                                                                     Med_HSC_HiChIP_loops=as.numeric(summary(HSC_HiChIP_loops)[3]),
                                                                     Q3_HSC_HiChIP_loops=as.numeric(summary(HSC_HiChIP_loops)[5]),
                                                                     max_HSC_HiChIP_loops=as.numeric(summary(HSC_HiChIP_loops)[6]),
                                                                     min_CMP_HiChIP_loops=as.numeric(summary(CMP_HiChIP_loops)[1]),
                                                                     Q1_CMP_HiChIP_loops=as.numeric(summary(CMP_HiChIP_loops)[2]),
                                                                     Med_CMP_HiChIP_loops=as.numeric(summary(CMP_HiChIP_loops)[3]),
                                                                     Q3_CMP_HiChIP_loops=as.numeric(summary(CMP_HiChIP_loops)[5]),
                                                                     max_CMP_HiChIP_loops=as.numeric(summary(CMP_HiChIP_loops)[6]),
                                                                     min_GMP_HiChIP_loops=as.numeric(summary(GMP_HiChIP_loops)[1]),
                                                                     Q1_GMP_HiChIP_loops=as.numeric(summary(GMP_HiChIP_loops)[2]),
                                                                     Med_GMP_HiChIP_loops=as.numeric(summary(GMP_HiChIP_loops)[3]),
                                                                     Q3_GMP_HiChIP_loops=as.numeric(summary(GMP_HiChIP_loops)[5]),
                                                                     max_GMP_HiChIP_loops=as.numeric(summary(GMP_HiChIP_loops)[6]),
                                                                     min_MEP_HiChIP_loops=as.numeric(summary(MEP_HiChIP_loops)[1]),
                                                                     Q1_MEP_HiChIP_loops=as.numeric(summary(MEP_HiChIP_loops)[2]),
                                                                     Med_MEP_HiChIP_loops=as.numeric(summary(MEP_HiChIP_loops)[3]),
                                                                     Q3_MEP_HiChIP_loops=as.numeric(summary(MEP_HiChIP_loops)[5]),
                                                                     max_MEP_HiChIP_loops=as.numeric(summary(MEP_HiChIP_loops)[6])), 
                                                                  by=key(ATAC_CHIP_seq.dt)])
  
  cat("ATAC_CHIP_seq_breaks_per_cell_type_0\n")
  cat(str(ATAC_CHIP_seq_breaks_per_cell_type))
  cat("\n")
  
  ##### toy$Y1 <- cut(toy$value,breaks = c(0,0.2,0.5,0.7,0.8,0.9,0.95,0.999999,Inf),right = FALSE) -----
  
  ATAC_CHIP_seq$BIN_HSC_HiChIP_loops<-cut(ATAC_CHIP_seq$HSC_HiChIP_loops, breaks=c(0,BIN_1,BIN_2,BIN_3,BIN_4), right =FALSE)
  ATAC_CHIP_seq$BIN_CMP_HiChIP_loops<-cut(ATAC_CHIP_seq$CMP_HiChIP_loops, breaks=c(0,BIN_1,BIN_2,BIN_3,BIN_4), right =FALSE)
  ATAC_CHIP_seq$BIN_GMP_HiChIP_loops<-cut(ATAC_CHIP_seq$GMP_HiChIP_loops, breaks=c(0,BIN_1,BIN_2,BIN_3,BIN_4), right =FALSE)
  ATAC_CHIP_seq$BIN_MEP_HiChIP_loops<-cut(ATAC_CHIP_seq$MEP_HiChIP_loops, breaks=c(0,BIN_1,BIN_2,BIN_3,BIN_4), right =FALSE)
  
  
  cat("ATAC_CHIP_seq_1\n")
  cat(str(ATAC_CHIP_seq))
  cat("\n")
  
 
  
  ############## SAVE -------------------------
  
  setwd(out)

  saveRDS(ATAC_CHIP_seq_breaks_per_cell_type,file="ATAC_CHIP_seq_breaks_per_cell_type.rds")
  saveRDS(ATAC_CHIP_seq,file="ATAC_CHIP_seq_BINNED.rds")
  
  
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
    make_option(c("--ATAC_CHIP_seq"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--BIN_1"), type="numeric", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--BIN_2"), type="numeric", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--BIN_3"), type="numeric", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--BIN_4"), type="numeric", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
       make_option(c("--type"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
        make_option(c("--out"), type="character", default=NULL, 
                metavar="filename", 
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
  
  Data_wrangling(opt)
  
}


###########################################################################

system.time( main() )
