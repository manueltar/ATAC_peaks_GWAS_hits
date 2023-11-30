
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
  
  #### READ Table_S1 ----
  
  Table_S1<-as.data.frame(fread(file=opt$Table_S1, sep="\t", header=T), stringsAsFactors=F)
  
  Table_S1<-unique(Table_S1[,c(which(colnames(Table_S1) == 'VAR'),which(colnames(Table_S1) == 'rs'))])
  
  cat("Table_S1_0\n")
  cat(str(Table_S1))
  cat("\n")
  cat(str(unique(Table_S1$VAR)))
  cat("\n")
 
  
  #### READ and transform RBC_traits ----
  
  RBC_traits = unlist(strsplit(opt$RBC_traits, split=","))
  
  cat("RBC_traits_\n")
  cat(sprintf(as.character(RBC_traits)))
  cat("\n")
  
  #### READ and transform Mega_traits ----
  
  Mega_traits = unlist(strsplit(opt$Mega_traits, split=","))
  
  cat("Mega_traits_\n")
  cat(sprintf(as.character(Mega_traits)))
  cat("\n")
  
  #### READ and transform GM_traits ----
  
  GM_traits = unlist(strsplit(opt$GM_traits, split=","))
  
  cat("GM_traits_\n")
  cat(sprintf(as.character(GM_traits)))
  cat("\n")
  
  #### READ and transform Lymph_traits ----
  
  Lymph_traits = unlist(strsplit(opt$Lymph_traits, split=","))
  
  cat("Lymph_traits_\n")
  cat(sprintf(as.character(Lymph_traits)))
  cat("\n")
  
  #### Read ATAC_data file ----
  
  ATAC_data<-readRDS(file=opt$ATAC_data)
  
  cat("ATAC_data_0\n")
  cat(str(ATAC_data))
  cat("\n")
  cat(str(unique(ATAC_data$VAR)))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(ATAC_data$variable))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(ATAC_data$variable)))))
  cat("\n")
  
  #### Read multi_ATAC_data file ----
  
  multi_ATAC_data<-readRDS(file=opt$multi_ATAC_data)
  
  cat("multi_ATAC_data_0\n")
  cat(str(multi_ATAC_data))
  cat("\n")
  cat(str(unique(multi_ATAC_data$VAR)))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(multi_ATAC_data$variable))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(multi_ATAC_data$variable)))))
  cat("\n")
  
  #### Read GWAS file ----
  
  MAF_GLOBAL<-as.data.frame(fread(file=opt$MAF_GLOBAL, sep="\t", header=T), stringsAsFactors=F)
  
  cat("MAF_GLOBAL_0\n")
  cat(str(MAF_GLOBAL))
  cat("\n")
  cat(str(unique(MAF_GLOBAL$VAR)))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(MAF_GLOBAL$variable))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(MAF_GLOBAL$variable)))))
  cat("\n")
  
  GWAS_GLOBAL_per_traits<-as.data.frame(fread(file=opt$GWAS_GLOBAL_per_traits, sep="\t", header=T), stringsAsFactors=F)
  
  cat("GWAS_GLOBAL_per_traits_0\n")
  cat(str(GWAS_GLOBAL_per_traits))
  cat("\n")
  cat(str(unique(GWAS_GLOBAL_per_traits$VAR)))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(GWAS_GLOBAL_per_traits$phenotype))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(GWAS_GLOBAL_per_traits$phenotype)))))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(GWAS_GLOBAL_per_traits$variable))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(GWAS_GLOBAL_per_traits$variable)))))
  cat("\n")
  
  GWAS_GLOBAL_per_traits$lineage<-NA
  
  GWAS_GLOBAL_per_traits$lineage[which(GWAS_GLOBAL_per_traits$phenotype%in%RBC_traits)]<-'RBC_lineage'
  GWAS_GLOBAL_per_traits$lineage[which(GWAS_GLOBAL_per_traits$phenotype%in%Mega_traits)]<-'Mega_lineage'
  GWAS_GLOBAL_per_traits$lineage[which(GWAS_GLOBAL_per_traits$phenotype%in%GM_traits)]<-'GM_lineage'
  GWAS_GLOBAL_per_traits$lineage[which(GWAS_GLOBAL_per_traits$phenotype%in%Lymph_traits)]<-'Lymph_lineage'
  
  
  GWAS_GLOBAL_per_traits_CLASSIFIED<-GWAS_GLOBAL_per_traits[!is.na(GWAS_GLOBAL_per_traits$lineage),]
 
  cat("GWAS_GLOBAL_per_traits_CLASSIFIED_0\n")
  cat(str(GWAS_GLOBAL_per_traits_CLASSIFIED))
  cat("\n")
  cat(str(unique(GWAS_GLOBAL_per_traits_CLASSIFIED$VAR)))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(GWAS_GLOBAL_per_traits_CLASSIFIED$phenotype))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(GWAS_GLOBAL_per_traits_CLASSIFIED$phenotype)))))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(GWAS_GLOBAL_per_traits_CLASSIFIED$variable))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(GWAS_GLOBAL_per_traits_CLASSIFIED$variable)))))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(GWAS_GLOBAL_per_traits_CLASSIFIED$lineage))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(GWAS_GLOBAL_per_traits_CLASSIFIED$lineage)))))
  cat("\n")
  
  #### Read Constraint_Z file ----
  
  Constraint_Z<-readRDS(file=opt$Constraint_Z)
  
  cat("Constraint_Z_0\n")
  cat(str(Constraint_Z))
  cat("\n")
  cat(str(unique(Constraint_Z$VAR)))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(Constraint_Z$variable))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(Constraint_Z$variable)))))
  cat("\n")
  
  
  #### Read NCBoost file ----
  
  NCBoost<-readRDS(file=opt$NCBoost)
  
  cat("NCBoost_0\n")
  cat(str(NCBoost))
  cat("\n")
  cat(str(unique(NCBoost$VAR)))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(NCBoost$variable))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(NCBoost$variable)))))
  cat("\n")
  
 
  
  #### Read SpliceAI file ----
  
  SpliceAI<-readRDS(file=opt$SpliceAI)
  
  cat("SpliceAI_0\n")
  cat(str(SpliceAI))
  cat("\n")
  cat(str(unique(SpliceAI$VAR)))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(SpliceAI$variable))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(SpliceAI$variable)))))
  cat("\n")
  
  SpliceAI_subset<-unique(SpliceAI[,-c(which(colnames(SpliceAI) == 'mean_value'),
                                                                which(colnames(SpliceAI) == 'sd_value'))])
  
  cat("SpliceAI_subset_0\n")
  cat(str(SpliceAI_subset))
  cat("\n")
  cat(str(unique(SpliceAI_subset$VAR)))
  cat("\n")
 
  
  #### Read CADD file ----
  
  CADD<-readRDS(file=opt$CADD)
  
  cat("CADD_0\n")
  cat(str(CADD))
  cat("\n")
  cat(str(unique(CADD$VAR)))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(CADD$variable))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(CADD$variable)))))
  cat("\n")
  
 
 
  
  #### Read oe_LOF file ----
  
  oe_LOF<-readRDS(file=opt$oe_LOF)
  
  cat("oe_LOF_0\n")
  cat(str(oe_LOF))
  cat("\n")
  cat(str(unique(oe_LOF$VAR)))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(oe_LOF$variable))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(oe_LOF$variable)))))
  cat("\n")
 
 
  
  #### Read COGS file ----
  
  COGS<-readRDS(file=opt$COGS)
  
  cat("COGS_0\n")
  cat(str(COGS))
  cat("\n")
  cat(str(unique(COGS$VAR)))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(COGS$variable))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(COGS$variable)))))
  cat("\n")
  
  COGS$lineage<-NA
  
  COGS$lineage[which(COGS$phenotype%in%RBC_traits)]<-'RBC_lineage'
  COGS$lineage[which(COGS$phenotype%in%Mega_traits)]<-'Mega_lineage'
  COGS$lineage[which(COGS$phenotype%in%GM_traits)]<-'GM_lineage'
  COGS$lineage[which(COGS$phenotype%in%Lymph_traits)]<-'Lymph_lineage'
  
  
  COGS_CLASSIFIED<-COGS[!is.na(COGS$lineage),]
  
  cat("COGS_CLASSIFIED_0\n")
  cat(str(COGS_CLASSIFIED))
  cat("\n")
  cat(str(unique(COGS_CLASSIFIED$VAR)))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(COGS_CLASSIFIED$phenotype))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(COGS_CLASSIFIED$phenotype)))))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(COGS_CLASSIFIED$variable))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(COGS_CLASSIFIED$variable)))))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(COGS_CLASSIFIED$lineage))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(COGS_CLASSIFIED$lineage)))))
  cat("\n")
 
  #### Read GENE_EXP file ----
  
  GENE_EXP<-readRDS(file=opt$GENE_EXP)
  
 
  cat("GENE_EXP_0\n")
  cat(str(GENE_EXP))
  cat("\n")
  cat(str(unique(GENE_EXP$VAR)))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(GENE_EXP$variable))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(GENE_EXP$variable)))))
  cat("\n")
  
  
  #### Read PCHiC file ----
  
  PCHiC<-readRDS(file=opt$PCHiC)
  
  cat("PCHiC_0\n")
  cat(str(PCHiC))
  cat("\n")
  cat(str(unique(PCHiC$VAR)))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(PCHiC$variable))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(PCHiC$variable)))))
  cat("\n")
  
  #### Read chromstates file ----
  
  chromstates<-readRDS(file=opt$chromstates)
  
  cat("chromstates_0\n")
  cat(str(chromstates))
  cat("\n")
  cat(str(unique(chromstates$VAR)))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(chromstates$variable))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(chromstates$variable)))))
  cat("\n")
  
  ### Read VAR_Prioritization_dB----
  
  
  VAR_Prioritization_dB<-as.data.frame(readRDS(file=opt$VAR_Prioritization_dB) , stringsAsFactors=F)
  
  cat("VAR_Prioritization_dB_0\n")
  cat(str(VAR_Prioritization_dB))
  cat("\n")
  cat(str(unique(VAR_Prioritization_dB$VAR)))
  cat("\n")

  
  VAR_Prioritization_dB_subset<-unique(VAR_Prioritization_dB[,c(which(colnames(VAR_Prioritization_dB) == 'VAR'),
                                                         which(colnames(VAR_Prioritization_dB) == 'Fig1_Annot_Category'))])
  
  cat("VAR_Prioritization_dB_subset_0\n")
  cat(str(VAR_Prioritization_dB_subset))
  cat("\n")
  cat(str(unique(VAR_Prioritization_dB_subset$VAR)))
  cat("\n")
  
  colnames(VAR_Prioritization_dB_subset)[which(colnames(VAR_Prioritization_dB_subset) == 'Fig1_Annot_Category')]<-'V2F_CLASSIF'

  cat("VAR_Prioritization_dB_subset_1\n")
  cat(str(VAR_Prioritization_dB_subset))
  cat("\n")
  cat(str(unique(VAR_Prioritization_dB_subset$VAR)))
  cat("\n")
  
  
  VAR_Prioritization_dB_subset<-merge(VAR_Prioritization_dB_subset,
                               Table_S1,
                               by="VAR")
  
  cat("VAR_Prioritization_dB_subset_2\n")
  cat(str(VAR_Prioritization_dB_subset))
  cat("\n")
  cat(str(unique(VAR_Prioritization_dB_subset$VAR)))
  cat("\n")
 

  ### Var_features ----
  
  VAR_features_df<-rbind(MAF_GLOBAL,Constraint_Z,NCBoost,CADD,ATAC_data,multi_ATAC_data,PCHiC,chromstates)
  
  VAR_features_df$variable<-factor(VAR_features_df$variable,
                               levels=c("MAF",
                                        "CADD_raw",'Constraint_Z',"NCBoost",
                                        "Rank_ATAC_erythroid_lineage","Rank_ATAC_mega_lineage","Rank_ATAC_gran_mono_lineage","Rank_ATAC_lymph_lineage","multi_lineage_ATAC","Rank_PCHiC","Rank_chromstates"),
                               ordered=T)
  
  VAR_features_df<-merge(VAR_features_df,
                         Table_S1,
                         by="VAR")
  
  cat("VAR_features_df_0\n")
  cat(str(VAR_features_df))
  cat("\n")
  cat(str(unique(VAR_features_df$VAR)))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(VAR_features_df$variable))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(VAR_features_df$variable)))))
  cat("\n")
  
  ### VAR_gene ----
  
  VAR_gene_df<-rbind(SpliceAI_subset,oe_LOF,GENE_EXP)
  
  VAR_gene_df$variable<-factor(VAR_gene_df$variable,
                                   levels=c("SpliceAI_DG","SpliceAI_DL","SpliceAI_AG","SpliceAI_AL",
                                            "oe_lof","Rank_GENE_EXP"),
                                   ordered=T)
  
  VAR_gene_df<-merge(VAR_gene_df,
                         Table_S1,
                         by="VAR")
  
  cat("VAR_gene_df_0\n")
  cat(str(VAR_gene_df))
  cat("\n")
  cat(str(unique(VAR_gene_df$VAR)))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(VAR_gene_df$variable))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(VAR_gene_df$variable)))))
  cat("\n")
  
  ### VAR_phenotype ----
  
  VAR_phenotype_df<-GWAS_GLOBAL_per_traits_CLASSIFIED
  
  VAR_phenotype_df$variable<-factor(VAR_phenotype_df$variable,
                               levels=c("PP","Absolute_effect_size","credset_size","ncondind"),
                               ordered=T)
  
  VAR_phenotype_df$lineage<-factor(VAR_phenotype_df$lineage,
                                    levels=c("RBC_lineage","Mega_lineage","GM_lineage","Lymph_lineage"),
                                    ordered=T)
  
  VAR_phenotype_df<-merge(VAR_phenotype_df,
                         Table_S1,
                         by="VAR")
  
  
  cat("VAR_phenotype_df_0\n")
  cat(str(VAR_phenotype_df))
  cat("\n")
  cat(str(unique(VAR_phenotype_df$VAR)))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(VAR_phenotype_df$lineage))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(VAR_phenotype_df$lineage)))))
  cat("\n")
  
  ### VAR_phenotype_gene ----
  
  VAR_phenotype_gene_df<-COGS_CLASSIFIED
  
  VAR_phenotype_gene_df$variable<-factor(VAR_phenotype_gene_df$variable,
                                    levels=c("COGS"),
                                    ordered=T)
  
  VAR_phenotype_gene_df$lineage<-factor(VAR_phenotype_gene_df$lineage,
                                   levels=c("RBC_lineage","Mega_lineage","GM_lineage","Lymph_lineage"),
                                   ordered=T)
  
  VAR_phenotype_gene_df<-merge(VAR_phenotype_gene_df,
                          Table_S1,
                          by="VAR")
  
  cat("VAR_phenotype_gene_df_0\n")
  cat(str(VAR_phenotype_gene_df))
  cat("\n")
  cat(str(unique(VAR_phenotype_gene_df$VAR)))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(VAR_phenotype_gene_df$lineage))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(VAR_phenotype_gene_df$lineage)))))
  cat("\n")

  
  ############## SAVE -------------------------
  
  setwd(out)
  
  saveRDS(VAR_Prioritization_dB_subset,file="VAR_V2F_CLASSIF.rds")
  
  saveRDS(VAR_features_df,file="VAR_features.rds")
  
  saveRDS(VAR_gene_df,file="VAR_gene.rds")
  
  saveRDS(VAR_phenotype_df,file="VAR_phenotype.rds")
  
  saveRDS(VAR_phenotype_gene_df,file="VAR_phenotype_gene.rds")
  
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
    make_option(c("--Table_S1"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--VAR_Prioritization_dB"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--ATAC_data"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--multi_ATAC_data"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--COGS"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--oe_LOF"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--CADD"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--NCBoost"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--Constraint_Z"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--SpliceAI"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--MAF_GLOBAL"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--GWAS_GLOBAL_per_traits"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--GENE_EXP"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--PCHiC"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--chromstates"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--RBC_traits"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--Mega_traits"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--GM_traits"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--Lymph_traits"), type="character", default=NULL, 
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
