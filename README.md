### path: /group/soranzo/manuel.tardaguila/CHIP_seq_Pimanda

$ bash ~/Scripts/Wraper_scripts/15_Pimanda_intersection_and_annotation_v2.sh /group/soranzo/manuel.tardaguila/CHIP_seq_Pimanda/ First_analysis

> summary(df$BIN_HSC_HiChIP_loops)
  [0,1)   [1,2)   [2,5) [5,100)
   9172    1405    1062     217
> summary(df$BIN_CMP_HiChIP_loops)
  [0,1)   [1,2)   [2,5) [5,100)
  11111     481     242      22
> summary(df$BIN_GMP_HiChIP_loops)
  [0,1)   [1,2)   [2,5) [5,100)
   8129    1279    1291    1157
> summary(df$BIN_MEP_HiChIP_loops)
  [0,1)   [1,2)   [2,5) [5,100)
   7695    1453    1696    1012

chr1_KI270713v1_random chr11_KI270830v1_alt chr11_KI270927v1_alt chr15_KI270905v1_alt chr16_KI270853v1_alt chr17_KI270909v1_alt chr22_KI270879v1_alt chr5_GL\
339449v2_alt chr7_KI270803v1_alt chrUn_GL000220v1 chrUn_KI270744v1
4 4 4 8 8 4 4 4 16 4 4

Selections

df<-readRDS(file="Annotated_HITS_RBC_lineage_PP.rds")
df_sel<-droplevels(df[which(df$cell_type == 'HSC' & df$HiChIP_loops == 'HSC_HiChIP_loops' & df$BIN == '>5' & df$value_Z_score >= 1.5),])
rs_sel<-unique(df_sel$rs)
 "rs609018"   "rs28429906"


df<-readRDS(file="Annotated_HITS_GM_lineage_PP.rds")
df_sel<-droplevels(df[which(df$cell_type == 'GMP' & df$HiChIP_loops == 'GMP_HiChIP_loops' & df$BIN == '>5' & df$value_Z_score >= 1),])
rs_sel<-unique(df_sel$rs)

check rs2384952

df<-readRDS(file="Annotated_HITS_Mega_lineage_PP.rds")
df_sel<-droplevels(df[which(df$cell_type == 'MEP' & df$HiChIP_loops == 'MEP_HiChIP_loops' & df$BIN == '>5' & df$value_Z_score >= 2),])
rs_sel<-unique(df_sel$rs)

check "rs75623217" vs "rs2789422"  "rs4714633"

655 different peaks
897 different variants
