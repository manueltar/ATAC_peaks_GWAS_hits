#!/bin/bash

eval "$(conda shell.bash hook)"



Rscripts_path=$(echo "/home/manuel.tardaguila/Scripts/R/")

MASTER_ROUTE=$1
analysis=$2

Master_path_analysis=$(echo "$MASTER_ROUTE""$analysis""/")

#rm -rf $Master_path_analysis
#mkdir -p $Master_path_analysis


 
Log_files=$(echo "$Master_path_analysis""/""Log_files/")

#rm -rf $Log_files
#mkdir -p $Log_files



# #### Prepare_variant_annotation_files ####

module load R/4.1.0

type=$(echo "Prepare_variant_annotation_files""_""$analysis")
outfile_Prepare_variant_annotation_files=$(echo "$Log_files""outfile_""$type"".log")
touch $outfile_Prepare_variant_annotation_files
echo -n "" > $outfile_Prepare_variant_annotation_files
name_Prepare_variant_annotation_files=$(echo "$type""_job")


Rscript_Prepare_variant_annotation_files=$(echo "$Rscripts_path""52_binder_of_Prepared_files_and_add_annotation_layer.R")

Table_S1=$(echo "/group/soranzo/manuel.tardaguila/CHIP_seq_Pimanda/Variant_annotations/Table_S1.tsv")
ATAC_data=$(echo "/group/soranzo/manuel.tardaguila/CHIP_seq_Pimanda/Variant_annotations/Prepared_file_ATAC.rds")
multi_ATAC_data=$(echo "/group/soranzo/manuel.tardaguila/CHIP_seq_Pimanda/Variant_annotations/Prepared_file_multi_lineage_ATAC.rds")
COGS=$(echo "/group/soranzo/manuel.tardaguila/CHIP_seq_Pimanda/Variant_annotations/Prepared_file_COGS.rds")
oe_LOF=$(echo "/group/soranzo/manuel.tardaguila/CHIP_seq_Pimanda/Variant_annotations/Prepared_file_GENE_PLOEUF.rds")
GENE_EXP=$(echo "/group/soranzo/manuel.tardaguila/CHIP_seq_Pimanda/Variant_annotations/Prepared_file_GENE_EXP.rds")
PCHiC=$(echo "/group/soranzo/manuel.tardaguila/CHIP_seq_Pimanda/Variant_annotations/Prepared_file_PCHiC.rds")
chromstates=$(echo "/group/soranzo/manuel.tardaguila/CHIP_seq_Pimanda/Variant_annotations/Prepared_file_chromstates.rds")
CADD=$(echo "/group/soranzo/manuel.tardaguila/CHIP_seq_Pimanda/Variant_annotations/Prepared_file_CADD.rds")
Constraint_Z=$(echo "/group/soranzo/manuel.tardaguila/CHIP_seq_Pimanda/Variant_annotations/Prepared_file_Constraint_Z.rds")
NCBoost=$(echo "/group/soranzo/manuel.tardaguila/CHIP_seq_Pimanda/Variant_annotations/Prepared_file_NCBoost.rds")
SpliceAI=$(echo "/group/soranzo/manuel.tardaguila/CHIP_seq_Pimanda/Variant_annotations/Prepared_file_SpliceAI.rds")
GWAS_GLOBAL_per_traits=$(echo "/group/soranzo/manuel.tardaguila/CHIP_seq_Pimanda/Variant_annotations/GWAS_GLOBAL_per_traits.tsv")
MAF_GLOBAL=$(echo "/group/soranzo/manuel.tardaguila/CHIP_seq_Pimanda/Variant_annotations/MAF_GLOBAL.tsv")
VAR_Prioritization_dB=$(echo "/group/soranzo/manuel.tardaguila/CHIP_seq_Pimanda/Variant_annotations/ER_Labelling_Initial_Selection.rds")


RBC_traits=$(echo 'hct,hgb,hlr,hlr_p,irf,mch,mchc,mcv,mrv,mscv,rbc,rdw_cv,ret,ret_p')
Mega_traits=$(echo 'mpv,pct,pdw,plt')
GM_traits=$(echo 'baso,eo,mono,neut')
Lymph_traits=$(echo 'lymph')

myjobid_Prepare_variant_annotation_files=$(sbatch --job-name=$name_Prepare_variant_annotation_files --output=$outfile_Prepare_variant_annotation_files --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=2 --mem-per-cpu=1024M --parsable --wrap="Rscript $Rscript_Prepare_variant_annotation_files --ATAC_data $ATAC_data --multi_ATAC_data $multi_ATAC_data --COGS $COGS --oe_LOF $oe_LOF --GENE_EXP $GENE_EXP --PCHiC $PCHiC --chromstates $chromstates --CADD $CADD --Table_S1 $Table_S1 --Constraint_Z $Constraint_Z --NCBoost $NCBoost --SpliceAI $SpliceAI --GWAS_GLOBAL_per_traits $GWAS_GLOBAL_per_traits --MAF_GLOBAL $MAF_GLOBAL --VAR_Prioritization_dB $VAR_Prioritization_dB --RBC_traits $RBC_traits --Mega_traits $Mega_traits --GM_traits $GM_traits --Lymph_traits $Lymph_traits --type $type --out $Master_path_analysis")
myjobid_seff_Prepare_variant_annotation_files=$(sbatch --dependency=afterany:$myjobid_Prepare_variant_annotation_files --open-mode=append --output=$outfile_Prepare_variant_annotation_files --job-name=$(echo "seff""_""$name_Prepare_variant_annotation_files") --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_Prepare_variant_annotation_files >> $outfile_Prepare_variant_annotation_files")



# #### Binarize_Pimanda ####

module load R/4.1.0

type=$(echo "Binarize_Pimanda""_""$analysis")
outfile_Binarize_Pimanda=$(echo "$Log_files""outfile_""$type"".log")
touch $outfile_Binarize_Pimanda
echo -n "" > $outfile_Binarize_Pimanda
name_Binarize_Pimanda=$(echo "$type""_job")


Rscript_Binarize_Pimanda=$(echo "$Rscripts_path""53_Binarize_ATAC_Pimanda.R")


ATAC_CHIP_seq=$(echo "/group/soranzo/manuel.tardaguila/CHIP_seq_Pimanda/cell_specific_peaks_HiChIP_loops.tsv")
BIN_1=$(echo '1')
BIN_2=$(echo '2')
BIN_3=$(echo '5')
BIN_4=$(echo '100')


myjobid_Binarize_Pimanda=$(sbatch --job-name=$name_Binarize_Pimanda --output=$outfile_Binarize_Pimanda --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=512M --parsable --wrap="Rscript $Rscript_Binarize_Pimanda --ATAC_CHIP_seq $ATAC_CHIP_seq --BIN_1 $BIN_1 --BIN_2 $BIN_2 --BIN_3 $BIN_3 --BIN_4 $BIN_4 --type $type --out $Master_path_analysis")
myjobid_seff_Binarize_Pimanda=$(sbatch --dependency=afterany:$myjobid_Binarize_Pimanda --open-mode=append --output=$outfile_Binarize_Pimanda --job-name=$(echo "seff""_""$name_Binarize_Pimanda") --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_Binarize_Pimanda >> $outfile_Binarize_Pimanda")





Annotation_file_string=$(echo "VAR_gene.rds,VAR_features.rds,VAR_V2F_CLASSIF.rds,VAR_phenotype.rds,VAR_phenotype_gene.rds")

a=($(echo "$Annotation_file_string" | tr "," '\n'))

declare -a arr

for i  in "${a[@]}"
    do
        Annotation_file_string_sel=${i}
        echo "$Annotation_file_string_sel"
	Annotation_type=$(echo $Annotation_file_string_sel|sed -r 's/.rds//g')
	echo "$Annotation_type"

	#### LiftOver ####

	module load R/4.1.0


	
	type=$(echo "LiftOver""_""$Annotation_type""_""$analysis")
	outfile_LiftOver=$(echo "$Log_files""outfile_""$type"".log")
	touch $outfile_LiftOver
	echo -n "" > $outfile_LiftOver
	name_LiftOver=$(echo "$type""_job")


	Rscript_LiftOver=$(echo "$Rscripts_path""54_liftOver_Annotation.R")


	Annotation_file=$(echo "/group/soranzo/manuel.tardaguila/CHIP_seq_Pimanda/First_analysis/""$Annotation_file_string_sel")



	myjobid_LiftOver=$(sbatch --dependency=afterany:$myjobid_Prepare_variant_annotation_files  --job-name=$name_LiftOver --output=$outfile_LiftOver --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=4 --mem-per-cpu=1024M --parsable --wrap="Rscript $Rscript_LiftOver --Annotation_file $Annotation_file --type $type --out $Master_path_analysis")
	myjobid_seff_LiftOver=$(sbatch --dependency=afterany:$myjobid_LiftOver --open-mode=append --output=$outfile_LiftOver --job-name=$(echo "seff""_""$name_LiftOver") --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_LiftOver >> $outfile_LiftOver")

	#### Intersect ####

	module load R/4.1.0


	type=$(echo "Intersect""_""$Annotation_type""_""$analysis")
	outfile_Intersect=$(echo "$Log_files""outfile_""$type"".log")
	touch $outfile_Intersect
	echo -n "" > $outfile_Intersect
	name_Intersect=$(echo "$type""_job")


	Rscript_Intersect=$(echo "$Rscripts_path""55_Intersect_Annotations.R")


	Annotation_file_LiftOver=$(echo "/group/soranzo/manuel.tardaguila/CHIP_seq_Pimanda/First_analysis/""$Annotation_type""_liftOver.rds")
	ATAC_CHIP_seq_BINNED=$(echo "/group/soranzo/manuel.tardaguila/CHIP_seq_Pimanda/First_analysis/ATAC_CHIP_seq_BINNED.rds")


	myjobid_Intersect=$(sbatch --dependency=afterany:$myjobid_LiftOver:$myjobid_Binarize_Pimanda --job-name=$name_Intersect --output=$outfile_Intersect --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=4 --mem-per-cpu=1024M --parsable --wrap="Rscript $Rscript_Intersect --Annotation_file_LiftOver $Annotation_file_LiftOver --ATAC_CHIP_seq_BINNED $ATAC_CHIP_seq_BINNED --type $type --out $Master_path_analysis")
	myjobid_seff_Intersect=$(sbatch --dependency=afterany:$myjobid_Intersect --open-mode=append --output=$outfile_Intersect --job-name=$(echo "seff""_""$name_Intersect") --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_Intersect >> $outfile_Intersect")

	if [ $Annotation_type = 'VAR_phenotype' ] || [ $Annotation_type = 'VAR_phenotype_gene' ] || [ $Annotation_type = 'VAR_V2F_CLASSIF' ]; then


	    if [ $Annotation_type = 'VAR_V2F_CLASSIF' ]; then

		echo "--------------->Hello_world_1a"

		#### graphical_summary 1b ####

		module load R/4.1.0


		type=$(echo "graphical_summary""_""$Annotation_type""_""$analysis")
		outfile_graphical_summary=$(echo "$Log_files""outfile_""$type"".log")
		touch $outfile_graphical_summary
		echo -n "" > $outfile_graphical_summary
		name_graphical_summary=$(echo "$type""_job")


		Rscript_graphical_summary=$(echo "$Rscripts_path""58_graphical_summary_V2F_classification.R")


		Annotation_file_LiftOver=$(echo "/group/soranzo/manuel.tardaguila/CHIP_seq_Pimanda/First_analysis/""$Annotation_type""_liftOver.rds")
		ATAC_CHIP_seq_BINNED=$(echo "/group/soranzo/manuel.tardaguila/CHIP_seq_Pimanda/First_analysis/ATAC_CHIP_seq_BINNED.rds")
		HITS=$(echo "/group/soranzo/manuel.tardaguila/CHIP_seq_Pimanda/First_analysis/""HITS_""$Annotation_type"".rds")

		#--dependency=afterany:$myjobid_Intersect
		
		myjobid_graphical_summary=$(sbatch --dependency=afterany:$myjobid_Intersect --job-name=$name_graphical_summary --output=$outfile_graphical_summary --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=2 --mem-per-cpu=1024M --parsable --wrap="Rscript $Rscript_graphical_summary --Annotation_file_LiftOver $Annotation_file_LiftOver --ATAC_CHIP_seq_BINNED $ATAC_CHIP_seq_BINNED --HITS $HITS --type $type --out $Master_path_analysis")
		myjobid_seff_graphical_summary=$(sbatch --dependency=afterany:$myjobid_graphical_summary --open-mode=append --output=$outfile_graphical_summary --job-name=$(echo "seff""_""$name_graphical_summary") --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_graphical_summary >> $outfile_graphical_summary")

	    else

		if [ $Annotation_type = 'VAR_phenotype_gene' ]; then

		   echo "--------------->Hello_world_1b1"

		    #### graphical_summary 1b ####

		    module load R/4.1.0


		    type=$(echo "graphical_summary""_""$Annotation_type""_""$analysis")
		    outfile_graphical_summary=$(echo "$Log_files""outfile_""$type"".log")
		    touch $outfile_graphical_summary
		    echo -n "" > $outfile_graphical_summary
		    name_graphical_summary=$(echo "$type""_job")


		    Rscript_graphical_summary=$(echo "$Rscripts_path""57_graphical_summary_var_phenotype_GENES_version.R")


		    Annotation_file_LiftOver=$(echo "/group/soranzo/manuel.tardaguila/CHIP_seq_Pimanda/First_analysis/""$Annotation_type""_liftOver.rds")
		    ATAC_CHIP_seq_BINNED=$(echo "/group/soranzo/manuel.tardaguila/CHIP_seq_Pimanda/First_analysis/ATAC_CHIP_seq_BINNED.rds")
		    HITS=$(echo "/group/soranzo/manuel.tardaguila/CHIP_seq_Pimanda/First_analysis/""HITS_""$Annotation_type"".rds")

		    myjobid_graphical_summary=$(sbatch --dependency=afterany:$myjobid_Intersect --job-name=$name_graphical_summary --output=$outfile_graphical_summary --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=2 --mem-per-cpu=1024M --parsable --wrap="Rscript $Rscript_graphical_summary --Annotation_file_LiftOver $Annotation_file_LiftOver --ATAC_CHIP_seq_BINNED $ATAC_CHIP_seq_BINNED --HITS $HITS --type $type --out $Master_path_analysis")
		    myjobid_seff_graphical_summary=$(sbatch --dependency=afterany:$myjobid_graphical_summary --open-mode=append --output=$outfile_graphical_summary --job-name=$(echo "seff""_""$name_graphical_summary") --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_graphical_summary >> $outfile_graphical_summary")


		else

   		    echo "--------------->Hello_world_1b2"

		    #### graphical_summary 1b ####

		    module load R/4.1.0


		    type=$(echo "graphical_summary""_""$Annotation_type""_""$analysis")
		    outfile_graphical_summary=$(echo "$Log_files""outfile_""$type"".log")
		    touch $outfile_graphical_summary
		    echo -n "" > $outfile_graphical_summary
		    name_graphical_summary=$(echo "$type""_job")


		    Rscript_graphical_summary=$(echo "$Rscripts_path""57_graphical_summary_var_phenotype.R")


		    Annotation_file_LiftOver=$(echo "/group/soranzo/manuel.tardaguila/CHIP_seq_Pimanda/First_analysis/""$Annotation_type""_liftOver.rds")
		    ATAC_CHIP_seq_BINNED=$(echo "/group/soranzo/manuel.tardaguila/CHIP_seq_Pimanda/First_analysis/ATAC_CHIP_seq_BINNED.rds")
		    HITS=$(echo "/group/soranzo/manuel.tardaguila/CHIP_seq_Pimanda/First_analysis/""HITS_""$Annotation_type"".rds")

		    myjobid_graphical_summary=$(sbatch --dependency=afterany:$myjobid_Intersect --job-name=$name_graphical_summary --output=$outfile_graphical_summary --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=2 --mem-per-cpu=1024M --parsable --wrap="Rscript $Rscript_graphical_summary --Annotation_file_LiftOver $Annotation_file_LiftOver --ATAC_CHIP_seq_BINNED $ATAC_CHIP_seq_BINNED --HITS $HITS --type $type --out $Master_path_analysis")
		    myjobid_seff_graphical_summary=$(sbatch --dependency=afterany:$myjobid_graphical_summary --open-mode=append --output=$outfile_graphical_summary --job-name=$(echo "seff""_""$name_graphical_summary") --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_graphical_summary >> $outfile_graphical_summary")

		fi

		
	     fi

	else



	    if [ $Annotation_type = 'VAR_gene' ]; then

 		echo "--------------->Hello_world_1c2"
		
		#### graphical_summary 1c1 ####

		module load R/4.1.0


		type=$(echo "graphical_summary""_""$Annotation_type""_""$analysis")
		outfile_graphical_summary=$(echo "$Log_files""outfile_""$type"".log")
		touch $outfile_graphical_summary
		echo -n "" > $outfile_graphical_summary
		name_graphical_summary=$(echo "$type""_job")


		Rscript_graphical_summary=$(echo "$Rscripts_path""56_graphical_summary_GENES_version.R")


		Annotation_file_LiftOver=$(echo "/group/soranzo/manuel.tardaguila/CHIP_seq_Pimanda/First_analysis/""$Annotation_type""_liftOver.rds")
		ATAC_CHIP_seq_BINNED=$(echo "/group/soranzo/manuel.tardaguila/CHIP_seq_Pimanda/First_analysis/ATAC_CHIP_seq_BINNED.rds")
		HITS=$(echo "/group/soranzo/manuel.tardaguila/CHIP_seq_Pimanda/First_analysis/""HITS_""$Annotation_type"".rds")

		myjobid_graphical_summary=$(sbatch --dependency=afterany:$myjobid_Intersect --job-name=$name_graphical_summary --output=$outfile_graphical_summary --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=4 --mem-per-cpu=1024M --parsable --wrap="Rscript $Rscript_graphical_summary --Annotation_file_LiftOver $Annotation_file_LiftOver --ATAC_CHIP_seq_BINNED $ATAC_CHIP_seq_BINNED --HITS $HITS --type $type --out $Master_path_analysis")
		myjobid_seff_graphical_summary=$(sbatch --dependency=afterany:$myjobid_graphical_summary --open-mode=append --output=$outfile_graphical_summary --job-name=$(echo "seff""_""$name_graphical_summary") --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_graphical_summary >> $outfile_graphical_summary")



	    else

		echo "--------------->Hello_world_1c1"
		
		#### graphical_summary 1c1 ####

		module load R/4.1.0


		type=$(echo "graphical_summary""_""$Annotation_type""_""$analysis")
		outfile_graphical_summary=$(echo "$Log_files""outfile_""$type"".log")
		touch $outfile_graphical_summary
		echo -n "" > $outfile_graphical_summary
		name_graphical_summary=$(echo "$type""_job")


		Rscript_graphical_summary=$(echo "$Rscripts_path""56_graphical_summary.R")


		Annotation_file_LiftOver=$(echo "/group/soranzo/manuel.tardaguila/CHIP_seq_Pimanda/First_analysis/""$Annotation_type""_liftOver.rds")
		ATAC_CHIP_seq_BINNED=$(echo "/group/soranzo/manuel.tardaguila/CHIP_seq_Pimanda/First_analysis/ATAC_CHIP_seq_BINNED.rds")
		HITS=$(echo "/group/soranzo/manuel.tardaguila/CHIP_seq_Pimanda/First_analysis/""HITS_""$Annotation_type"".rds")

		myjobid_graphical_summary=$(sbatch --dependency=afterany:$myjobid_Intersect --job-name=$name_graphical_summary --output=$outfile_graphical_summary --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=4 --mem-per-cpu=1024M --parsable --wrap="Rscript $Rscript_graphical_summary --Annotation_file_LiftOver $Annotation_file_LiftOver --ATAC_CHIP_seq_BINNED $ATAC_CHIP_seq_BINNED --HITS $HITS --type $type --out $Master_path_analysis")
		myjobid_seff_graphical_summary=$(sbatch --dependency=afterany:$myjobid_graphical_summary --open-mode=append --output=$outfile_graphical_summary --job-name=$(echo "seff""_""$name_graphical_summary") --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_graphical_summary >> $outfile_graphical_summary")


	    fi
	    


	    
	fi

    done






