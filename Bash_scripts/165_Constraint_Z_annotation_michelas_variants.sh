#!/bin/bash>

MASTER_ROUTE=$1
output_dir=$MASTER_ROUTE

MASTER_ROUTE=$1
analysis=$2


Rscripts_path=$(echo "/home/manuel.tardaguila/Scripts/R/")
module load R/4.1.0


bashrc_file=$(echo "/home/manuel.tardaguila/.bashrc")

source $bashrc_file
eval "$(conda shell.bash hook)"


output_dir=$(echo "$MASTER_ROUTE""/""$analysis""/")

Log_files=$(echo "$output_dir""Log_files""/")

rm -rf $Log_files
mkdir -p $Log_files
 
#######################################################################binder_of_scores_AbC ###################################


type=$(echo "binder_of_scores_AbC""_""$analysis")
outfile_binder_of_scores_AbC=$(echo "$Log_files""outfile_1_""$type"".log")
touch $outfile_binder_of_scores_AbC
echo -n "" > $outfile_binder_of_scores_AbC
name_binder_of_scores_AbC=$(echo "$type""_job")
seff_name=$(echo "seff""_""$type")

Rscript_binder_of_scores_AbC=$(echo "$Rscripts_path""458_Constraint_Z_score_annotation_Michelas_variants.R")



Table_of_variants=$(echo "/group/soranzo/manuel.tardaguila/CH/ALL_variants_in_CS/config_file_assigned_ref_and_alt.tsv")
Constraint_Z=$(echo "/group/soranzo/manuel.tardaguila/Paper_bits/AbC_Engreitz/constraint_z_genome_1kb.qc.download.txt.gz")


myjobid_binder_of_scores_AbC=$(sbatch --job-name=$name_binder_of_scores_AbC --output=$outfile_binder_of_scores_AbC --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=4 --mem-per-cpu=4096M --parsable --wrap="Rscript $Rscript_binder_of_scores_AbC --Table_of_variants $Table_of_variants --Constraint_Z $Constraint_Z --type $type --out $output_dir")
myjobid_seff_binder_of_scores_AbC=$(sbatch --dependency=afterany:$myjobid_binder_of_scores_AbC --open-mode=append --output=$outfile_binder_of_scores_AbC --job-name=$seff_name --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_binder_of_scores_AbC >> $outfile_binder_of_scores_AbC")
 

