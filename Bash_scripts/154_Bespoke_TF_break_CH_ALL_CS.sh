#!/bin/bash

MASTER_ROUTE=$1
analysis=$2
TF_REF=$3
TF_ALT=$4


Rscripts_path=$(echo "/home/manuel.tardaguila/Scripts/R/")
module load R/4.1.0

eval "$(conda shell.bash hook)"


output_dir=$(echo "$MASTER_ROUTE""$analysis""/")

Log_files=$(echo "$output_dir""/""Log_files/")

#### bespoke_consultation

eval "$(conda shell.bash hook)"

type=$(echo "$TF_REF""_""$TF_ALT""_""bespoke_consultation")
outfile_bespoke_consultation=$(echo "$Log_files""outfile_8_""$type"".log")
touch $outfile_bespoke_consultation
echo -n "" > $outfile_bespoke_consultation
name_bespoke_consultation=$(echo "$type""_job")
seff_name=$(echo "seff""_""$type")

Rscript_bespoke_consultation=$(echo "$Rscripts_path""449_Inquire_for_TF_changes.R")


Final_table_TF_motif_prediction=$(echo "/group/soranzo/manuel.tardaguila/CH/ALL_variants_in_CS/TF_motif_analysis/Final_table_TF_motif_prediction.tsv")
CHIP_variants_file=$(echo "/group/soranzo/manuel.tardaguila/CH/ALL_variants_in_CS/CHIP_variants.tsv")

# $dependency_string_array

myjobid_bespoke_consultation=$(sbatch --job-name=$name_bespoke_consultation --output=$outfile_bespoke_consultation --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=2 --mem-per-cpu=4096M --parsable --wrap="Rscript $Rscript_bespoke_consultation --Final_table_TF_motif_prediction $Final_table_TF_motif_prediction --TF_REF $TF_REF --TF_ALT $TF_ALT --CHIP_variants_file $CHIP_variants_file --type $type --out $output_dir")
myjobid_seff_bespoke_consultation=$(sbatch --dependency=afterany:$myjobid_bespoke_consultation --open-mode=append --output=$outfile_bespoke_consultation --job-name=$seff_name --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_bespoke_consultation >> $outfile_bespoke_consultation")	



