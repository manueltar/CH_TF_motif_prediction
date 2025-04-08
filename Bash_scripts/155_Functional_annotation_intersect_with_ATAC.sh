#!/bin/bash

MASTER_ROUTE=$1
analysis=$2

Rscripts_path=$(echo "/home/manuel.tardaguila/Scripts/R/")
module load R/4.1.0

eval "$(conda shell.bash hook)"


output_dir=$(echo "$MASTER_ROUTE""$analysis""/")

Log_files=$(echo "$output_dir""/""Log_files/")

rm -rf $Log_files
mkdir -p $Log_files

#### Intersect_ATAC ####


type=$(echo "Intersect_ATAC")
outfile_Intersect_ATAC=$(echo "$Log_files""outfile_1_""$type"".log")
touch $outfile_Intersect_ATAC
echo -n "" > $outfile_Intersect_ATAC
name_Intersect_ATAC=$(echo "$type""_job")


Rscript_Intersect_ATAC=$(echo "$Rscripts_path""450_Intersect_ATAC_funtional_annotation.R")


ATAC_peaks=$(echo "/group/soranzo/manuel.tardaguila/Rebuild_sit_plots/Reference_files/29August2017_EJCsamples_allReads_500bp.bed")
ATAC_counts=$(echo "/group/soranzo/manuel.tardaguila/Rebuild_sit_plots/Reference_files/29August2017_EJCsamples_allReads_500bp.counts.txt")
Table_of_variants=$(echo "/group/soranzo/manuel.tardaguila/CH/ALL_variants_in_CS/config_file_assigned_ref_and_alt.tsv")

myjobid_Intersect_ATAC=$(sbatch  --job-name=$name_Intersect_ATAC --output=$outfile_Intersect_ATAC --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=3 --mem-per-cpu=1024M --parsable --wrap="Rscript $Rscript_Intersect_ATAC --ATAC_peaks $ATAC_peaks --ATAC_counts $ATAC_counts --Table_of_variants $Table_of_variants  --type $type --out $output_dir")
myjobid_seff_Intersect_ATAC=$(sbatch --dependency=afterany:$myjobid_Intersect_ATAC --open-mode=append --output=$outfile_Intersect_ATAC --job-name=$(echo "seff""_""$name_Intersect_ATAC") --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_Intersect_ATAC >> $outfile_Intersect_ATAC")

