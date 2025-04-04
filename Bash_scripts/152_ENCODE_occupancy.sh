#!/bin/bash
 
eval "$(conda shell.bash hook)"
 

Rscripts_path=$(echo "/home/manuel.tardaguila/Scripts/R/")

MASTER_ROUTE=$1
analysis=$2


output_dir=$(echo "$MASTER_ROUTE""$analysis""/")

Master_path_analysis=$(echo "$output_dir")

 
Log_files=$(echo "$Master_path_analysis""/""Log_files/")

table_array=$(echo 'FULL,NoMosaic,tet2,DNMT')
a=($(echo "$table_array" | tr "," '\n'))


declare -a global=()

global_length=${#a[@]}

for (( iteration=0; iteration<${global_length}; iteration=iteration+1 ));
do

    table_array_sel=${a[$iteration]}

    echo "------------------------------------------------------------------------------------------------------------->$table_array_sel"


    ####################### Peaks are in hg19 liftOver
    ####################### Peaks are in hg19 liftOver
    ####################### Peaks are in hg19 liftOver

    type=$(echo "$table_array_sel""_""ENCODE_Match_TF_occupancy")
    outfile_ENCODE_Match_TF_occupancy=$(echo "$Log_files""outfile_3_""$type"".log")
    touch $outfile_ENCODE_Match_TF_occupancy
    echo -n "" > $outfile_ENCODE_Match_TF_occupancy
    name_ENCODE_Match_TF_occupancy=$(echo "$type""_job")
    seff_name=$(echo "seff""_""$type")


    Rscript_ENCODE_Match_TF_occupancy=$(echo "$Rscripts_path""447_ENCODE_TF_occupancy.R")
    
    
    Final_table_TF_motif_prediction=$(echo "/group/soranzo/manuel.tardaguila/CH/""Final_table_TF_motif_prediction_with_chipseq_support_""$table_array_sel"".tsv")
    ENCODE_file=$(echo "/group/soranzo/manuel.tardaguila/CH/reference_files/ENCODE/CHIPseq/ENCODE_CHIP_data.tsv.gz")
    out2=$(echo "/group/soranzo/manuel.tardaguila/CH/")

    myjobid_ENCODE_Match_TF_occupancy=$(sbatch --output=$outfile_ENCODE_Match_TF_occupancy --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=4 --mem-per-cpu=1024 --parsable --job-name $name_ENCODE_Match_TF_occupancy --wrap="Rscript $Rscript_ENCODE_Match_TF_occupancy --ENCODE_file $ENCODE_file --Final_table_TF_motif_prediction $Final_table_TF_motif_prediction --table_sel $table_array_sel --type $type --out $output_dir --out2 $out2")
    myjobid_seff_ENCODE_Match_TF_occupancy=$(sbatch --dependency=afterany:$myjobid_ENCODE_Match_TF_occupancy --open-mode=append --output=$outfile_ENCODE_Match_TF_occupancy --job-name=$seff_name --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_ENCODE_Match_TF_occupancy >> $outfile_ENCODE_Match_TF_occupancy")




done








