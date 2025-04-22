#!/bin/bash
 
eval "$(conda shell.bash hook)"
 

Rscripts_path=$(echo "/home/manuel.tardaguila/Scripts/R/")

MASTER_ROUTE=$1
analysis=$2


output_dir=$(echo "$MASTER_ROUTE""$analysis""/")

Master_path_analysis=$(echo "$output_dir")

 
Log_files=$(echo "$Master_path_analysis""/""Log_files/")

rm -rf $Log_files
mkdir -p $Log_files

####################### Peaks are in hg19 liftOver
####################### Peaks are in hg19 liftOver
####################### Peaks are in hg19 liftOver

type=$(echo "ENCODE_Match_TF_occupancy")
outfile_ENCODE_Match_TF_occupancy=$(echo "$Log_files""outfile_1_""$type"".log")
touch $outfile_ENCODE_Match_TF_occupancy
echo -n "" > $outfile_ENCODE_Match_TF_occupancy
name_ENCODE_Match_TF_occupancy=$(echo "$type""_job")
seff_name=$(echo "seff""_""$type")


Rscript_ENCODE_Match_TF_occupancy=$(echo "$Rscripts_path""452_ENCODE_CHIPseq_controls.R")

ENCODE_file=$(echo "/group/soranzo/manuel.tardaguila/CH/reference_files/ENCODE/CHIPseq/ENCODE_CHIP_data.tsv.gz")
selected_cell_types=$(echo "hESC_H1,K562")
TF_to_search=$(echo "BCL11A,GATA6,GATA1")



myjobid_ENCODE_Match_TF_occupancy=$(sbatch --output=$outfile_ENCODE_Match_TF_occupancy --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=4096 --parsable --job-name $name_ENCODE_Match_TF_occupancy --wrap="Rscript $Rscript_ENCODE_Match_TF_occupancy --ENCODE_file $ENCODE_file --selected_cell_types $selected_cell_types --TF_to_search $TF_to_search --type $type --out $output_dir")
myjobid_seff_ENCODE_Match_TF_occupancy=$(sbatch --dependency=afterany:$myjobid_ENCODE_Match_TF_occupancy --open-mode=append --output=$outfile_ENCODE_Match_TF_occupancy --job-name=$seff_name --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_ENCODE_Match_TF_occupancy >> $outfile_ENCODE_Match_TF_occupancy")

#### bedtools_getfasta #############################


type=$(echo "bedtools_getfasta""_""$analysis")
outfile_bedtools_getfasta=$(echo "$Log_files""outfile_2_""$type"".log")
touch $outfile_bedtools_getfasta
echo -n "" > $outfile_bedtools_getfasta
name_bedtools_getfasta=$(echo "$type""_job")
seff_name=$(echo "seff""_""$type")

module load bedtools2/2.31.0

reference_genome=$(echo "/processing_data/reference_datasets/iGenomes/2022.1/Homo_sapiens/NCBI/GRCh38/Sequence/WholeGenomeFasta/genome.fa")
input_bed=$(echo "$output_dir""controls"".bed")
output_fasta=$(echo "$output_dir""intervals_controls"".fasta")

myjobid_bedtools_getfasta=$(sbatch --dependency=afterany:$myjobid_ENCODE_Match_TF_occupancy --job-name=$name_bedtools_getfasta --output=$outfile_bedtools_getfasta --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=512M --parsable --wrap="bedtools getfasta -fi $reference_genome -bed $input_bed > $output_fasta")
myjobid_seff_bedtools_getfasta=$(sbatch --dependency=afterany:$myjobid_bedtools_getfasta --open-mode=append --output=$outfile_bedtools_getfasta --job-name=$seff_name --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_bedtools_getfasta >> $outfile_bedtools_getfasta")


#### printer_fasta_files #############################


type=$(echo "printer_fasta_files""_""$analysis")
outfile_printer_fasta_files=$(echo "$Log_files""outfile_3_""$type"".log")
touch $outfile_printer_fasta_files
echo -n "" > $outfile_printer_fasta_files
name_printer_fasta_files=$(echo "$type""_job")
seff_name=$(echo "seff""_""$type")

Rscript_printer_fasta_files=$(echo "$Rscripts_path""462_Printer_of_input_fasta_files_for_CHIP_ctrls.R")

input_bed=$(echo "$output_dir""controls"".bed")
input_fasta=$(echo "$output_dir""intervals_controls"".fasta")

# --dependency=afterany:$myjobid_bedtools_getfasta

myjobid_printer_fasta_files=$(sbatch --dependency=afterany:$myjobid_bedtools_getfasta --job-name=$name_printer_fasta_files --output=$outfile_printer_fasta_files --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=2 --mem-per-cpu=1024M --parsable --wrap="Rscript $Rscript_printer_fasta_files --input_bed $input_bed --input_fasta $input_fasta --type $type --out $output_dir")
myjobid_seff_printer_fasta_files=$(sbatch --dependency=afterany:$myjobid_printer_fasta_files --open-mode=append --output=$outfile_printer_fasta_files --job-name=$seff_name --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_printer_fasta_files >> $outfile_printer_fasta_files")


conda activate my_gimme_scan


reference_TFmotif_collections=$(echo "HOMER,JASPAR2020_vertebrates")

b=($(echo "$reference_TFmotif_collections" | tr "," '\n'))

declare -a arr=()
declare -a arr2=()

for i  in "${b[@]}"
do
    reference_TFmotif_collections_sel=${i}
    #        echo "$reference_TFmotif_collections_sel"

    type=$(echo "$reference_TFmotif_collections_sel""_""gimme_scan")
    outfile_gimme_scan_REF=$(echo "$Log_files""outfile_4_""$type"".log")
    touch $outfile_gimme_scan_REF
    echo -n "" > $outfile_gimme_scan_REF
    name_gimme_scan_REF=$(echo "$type""_job")
    seff_name=$(echo "seff""_""$type")

    output_file=$(echo "$output_dir""/""$reference_TFmotif_collections_sel""_controls"".bed")
    fasta_file_for_prediction=$(echo "$output_dir""/""TF_search"".fasta")

    # --dependency=afterany:$myjobid_printer_fasta_files
    
    myjobid_gimme_scan_REF=$(sbatch --dependency=afterany:$myjobid_printer_fasta_files --job-name=$name_gimme_scan_REF --output=$outfile_gimme_scan_REF --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=2 --mem-per-cpu=1024 --parsable --wrap="gimme scan $fasta_file_for_prediction -p $reference_TFmotif_collections_sel -c 0.85 -n 20 -b > $output_file")
    myjobid_seff_gimme_scan_REF=$(sbatch --dependency=afterany:$myjobid_gimme_scan_REF --open-mode=append --output=$outfile_gimme_scan_REF --job-name=$seff_name --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_gimme_scan_REF >> $outfile_gimme_scan_REF")

    arr[${#arr[@]}]="$myjobid_gimme_scan_REF"


done

conda deactivate

done_string=$(echo "--dependency=afterany:"""""${arr[@]}"""")
dependency_string=$(echo $done_string|sed -r 's/ /:/g')
echo "$dependency_string"

#### Merge results from the TF motif predictors

type=$(echo "collect_TF_motif_results")
outfile_collect_TF_motif_results=$(echo "$Log_files""outfile_5_""$type"".log")
touch $outfile_collect_TF_motif_results
echo -n "" > $outfile_collect_TF_motif_results
name_collect_TF_motif_results=$(echo "$type""_job")
seff_name=$(echo "seff""_""$type")

Rscript_collect_TF_motif_results=$(echo "$Rscripts_path""463_TFmotif_collection_CHIP_ctrls.R")


input_file_HOMER=$(echo "$output_dir""HOMER""_controls"".bed")
input_file_JASPAR=$(echo "$output_dir""JASPAR2020_vertebrates""_controls"".bed")


# $dependency_string 

myjobid_collect_TF_motif_results=$(sbatch $dependency_string --job-name=$name_collect_TF_motif_results --output=$outfile_collect_TF_motif_results --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=4 --mem-per-cpu=4096M --parsable --wrap="Rscript $Rscript_collect_TF_motif_results --input_file_HOMER $input_file_HOMER --input_file_JASPAR $input_file_JASPAR --type $type --out $output_dir")
myjobid_seff_collect_TF_motif_results=$(sbatch --dependency=afterany:$myjobid_collect_TF_motif_results --open-mode=append --output=$outfile_collect_TF_motif_results --job-name=$seff_name --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_collect_TF_motif_results >> $outfile_collect_TF_motif_results")	

#### Final_table

type=$(echo "Final_table")
outfile_Final_table=$(echo "$Log_files""outfile_6_""$type"".log")
touch $outfile_Final_table
echo -n "" > $outfile_Final_table
name_Final_table=$(echo "$type""_job")
seff_name=$(echo "seff""_""$type")

Rscript_Final_table=$(echo "$Rscripts_path""464_select_controls_and_merge_them.R")


top_file=$(echo "$output_dir""ENCODE_top50_controls_BCL11A_GATA6_GATA1.tsv")
collect_TF_motif_results=$(echo "$output_dir""collect_TF_motif_results.tsv")
TF_search=$(echo "$output_dir""TF_search.fasta")

# --dependency=afterany:$myjobid_collect_TF_motif_results

myjobid_Final_table=$(sbatch --dependency=afterany:$myjobid_collect_TF_motif_results --job-name=$name_Final_table --output=$outfile_Final_table --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=4096M --parsable --wrap="Rscript $Rscript_Final_table --top_file $top_file --collect_TF_motif_results $collect_TF_motif_results --TF_search $TF_search --type $type --out $output_dir")
myjobid_seff_Final_table=$(sbatch --dependency=afterany:$myjobid_Final_table --open-mode=append --output=$outfile_Final_table --job-name=$seff_name --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_Final_table >> $outfile_Final_table")	



