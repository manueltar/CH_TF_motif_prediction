#!/bin/bash

MASTER_ROUTE=$1
analysis=$2
upstream_span=$3
downstream_span=$4
spanning_of_motif=$upstream_span


Rscripts_path=$(echo "/home/manuel.tardaguila/Scripts/R/")
module load R/4.1.0

eval "$(conda shell.bash hook)"


output_dir=$(echo "$MASTER_ROUTE""$analysis""/")

Log_files=$(echo "$output_dir""/""Log_files/")

rm -rf $Log_files
mkdir -p $Log_files

########################################################################################
table_array=$(echo 'FULL,NoMosaic,tet2,DNMT')


a=($(echo "$table_array" | tr "," '\n'))


declare -a global=()

global_length=${#a[@]}

for (( iteration=0; iteration<${global_length}; iteration=iteration+1 ));
do

    table_array_sel=${a[$iteration]}

    echo "------------------------------------------------------------------------------------------------------------->$table_array_sel"

    file=$(echo "$MASTER_ROUTE""annotated_results_""$table_array_sel""_PP01_MAF.txt")

    echo "---------------------------->$file"

    #### print_bed #############################


    type=$(echo "print_bed""_""$table_array_sel""_""$analysis")
    outfile_print_bed=$(echo "$Log_files""outfile_1_""$type"".log")
    touch $outfile_print_bed
    echo -n "" > $outfile_print_bed
    name_print_bed=$(echo "$type""_job")
    seff_name=$(echo "seff""_""$type")

    Rscript_print_bed=$(echo "$Rscripts_path""426_CH_predict_motifs_print_bed_v2.R")


 
    myjobid_print_bed=$(sbatch --job-name=$name_print_bed --output=$outfile_print_bed --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=2 --mem-per-cpu=1024M --parsable --wrap="Rscript $Rscript_print_bed --michelas_table $file --type $type --upstream_span $upstream_span --downstream_span $downstream_span --table_sel $table_array_sel --out $output_dir")
    myjobid_seff_print_bed=$(sbatch --dependency=afterany:$myjobid_print_bed --open-mode=append --output=$outfile_print_bed --job-name=$seff_name --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_print_bed >> $outfile_print_bed")


    #### bedtools_getfasta #############################


    type=$(echo "bedtools_getfasta""_""$table_array_sel""_""$analysis")
    outfile_bedtools_getfasta=$(echo "$Log_files""outfile_2_""$type"".log")
    touch $outfile_bedtools_getfasta
    echo -n "" > $outfile_bedtools_getfasta
    name_bedtools_getfasta=$(echo "$type""_job")
    seff_name=$(echo "seff""_""$type")

    module load bedtools2/2.31.0

    reference_genome=$(echo "/processing_data/reference_datasets/iGenomes/2022.1/Homo_sapiens/NCBI/GRCh38/Sequence/WholeGenomeFasta/genome.fa")
    input_bed=$(echo "$output_dir""VARS_""$table_array_sel"".bed")
    output_fasta=$(echo "$output_dir""intervals_""$table_array_sel"".fasta")

    myjobid_bedtools_getfasta=$(sbatch --dependency=afterany:$myjobid_print_bed --job-name=$name_bedtools_getfasta --output=$outfile_bedtools_getfasta --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=512M --parsable --wrap="bedtools getfasta -fi $reference_genome -bed $input_bed > $output_fasta")
    myjobid_seff_bedtools_getfasta=$(sbatch --dependency=afterany:$myjobid_bedtools_getfasta --open-mode=append --output=$outfile_bedtools_getfasta --job-name=$seff_name --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_bedtools_getfasta >> $outfile_bedtools_getfasta")


    #### printer_fasta_files #############################


    type=$(echo "printer_fasta_files""_""$table_array_sel""_""$analysis")
    outfile_printer_fasta_files=$(echo "$Log_files""outfile_3_""$type"".log")
    touch $outfile_printer_fasta_files
    echo -n "" > $outfile_printer_fasta_files
    name_printer_fasta_files=$(echo "$type""_job")
    seff_name=$(echo "seff""_""$type")

    Rscript_printer_fasta_files=$(echo "$Rscripts_path""427_Printer_of_input_fasta_files_v2.R")

    input_bed=$(echo "$output_dir""VARS_""$table_array_sel"".bed")
    input_fasta=$(echo "$output_dir""intervals_""$table_array_sel"".fasta")


    myjobid_printer_fasta_files=$(sbatch --dependency=afterany:$myjobid_bedtools_getfasta --job-name=$name_printer_fasta_files --output=$outfile_printer_fasta_files --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=2 --mem-per-cpu=1024M --parsable --wrap="Rscript $Rscript_printer_fasta_files --input_bed $input_bed --input_fasta $input_fasta --upstream_span $upstream_span --downstream_span $downstream_span --table_sel $table_array_sel --type $type --out $output_dir")
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

	type=$(echo "$table_array_sel""_""REF_""$reference_TFmotif_collections_sel""_""gimme_scan")
	outfile_gimme_scan_REF=$(echo "$Log_files""outfile_4_""$type"".log")
	touch $outfile_gimme_scan_REF
	echo -n "" > $outfile_gimme_scan_REF
	name_gimme_scan_REF=$(echo "$type""_job")
	seff_name=$(echo "seff""_""$type")

	output_file=$(echo "$output_dir""/""$reference_TFmotif_collections_sel""_REF_""$table_array_sel"".bed")
	fasta_file_for_prediction=$(echo "$output_dir""/""TF_search_REF_Allele_""$table_array_sel"".fasta")

	# --dependency=afterany:$myjobid_Print_the_sequences_to_search_for_motifs
	
	myjobid_gimme_scan_REF=$(sbatch --dependency=afterany:$myjobid_printer_fasta_files --job-name=$name_gimme_scan_REF --output=$outfile_gimme_scan_REF --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=2 --mem-per-cpu=1024 --parsable --wrap="gimme scan $fasta_file_for_prediction -p $reference_TFmotif_collections_sel -c 0.85 -n 20 -b > $output_file")
	myjobid_seff_gimme_scan_REF=$(sbatch --dependency=afterany:$myjobid_gimme_scan_REF --open-mode=append --output=$outfile_gimme_scan_REF --job-name=$seff_name --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_gimme_scan_REF >> $outfile_gimme_scan_REF")

	arr[${#arr[@]}]="$myjobid_gimme_scan_REF"


	################################################################################################################################ REF - ALT split #####################################################################################################################################################################
	################################################################################################################################ REF - ALT split #####################################################################################################################################################################
	################################################################################################################################ REF - ALT split #####################################################################################################################################################################

	reference_TFmotif_collections_sel=${i}
#	echo "$reference_TFmotif_collections_sel"

	type=$(echo "$table_array_sel""_""ALT_""$reference_TFmotif_collections_sel""_""gimme_scan")
	outfile_gimme_scan_ALT=$(echo "$Log_files""outfile_6_""$type"".log")
	touch $outfile_gimme_scan_ALT
	echo -n "" > $outfile_gimme_scan_ALT
	name_gimme_scan_ALT=$(echo "$type""_job")
	seff_name=$(echo "seff""_""$type")

	output_file=$(echo "$output_dir""/""$reference_TFmotif_collections_sel""_ALT_""$table_array_sel"".bed")
	fasta_file_for_prediction=$(echo "$output_dir""/""TF_search_ALT_Allele_""$table_array_sel"".fasta")

	# --dependency=afterany:$myjobid_Print_the_sequences_to_search_for_motifs
	
	myjobid_gimme_scan_ALT=$(sbatch --dependency=afterany:$myjobid_printer_fasta_files --job-name=$name_gimme_scan_ALT --output=$outfile_gimme_scan_ALT --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=2 --mem-per-cpu=1024 --parsable --wrap="gimme scan $fasta_file_for_prediction -p $reference_TFmotif_collections_sel -c 0.85 -n 20 -b > $output_file")
	myjobid_seff_gimme_scan_ALT=$(sbatch --dependency=afterany:$myjobid_gimme_scan_ALT --open-mode=append --output=$outfile_gimme_scan_ALT --job-name=$seff_name --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_gimme_scan_ALT >> $outfile_gimme_scan_ALT")


	arr2[${#arr2[@]}]="$myjobid_gimme_scan_ALT"
	

    done

    conda deactivate

    done_string=$(echo "--dependency=afterany:"""""${arr[@]}"""")
   # echo "REF: $done_string"

    done_string=$(echo "--dependency=afterany:"""""${arr2[@]}"""")
   # echo "ALT: $done_string"



    concatenated=( ${arr[*]} ${arr2[*]} )

    done_string=$(echo "--dependency=afterany:"""""${concatenated[@]}"""")


  #  echo "concatenated: $done_string"

    dependency_string=$(echo $done_string|sed -r 's/ /:/g')

    echo "$dependency_string"


    allele_array=$(echo "REF,ALT")

    b=($(echo "$allele_array" | tr "," '\n'))


    declare -a arr2=()

    for i  in "${b[@]}"
    do
        allele_array_sel=${i}
        echo "$allele_array_sel"

	#### Merge results from the TF motif predictors

	eval "$(conda shell.bash hook)"

	type=$(echo "collect_TF_motif_results""_""$allele_array_sel""_""$table_array_sel")
	outfile_collect_TF_motif_results=$(echo "$Log_files""outfile_8_""$type"".log")
	touch $outfile_collect_TF_motif_results
	echo -n "" > $outfile_collect_TF_motif_results
	name_collect_TF_motif_results=$(echo "$type""_job")
	seff_name=$(echo "seff""_""$type")

	Rscript_collect_TF_motif_results=$(echo "$Rscripts_path""428_TFmotif_collection_only_v2_michelas_snps_v2.R")


	input_file_HOMER=$(echo "$output_dir""HOMER""_""$allele_array_sel""_""$table_array_sel"".bed")
	input_file_JASPAR=$(echo "$output_dir""JASPAR2020_vertebrates""_""$allele_array_sel""_""$table_array_sel"".bed")


	# $dependency_string


	myjobid_collect_TF_motif_results=$(sbatch $dependency_string --job-name=$name_collect_TF_motif_results --output=$outfile_collect_TF_motif_results --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=4 --mem-per-cpu=4096M --parsable --wrap="Rscript $Rscript_collect_TF_motif_results --input_file_HOMER $input_file_HOMER --input_file_JASPAR $input_file_JASPAR --spanning_of_motif $spanning_of_motif --allele $allele_array_sel --table_sel $table_array_sel --type $type --out $output_dir")
	myjobid_seff_collect_TF_motif_results=$(sbatch --dependency=afterany:$myjobid_collect_TF_motif_results --open-mode=append --output=$outfile_collect_TF_motif_results --job-name=$seff_name --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_collect_TF_motif_results >> $outfile_collect_TF_motif_results")	

	echo "->>>$myjobid_collect_TF_motif_results"
	arr2[${#arr2[@]}]="$myjobid_collect_TF_motif_results"


    done


    done_string=$(echo "--dependency=afterany:"""""${arr2[@]}"""")
    dependency_string_array=$(echo $done_string|sed -r 's/ /:/g')

    echo "FINAL: $dependency_string_array"

    #### Merge results from the TF motif predictors

    eval "$(conda shell.bash hook)"

    type=$(echo "final_TF_table")
    outfile_final_TF_table=$(echo "$Log_files""outfile_9_""$type"".log")
    touch $outfile_final_TF_table
    echo -n "" > $outfile_final_TF_table
    name_final_TF_table=$(echo "$type""_job")
    seff_name=$(echo "seff""_""$type")

    Rscript_final_TF_table=$(echo "$Rscripts_path""429_Final_table_of_TF_motifs_predicted_v2.R")


    input_REF=$(echo "$output_dir""collect_TF_motif_results_REF_""$table_array_sel"".tsv")
    input_ALT=$(echo "$output_dir""collect_TF_motif_results_ALT_""$table_array_sel"".tsv")

    # $dependency_string_array

    myjobid_final_TF_table=$(sbatch $dependency_string_array --job-name=$name_final_TF_table --output=$outfile_final_TF_table --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=4 --mem-per-cpu=4096M --parsable --wrap="Rscript $Rscript_final_TF_table --input_REF $input_REF --input_ALT $input_ALT --table_sel $table_array_sel --type $type --out $output_dir")
    myjobid_seff_final_TF_table=$(sbatch --dependency=afterany:$myjobid_final_TF_table --open-mode=append --output=$outfile_final_TF_table --job-name=$seff_name --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_final_TF_table >> $outfile_final_TF_table")	


done



