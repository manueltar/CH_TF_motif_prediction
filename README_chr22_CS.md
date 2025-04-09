# 1 copy all the files of a fine-mapping from ssu. In this case all of the ones in chr22

$ find /ssu/bsssu/michelas_finemapping/results/finemap_w_additional_info/ -name "*chr22*" -exec cp {} /group/soranzo/manuel.tardaguila/CH/ALL_variants_in_CS/ \;

# 2 run this bash script to map the GWAS traits and solve a0, a1 into ref and alt for all the variants in the CS

$ bash ~/Scripts/Wraper_scripts/149_ALL_CS_variants_TF_motif.sh /group/soranzo/manuel.tardaguila/CH/ ALL_variants_in_CS


# 3 run this to get all the TF motifs

$ bash ~/Scripts/Wraper_scripts/150_TF_motif_prediction_CS.sh /group/soranzo/manuel.tardaguila/CH/ALL_variants_in_CS/ TF_motif_analysis 15 15

# 4 to get specific breakups of motif, identify TF in REF and TF in ALT of interest and run:

$ bash ~/Scripts/Wraper_scripts/154_Bespoke_TF_break_CH_ALL_CS.sh /group/soranzo/manuel.tardaguila/CH/ALL_variants_in_CS/ TF_motif_analysis <TF_REF> <TF_ALT>

TF_REF = BCL11A
TF_ALT = GATA6

# 5 find overlap with ATAC data with

$ bash ~/Scripts/Wraper_scripts/155_Functional_annotation_intersect_with_ATAC.sh /group/soranzo/manuel.tardaguila/CH/ALL_variants_in_CS/ ATACseq_analysis

