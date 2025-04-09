####################################################	Functional Annotation of the snps in the 95% CS of the CHIP GWAS ########################################################################################################
####################################################	Functional Annotation of the snps in the 95% CS of the CHIP GWAS ########################################################################################################
####################################################	Functional Annotation of the snps in the 95% CS of the CHIP GWAS ########################################################################################################

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

# 6 to get TF occupancy from Blood CHIP-studies:

Download the peaks file from the ChiP_ATLAS: https://github.com/inutano/chip-atlas/wiki#downloads_doc

$ wget https://chip-atlas.dbcls.jp/data/hg19/allPeaks_light/allPeaks_light.hg19.10.bed.gz

and run

$ nohup bash ~/Scripts/Wraper_scripts/156_CHIP_atalas_TF_occupancy_CS.sh /group/soranzo/manuel.tardaguila/CH/ALL_variants_in_CS/TF_motif_analysis/ TF_occupancy_filter &

# 5 Find overlap with ATAC data with. ATAC-seq data set: (bulk) https://github.com/caleblareau/singlecell_bloodtraits/tree/master/data/bulk/ATAC 29August2017_EJCsamples_allReads_500bp.bed and 29August2017_EJCsamples_allReads_500bp.counts.txt from "Interrogation of human hematopoiesis at single-cell and single-variant resolution" (Ulirsch et al 2019). It is in hg19.

$ bash ~/Scripts/Wraper_scripts/155_Functional_annotation_intersect_with_ATAC.sh /group/soranzo/manuel.tardaguila/CH/ALL_variants_in_CS/ ATACseq_analysis

# 6 Find overlap with DNase data from ENCODE in K-562, HL60 and H1hesc cells

$ bash ~/Scripts/Wraper_scripts/157_Functional_annotation_intersect_with_DNase.sh /group/soranzo/manuel.tardaguila/CH/ALL_variants_in_CS/ DNase_ENCODE
