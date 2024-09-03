
# get dosage for NK receptor SNPs
# not imputation accuracy to get genotypes easily at the next step

plink2 \
	--pfile /home/nihtiju/work/NK_receptors/results/extract_range/all_datasets_plink_pgen_common_variants_filtered_ref_changed_newSNPsLDfiltered \
	--make-bed \
	--out /home/nihtiju/work/NK_receptors/results/survival/plink1

plink2 \
	--bfile /home/nihtiju/work/NK_receptors/results/survival/plink1 \
	--recode A \
	--out /home/nihtiju/work/NK_receptors/results/survival/plink1_dosage

