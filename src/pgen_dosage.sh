
# get dosage for NK receptor SNPs

plink2 \
	--pfile ./results/pgen_NK_SNPs/all_datasets_plink_pgen_common_variants_filtered_ref_changed \
	--recode A \
	--out ./results/pgen_NK_SNPs/all_datasets_plink_pgen_common_variants_filtered_ref_changed_dosage

