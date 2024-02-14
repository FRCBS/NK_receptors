
# get dosage for NK receptor SNPs
# not imputation accuracy to get genotypes easily at the next step

plink2 \
	--pfile ./results/pgen_NK_SNPs/all_datasets_plink_pgen_common_variants_filtered_ref_changed \
	--make-bed \
	--out ./results/survival/plink1

plink2 \
	--bfile ./results/survival/plink1 \
	--recode A \
	--out ./results/survival/dosage

