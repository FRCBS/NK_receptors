# harmonize the reference allele for all datasets, changes in the association analysis for some datasets if not done

all_datasets=("./results/pgen_NK_SNPs/newcastle" "./results/pgen_NK_SNPs/poland" "./results/pgen_NK_SNPs/katalonia" "./results/pgen_NK_SNPs/finns" "./results/pgen_NK_SNPs/all_datasets")


for i in 0 1 2 3 4
do
	
	DATASET=${all_datasets[i]}_plink_pgen_common_variants_filtered
		
	plink2 \
		--pfile ${DATASET} \
		--ref-allele force ./results/change_ref_allele.txt 1 2 \
		--make-pgen \
		--out ${DATASET}_ref_changed

	
done

