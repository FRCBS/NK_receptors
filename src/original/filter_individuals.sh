# filter out individuals

all_datasets=("./results/pgen_NK_SNPs/newcastle" "./results/pgen_NK_SNPs/poland" "./results/pgen_NK_SNPs/katalonia" "./results/pgen_NK_SNPs/finns" "./results/pgen_NK_SNPs/all_datasets")


for i in 0 1 2 3 4
do
	
	DATASET=${all_datasets[i]}_plink_pgen_common_variants
	
	plink2 \
		--pfile ${DATASET} \
		--remove ./results/pgen_NK_SNPs/individuals_to_remove.txt \
		--make-pgen \
		--out ${DATASET}_filtered

	
done  

