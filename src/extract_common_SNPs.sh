# extract common SNPs from the individual populations & the merged population

all_datasets=("./results/pgen_NK_SNPs/newcastle" "./results/pgen_NK_SNPs/poland" "./results/pgen_NK_SNPs/katalonia" "./results/pgen_NK_SNPs/finns" "./results/pgen_NK_SNPs/all_datasets")

for i in 0 1 2 3 4
do
	
	DATASET=${all_datasets[i]}_plink_pgen
	
	plink2 \
		--pfile ${DATASET} \
		--extract ./results/common_SNPs_all_imputed.txt \
		--make-pgen \
		--out ${DATASET}_common_variants

	
done

