# filter out individuals

all_datasets=("/home/nihtiju/work/NK_receptors/results/extract_range/newcastle" "/home/nihtiju/work/NK_receptors/results/extract_range/poland" "/home/nihtiju/work/NK_receptors/results/extract_range/katalonia" "/home/nihtiju/work/NK_receptors/results/extract_range/finns" "/home/nihtiju/work/NK_receptors/results/extract_range/all_datasets")


for i in 0 1 2 3 4
do
	
	DATASET=${all_datasets[i]}_plink_pgen_common_variants
	
	plink2 \
		--pfile ${DATASET} \
		--remove ./results/extract_range/individuals_to_remove.txt \
		--make-pgen \
		--out ${DATASET}_filtered

	
done  

