# extract common SNPs from the individual populations & the merged population

all_datasets=("/home/nihtiju/work/NK_receptors/results/extract_range/newcastle" "/home/nihtiju/work/NK_receptors/results/extract_range/poland" "/home/nihtiju/work/NK_receptors/results/extract_range/katalonia" "/home/nihtiju/work/NK_receptors/results/extract_range/finns" "/home/nihtiju/work/NK_receptors/results/extract_range/all_datasets")

for i in 0 1 2 3 4
do
	
	DATASET=${all_datasets[i]}_plink_pgen
	
	plink2 \
		--pfile ${DATASET} \
		--extract /home/nihtiju/work/NK_receptors/results/extract_range/common_SNPs_all_imputed.txt \
		--make-pgen \
		--out ${DATASET}_common_variants

	
done

