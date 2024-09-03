# harmonize the reference allele for all datasets, changes in the association analysis for some datasets if not done

all_datasets=("/home/nihtiju/work/NK_receptors/results/extract_range/newcastle" "/home/nihtiju/work/NK_receptors/results/extract_range/poland" "/home/nihtiju/work/NK_receptors/results/extract_range/katalonia" "/home/nihtiju/work/NK_receptors/results/extract_range/finns" "/home/nihtiju/work/NK_receptors/results/extract_range/all_datasets")


for i in 0 1 2 3 4
do
	
	DATASET=${all_datasets[i]}_plink_pgen_common_variants_filtered
		
	plink2 \
		--pfile ${DATASET} \
		--ref-allele force /home/nihtiju/work/NK_receptors/results/extract_range/change_ref_allele.txt 1 2 \
		--make-pgen \
		--out ${DATASET}_ref_changed

	
done

