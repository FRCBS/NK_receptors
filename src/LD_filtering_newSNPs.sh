####################################################################################################
# extract new variants into separate file

plink2 \
	--pfile /home/nihtiju/work/NK_receptors/results/extract_range/all_datasets_plink_pgen_common_variants_filtered_ref_changed \
	--extract /home/nihtiju/work/NK_receptors/results/extract_range/SNP_names_filtered.txt \
	--make-pgen \
	--out /home/nihtiju/work/NK_receptors/results/extract_range/all_datasets_plink_pgen_common_variants_filtered_ref_changed_new

####################################################################################################
# LD-prune the new SNPs

# find LD-based variants
plink2 \
	--pfile /home/nihtiju/work/NK_receptors/results/extract_range/all_datasets_plink_pgen_common_variants_filtered_ref_changed_new \
	--indep-pairwise 500kb 0.2 \
	--out /home/nihtiju/work/NK_receptors/results/extract_range/all_datasets_plink_pgen_common_variants_filtered_ref_changed_new_LD_pruning


# exclude the LD variants from all files
all_datasets=("all_datasets" "finns" "newcastle" "katalonia" "poland")

for i in 0 1 2 3 4
do
	
	DATASET=/home/nihtiju/work/NK_receptors/results/extract_range/${all_datasets[i]}_plink_pgen_common_variants_filtered_ref_changed
	
	
	plink2 \
		--pfile ${DATASET} \
		--exclude /home/nihtiju/work/NK_receptors/results/extract_range/all_datasets_plink_pgen_common_variants_filtered_ref_changed_new_LD_pruning.prune.out \
		--make-pgen \
		--out ${DATASET}_newSNPsLDfiltered
	
done




