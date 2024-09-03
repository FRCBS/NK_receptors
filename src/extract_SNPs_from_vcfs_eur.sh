# extract NK SNPs from vcf files

all_datasets=("/home/nihtiju/work/HSCT_predictor/results/imputation/newcastle_donors/newcastle_donors_hg38_nonreference_flipped_ambiguous_flipped_vcf_cleaned_vcf" "/home/nihtiju/work/HSCT_predictor/results/imputation/poland/poland_hg38_nonreference_flipped_ambiguous_flipped_vcf_cleaned_vcf" "/home/nihtiju/work/HSCT_predictor/results/imputation/katalonia/imputation_salla/IC2_3S_genotyped_cohort_hg38_cleaned_vcf" "/home/nihtiju/work/HSCT_predictor/results/imputation/newcastle/imputation_salla/VPU_ILLUMINA_NOV_2019_hg38_cleaned_vcf")
save_here=("/home/nihtiju/work/NK_receptors/results/extract_range/newcastle_donors" "/home/nihtiju/work/NK_receptors/results/extract_range/poland" "/home/nihtiju/work/NK_receptors/results/extract_range/katalonia" "/home/nihtiju/work/NK_receptors/results/extract_range/newcastle")

for i in 0 1 2 3
do
	
	DATASET=${all_datasets[i]}
	DATASET=${DATASET}_imputed_info
	save=${save_here[i]}
	
	for CHR in 1 3 4 5 6 7 12 16 18 19
	do
		    
		# extract NK SNPs
		
		bcftools view --include ID==@/home/nihtiju/work/NK_receptors/results/extract_range/SNP_names_old_and_new.txt ${DATASET}_chr${CHR}.vcf.gz -Oz -o ${save}_vcf_NK_SNPs_chr${CHR}.vcf.gz
		
		# make index files
		bcftools index -t ${save}_vcf_NK_SNPs_chr${CHR}.vcf.gz
		    
	done
	
done






