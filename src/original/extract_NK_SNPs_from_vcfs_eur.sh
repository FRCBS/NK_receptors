# extract NK SNPs from vcf files

all_datasets=("/home/nithiju/work/HSCT_predictor/results/imputation/newcastle_donors/newcastle_donors_hg38_nonreference_flipped_ambiguous_flipped_vcf_cleaned_vcf" "/home/nithiju/work/HSCT_predictor/results/imputation/poland/poland_hg38_nonreference_flipped_ambiguous_flipped_vcf_cleaned_vcf" "/home/nithiju/work/HSCT_predictor/results/imputation/katalonia/imputation_salla/IC2_3S_genotyped_cohort_hg38_cleaned_vcf" "/home/nithiju/work/HSCT_predictor/results/imputation/newcastle/imputation_salla/VPU_ILLUMINA_NOV_2019_hg38_cleaned_vcf")
save_here=("./results/pgen_NK_SNPs/newcastle_donors" "./results/pgen_NK_SNPs/poland" "./results/pgen_NK_SNPs/katalonia" "./results/pgen_NK_SNPs/newcastle")

for i in 0 1 2 3
do
	
	DATASET=${all_datasets[i]}
	DATASET=${DATASET}_imputed_info
	save=${save_here[i]}
	
	for CHR in 1 2 6 11 12 16 18 19
	do
		    
		# extract NK SNPs
		
		bcftools view --include ID==@./results/SNP_names_receptors_all.txt ${DATASET}_chr${CHR}.vcf.gz -Oz -o ${save}_vcf_NK_SNPs_chr${CHR}.vcf.gz
		
		# make index files
		bcftools index -t ${save}_vcf_NK_SNPs_chr${CHR}.vcf.gz
		    
	done
	
done








