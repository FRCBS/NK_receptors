# extract NK SNPs from vcf files

all_datasets=("/media/volume/datavolume/imputation/ic3/ic3_b38" "/media/volume/datavolume/imputation/ic1/ic1_b38" "/media/volume/datavolume/imputation/register/register_b38" "/media/volume/datavolume/imputation/McGill/McGill_b38_chr_and_positions_cleaned_no_duplicates_vcf")

save_here=("/media/volume/datavolume/NK_cell_markers_and_relapse/NK_SNPs_pgen/ic3" "/media/volume/datavolume/NK_cell_markers_and_relapse/NK_SNPs_pgen/ic1" "/media/volume/datavolume/NK_cell_markers_and_relapse/NK_SNPs_pgen/register" "/media/volume/datavolume/NK_cell_markers_and_relapse/NK_SNPs_pgen/mcgill")

for i in 0 1 2 3
do
	
	DATASET=${all_datasets[i]}
	DATASET=${DATASET}_imputed_info
	save=${save_here[i]}
	
	for CHR in 1 2 6 11 12 16 18 19
	do
		    
		# extract NK SNPs
		
		bcftools view --include ID==@/media/volume/datavolume/NK_cell_markers_and_relapse/SNP_names_receptors_all.txt ${DATASET}_chr${CHR}.vcf.gz -Oz -o ${save}_vcf_NK_SNPs_chr${CHR}.vcf.gz
		
		# make index files
		bcftools index -t ${save}_vcf_NK_SNPs_chr${CHR}.vcf.gz
		    
	done
	
done








