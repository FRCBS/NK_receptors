# merge CHRs for vcf files with extracted NK SNPs

all_datasets=("/media/volume/datavolume/NK_cell_markers_and_relapse/NK_SNPs_pgen/ic3" "/media/volume/datavolume/NK_cell_markers_and_relapse/NK_SNPs_pgen/ic1" "/media/volume/datavolume/NK_cell_markers_and_relapse/NK_SNPs_pgen/register" "/media/volume/datavolume/NK_cell_markers_and_relapse/NK_SNPs_pgen/mcgill")


for i in 0 1 2 3
do
	
	DATASET=${all_datasets[i]}
	
	bcftools concat ${DATASET}_vcf_NK_SNPs_chr1.vcf.gz ${DATASET}_vcf_NK_SNPs_chr2.vcf.gz ${DATASET}_vcf_NK_SNPs_chr6.vcf.gz ${DATASET}_vcf_NK_SNPs_chr11.vcf.gz ${DATASET}_vcf_NK_SNPs_chr12.vcf.gz ${DATASET}_vcf_NK_SNPs_chr16.vcf.gz ${DATASET}_vcf_NK_SNPs_chr18.vcf.gz ${DATASET}_vcf_NK_SNPs_chr19.vcf.gz -Oz -o ${DATASET}_merged_chrs.vcf.gz

	bcftools index -t ${DATASET}_merged_chrs.vcf.gz
	
done



