# merge CHRs for vcf files with extracted NK SNPs

all_datasets=("./results/pgen_NK_SNPs/newcastle_donors" "./results/pgen_NK_SNPs/poland" "./results/pgen_NK_SNPs/katalonia" "./results/pgen_NK_SNPs/newcastle")

for i in 0 1 2 3
do
	
	DATASET=${all_datasets[i]}
	
	bcftools concat ${DATASET}_vcf_NK_SNPs_chr1.vcf.gz ${DATASET}_vcf_NK_SNPs_chr2.vcf.gz ${DATASET}_vcf_NK_SNPs_chr6.vcf.gz ${DATASET}_vcf_NK_SNPs_chr11.vcf.gz ${DATASET}_vcf_NK_SNPs_chr12.vcf.gz ${DATASET}_vcf_NK_SNPs_chr16.vcf.gz ${DATASET}_vcf_NK_SNPs_chr18.vcf.gz ${DATASET}_vcf_NK_SNPs_chr19.vcf.gz -Oz -o ${DATASET}_merged_chrs.vcf.gz

	bcftools index -t ${DATASET}_merged_chrs.vcf.gz
	
done






