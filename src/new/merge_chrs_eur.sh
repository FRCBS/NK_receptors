# merge CHRs for vcf files with extracted NK SNPs

all_datasets=("/home/nihtiju/work/NK_receptors/results/extract_range/newcastle_donors" "/home/nihtiju/work/NK_receptors/results/extract_range/poland" "/home/nihtiju/work/NK_receptors/results/extract_range/katalonia" "/home/nihtiju/work/NK_receptors/results/extract_range/newcastle")

for i in 0 1 2 3
do
	
	DATASET=${all_datasets[i]}
	
	bcftools concat ${DATASET}_vcf_NK_SNPs_chr1.vcf.gz ${DATASET}_vcf_NK_SNPs_chr3.vcf.gz ${DATASET}_vcf_NK_SNPs_chr4.vcf.gz ${DATASET}_vcf_NK_SNPs_chr5.vcf.gz ${DATASET}_vcf_NK_SNPs_chr6.vcf.gz ${DATASET}_vcf_NK_SNPs_chr7.vcf.gz ${DATASET}_vcf_NK_SNPs_chr12.vcf.gz ${DATASET}_vcf_NK_SNPs_chr16.vcf.gz ${DATASET}_vcf_NK_SNPs_chr18.vcf.gz ${DATASET}_vcf_NK_SNPs_chr19.vcf.gz -Oz -o ${DATASET}_merged_chrs.vcf.gz

	bcftools index -t ${DATASET}_merged_chrs.vcf.gz
	
done





