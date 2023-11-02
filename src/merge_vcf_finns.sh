# merge finnish datasets (vcf files with extracted NK SNPs)

bcftools merge ./results/pgen_NK_SNPs/ic3_merged_chrs.vcf.gz ./results/pgen_NK_SNPs/ic1_merged_chrs.vcf.gz ./results/pgen_NK_SNPs/register_merged_chrs.vcf.gz ./results/pgen_NK_SNPs/mcgill_merged_chrs.vcf.gz --force-samples -Oz -o ./results/pgen_NK_SNPs/merged_finns.vcf.gz

bcftools index -t ./results/pgen_NK_SNPs/merged_finns.vcf.gz



