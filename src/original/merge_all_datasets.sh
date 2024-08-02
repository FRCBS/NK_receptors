# merge all datasets (vcf files with extracted NK SNPs)
# leaving in only common SNPs between te datasets later when the files are in plink format

bcftools merge ./results/pgen_NK_SNPs/newcastle_donors_merged_chrs.vcf.gz ./results/pgen_NK_SNPs/poland_merged_chrs.vcf.gz ./results/pgen_NK_SNPs/katalonia_merged_chrs.vcf.gz ./results/pgen_NK_SNPs/newcastle_merged_chrs.vcf.gz ./results/pgen_NK_SNPs/merged_finns.vcf.gz -Oz -o ./results/pgen_NK_SNPs/merged_datasets.vcf.gz

bcftools index -t ./results/pgen_NK_SNPs/merged_datasets.vcf.gz





