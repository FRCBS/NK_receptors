# merge finnish datasets (vcf files with extracted NK SNPs)

bcftools merge /home/nihtiju/work/NK_receptors/results/extract_range/ic3_merged_chrs.vcf.gz /home/nihtiju/work/NK_receptors/results/extract_range/ic1_merged_chrs.vcf.gz /home/nihtiju/work/NK_receptors/results/extract_range/register_merged_chrs.vcf.gz /home/nihtiju/work/NK_receptors/results/extract_range/mcgill_merged_chrs.vcf.gz --force-samples -Oz -o /home/nihtiju/work/NK_receptors/results/extract_range/merged_finns.vcf.gz

bcftools index -t /home/nihtiju/work/NK_receptors/results/extract_range/merged_finns.vcf.gz



