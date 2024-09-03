# merge all datasets (vcf files with extracted NK SNPs)
# leaving in only common SNPs between te datasets later when the files are in plink format

bcftools merge /home/nihtiju/work/NK_receptors/results/extract_range/newcastle_donors_merged_chrs.vcf.gz /home/nihtiju/work/NK_receptors/results/extract_range/poland_merged_chrs.vcf.gz /home/nihtiju/work/NK_receptors/results/extract_range/katalonia_merged_chrs.vcf.gz /home/nihtiju/work/NK_receptors/results/extract_range/newcastle_merged_chrs.vcf.gz /home/nihtiju/work/NK_receptors/results/extract_range/merged_finns.vcf.gz -Oz -o /home/nihtiju/work/NK_receptors/results/extract_range/merged_datasets.vcf.gz

bcftools index -t /home/nihtiju/work/NK_receptors/results/extract_range/merged_datasets.vcf.gz





