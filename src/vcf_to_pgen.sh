# files from vcf to pgen format
# files have extracted NK SNPs

all_datasets=("/home/nihtiju/work/NK_receptors/tmp/results/extract_range/newcastle_donors_merged_chrs" "/home/nihtiju/work/NK_receptors/tmp/results/extract_range/poland_merged_chrs" "/home/nihtiju/work/NK_receptors/tmp/results/extract_range/katalonia_merged_chrs" "/home/nihtiju/work/NK_receptors/tmp/results/extract_range/merged_finns" "/home/nihtiju/work/NK_receptors/tmp/results/extract_range/merged_datasets")
save_here=("/home/nihtiju/work/NK_receptors/tmp/results/extract_range/newcastle" "/home/nihtiju/work/NK_receptors/tmp/results/extract_range/poland" "/home/nihtiju/work/NK_receptors/tmp/results/extract_range/katalonia" "/home/nihtiju/work/NK_receptors/tmp/results/extract_range/finns" "/home/nihtiju/work/NK_receptors/tmp/results/extract_range/all_datasets")


for i in 0 1 2 3 4 
do
	
	DATASET=${all_datasets[i]}
	save=${save_here[i]}
	
	# turn vcf files to plink2 format
	plink2 \
	    --vcf ${DATASET}.vcf.gz dosage=GP-force \
	    --double-id \
	    --make-pgen \
	    --out ${save}_plink_pgen
	
done










