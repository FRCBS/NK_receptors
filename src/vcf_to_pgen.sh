# files from vcf to pgen format
# files have extracted NK SNPs

all_datasets=("./results/pgen_NK_SNPs/newcastle_donors_merged_chrs" "./results/pgen_NK_SNPs/poland_merged_chrs" "./results/pgen_NK_SNPs/katalonia_merged_chrs" "./results/pgen_NK_SNPs/merged_finns" "./results/pgen_NK_SNPs/merged_datasets")
save_here=("./results/pgen_NK_SNPs/newcastle" "./results/pgen_NK_SNPs/poland" "./results/pgen_NK_SNPs/katalonia" "./results/pgen_NK_SNPs/finns" "./results/pgen_NK_SNPs/all_datasets")


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










