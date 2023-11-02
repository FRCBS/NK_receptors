# leave SNPs with promising results for Leena's donors to see their genotypes

plink1.9 \
	--bfile ./results/leena_donors/leena_donors_chr1 \
	--extract ./results/leena_donors/snps_names.txt \
	--make-bed \
	--out ./results/leena_donors/leena_donors_chr1_extracted
	
plink1.9 \
	--bfile ./results/leena_donors/leena_donors_chr18 \
	--extract ./results/leena_donors/snps_names.txt \
	--make-bed \
	--out ./results/leena_donors/leena_donors_chr18_extracted

