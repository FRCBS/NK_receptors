# dosage for leena's donors for promising snps

plink1.9 \
	--bfile ./results/leena_donors/leena_donors_chr1_extracted  \
	--recode A \
	--make-bed \
	--out ./results/leena_donors/leena_donors_chr1_dosage

plink1.9 \
	--bfile ./results/leena_donors/leena_donors_chr18_extracted  \
	--recode A \
	--make-bed \
	--out ./results/leena_donors/leena_donors_chr18_dosage

