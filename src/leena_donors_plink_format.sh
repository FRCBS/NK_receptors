
# genotype files for Leena's donors: form vcf to plink format


for CHR in 1 2 6 11 12 16 18 19
do

	plink1.9 \
    		--vcf ./data/leena_donors/BB_finngen_R10_chr${CHR}_BLOOD_SERVICE_BB_extracted.003_2019.vcf.gz \
    		--make-bed \
    		--out ./results/leena_donors/leena_donors_chr${CHR}

done









