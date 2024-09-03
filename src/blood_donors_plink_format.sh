
# genotype files for Leena's blood donors: form vcf to plink format


for CHR in {1..23}
do

	plink1.9 \
    		--vcf /home/nihtiju/work/NK_receptors/data/blood_donors/BB_finngen_R12_chr${CHR}_BLOOD_SERVICE_BB_extracted.003_2019.vcf.gz \
    		--make-bed \
    		--out /home/nihtiju/work/NK_receptors/results/blood_donors/blood_donors_chr${CHR}

done









