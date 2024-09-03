# leave SNPs with promising results for Leena's donors to see their genotypes

for CHR in 1 6 7 12 18
do

	plink1.9 \
    		--bfile /home/nihtiju/work/NK_receptors/results/blood_donors/blood_donors_chr${CHR} \
    		--extract /home/nihtiju/work/NK_receptors/results/blood_donors/snps_names.txt \
    		--make-bed \
    		--out /home/nihtiju/work/NK_receptors/results/blood_donors/blood_donors_chr${CHR}_extracted

done

