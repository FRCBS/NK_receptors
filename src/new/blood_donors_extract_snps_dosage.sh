# dosage for leena's donors for promising snps

for CHR in 1 6 7 12 18
do

	plink1.9 \
    		--bfile /home/nihtiju/work/NK_receptors/results/blood_donors/blood_donors_chr${CHR}_extracted \
    		--recode A \
    		--out /home/nihtiju/work/NK_receptors/results/blood_donors/blood_donors_chr${CHR}_extracted_dosage

done

