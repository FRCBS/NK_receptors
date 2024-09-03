#---------------------------------------------------

# frequencies with plink2 for HSCT data

#---------------------------------------------------


all_datasets=("all_datasets" "finns" "katalonia" "newcastle" "poland")

for i in 0 1 2 3 4
do
	
	DATASET=/home/nihtiju/work/NK_receptors/results/extract_range/${all_datasets[i]}_plink_pgen_common_variants_filtered_ref_changed_newSNPsLDfiltered
	
	# allele freq
	plink2 \
		--pfile ${DATASET} \
		--keep ./results/ID_lists/train_sample_ids_donors.txt \
		--freq \
		--out /home/nihtiju/work/NK_receptors/results/frequencies/${all_datasets[i]}_allelefreq_plink2_train

	plink2 \
		--pfile ${DATASET} \
		--keep ./results/ID_lists/test_sample_ids_donors.txt \
		--freq \
		--out /home/nihtiju/work/NK_receptors/results/frequencies/${all_datasets[i]}_allelefreq_plink2_test
	
	# genotype freq
	plink2 \
		--pfile ${DATASET} \
		--keep ./results/ID_lists/train_sample_ids_donors.txt \
		--geno-counts \
		--out /home/nihtiju/work/NK_receptors/results/frequencies/${all_datasets[i]}_genotypefreq_plink2_train

	plink2 \
		--pfile ${DATASET} \
		--keep ./results/ID_lists/test_sample_ids_donors.txt \
		--geno-counts \
		--out /home/nihtiju/work/NK_receptors/results/frequencies/${all_datasets[i]}_genotypefreq_plink2_test
done

#---------------------------------------------------

# frequencies with plink1.9 for blood donor data

#---------------------------------------------------


for i in 1 6 7 12 18
do
	
	
	plink1.9 \
		--bfile /home/nihtiju/work/NK_receptors/results/blood_donors/blood_donors_chr${i}_extracted \
		--keep /home/nihtiju/work/NK_receptors/results/frequencies/in_vitro_keep.txt \
		--freq \
		--out /home/nihtiju/work/NK_receptors/results/frequencies/blood_donors_chr_${i}_allelefreq
	plink1.9 \
		--bfile /home/nihtiju/work/NK_receptors/results/blood_donors/blood_donors_chr${i}_extracted \
		--keep /home/nihtiju/work/NK_receptors/results/frequencies/in_vitro_keep.txt \
		--freqx \
		--out /home/nihtiju/work/NK_receptors/results/frequencies/blood_donors_chr_${i}_genotypefreq

done









