# association for NK receptor SNPs using pgen files

# for training set

# association for NK receptor SNPs using pgen files

function association_binary() {

	genotype_file=$1
	keep_individuals=$2
	save_to=$3
	keep_individuals_relapse=$4
	
    	# association for all aGvHD
    	plink2 \
		--pfile ${genotype_file} \
		--keep ${keep_individuals} \
		--glm omit-ref \
		--1 \
		--ci 0.95 \
		--adjust \
		--covar-variance-standardize \
		--covar /home/nihtiju/work/NK_receptors/results/pheno_and_covars/covars_donors_aGvHD_relapse_all_binary.txt \
		--pheno ./results/pheno_and_covars/pheno_separate_aGvHD_plink.txt \
		--out ${save_to}_aGvHD
	
    	# association for cGvHD
    	plink2 \
		--pfile ${genotype_file} \
		--keep ${keep_individuals} \
		--glm omit-ref \
		--1 \
		--ci 0.95 \
		--adjust \
		--covar-variance-standardize \
		--covar /home/nihtiju/work/NK_receptors/results/pheno_and_covars/covars_donors_cGvHD_all_binary.txt \
		--pheno ./results/pheno_and_covars/pheno_separate_cGvHD_plink.txt \
		--out ${save_to}_cGvHD
	
	# association for relapse
    	plink2 \
		--pfile ${genotype_file} \
		--keep ${keep_individuals_relapse} \
		--glm omit-ref \
		--1 \
		--ci 0.95 \
		--adjust \
		--covar-variance-standardize \
		--covar /home/nihtiju/work/NK_receptors/results/pheno_and_covars/covars_donors_aGvHD_relapse_all_binary.txt \
		--pheno ./results/pheno_and_covars/pheno_separate_relapse_plink.txt \
		--out ${save_to}_relapse
    	
}


function association_newcastle_binary() {

	genotype_file=$1
	keep_individuals=$2
	save_to=$3
	keep_individuals_relapse=$4
	
    	# association for aGvHD
    	plink2 \
	--pfile ${genotype_file} \
	--keep ${keep_individuals} \
	--glm omit-ref \
	--1 \
	--ci 0.95 \
	--adjust \
	--covar-variance-standardize \
	--covar /home/nihtiju/work/NK_receptors/results/pheno_and_covars/covars_donors_aGvHD_relapse_all_binary_newcastle.txt \
	--pheno ./results/pheno_and_covars/pheno_separate_aGvHD_plink.txt \
	--out ${save_to}_aGvHD
	
    	# association for cGvHD
    	plink2 \
	--pfile ${genotype_file} \
	--keep ${keep_individuals} \
	--glm omit-ref \
	--1 \
	--ci 0.95 \
	--adjust \
	--covar-variance-standardize \
	--covar /home/nihtiju/work/NK_receptors/results/pheno_and_covars/covars_donors_cGvHD_all_binary_newcastle.txt \
	--pheno ./results/pheno_and_covars/pheno_separate_cGvHD_plink.txt \
	--out ${save_to}_cGvHD
	
	# association for relapse
    	plink2 \
	--pfile ${genotype_file} \
	--keep ${keep_individuals_relapse} \
	--glm omit-ref \
	--1 \
	--ci 0.95 \
	--adjust \
	--covar-variance-standardize \
	--covar /home/nihtiju/work/NK_receptors/results/pheno_and_covars/covars_donors_aGvHD_relapse_all_binary_newcastle.txt \
	--pheno ./results/pheno_and_covars/pheno_separate_relapse_plink.txt \
	--out ${save_to}_relapse
    	
}

function association_poland_binary() {

	genotype_file=$1
	keep_individuals=$2
	save_to=$3
	keep_individuals_relapse=$4
	
    	# association for aGvHD
    	plink2 \
	--pfile ${genotype_file} \
	--keep ${keep_individuals} \
	--glm omit-ref \
	--1 \
	--ci 0.95 \
	--adjust \
	--covar-variance-standardize \
	--covar /home/nihtiju/work/NK_receptors/results/pheno_and_covars/covars_donors_aGvHD_relapse_all_binary_poland.txt \
	--pheno ./results/pheno_and_covars/pheno_separate_aGvHD_POLAND_plink.txt \
	--out ${save_to}_aGvHD
	
    	# association for cGvHD
    	plink2 \
	--pfile ${genotype_file} \
	--keep ${keep_individuals} \
	--glm omit-ref \
	--1 \
	--ci 0.95 \
	--adjust \
	--covar-variance-standardize \
	--covar /home/nihtiju/work/NK_receptors/results/pheno_and_covars/covars_donors_cGvHD_all_binary_poland.txt \
	--pheno ./results/pheno_and_covars/pheno_separate_cGvHD_POLAND_plink.txt \
	--out ${save_to}_cGvHD
	
	# association for relapse
    	plink2 \
	--pfile ${genotype_file} \
	--keep ${keep_individuals_relapse} \
	--glm omit-ref \
	--1 \
	--ci 0.95 \
	--adjust \
	--covar-variance-standardize \
	--covar /home/nihtiju/work/NK_receptors/results/pheno_and_covars/covars_donors_aGvHD_relapse_all_binary_poland.txt \
	--pheno ./results/pheno_and_covars/pheno_separate_relapse_plink.txt \
	--out ${save_to}_relapse
    	
}

#-------------------------------------------------------------------------------------------------------------------------------------------------------

# all populations: 
association_binary "/home/nihtiju/work/NK_receptors/results/extract_range/all_datasets_plink_pgen_common_variants_filtered_ref_changed_newSNPsLDfiltered" "./results/ID_lists/train_sample_ids_donors.txt" "/home/nihtiju/work/NK_receptors/results/association_LD_newOnly/train/all_train/association_all_train" "./results/ID_lists/donors_train_relapse.txt"

# fin:
association_binary "/home/nihtiju/work/NK_receptors/results/extract_range/finns_plink_pgen_common_variants_filtered_ref_changed_newSNPsLDfiltered" "./results/ID_lists/train_sample_ids_donors.txt" "/home/nihtiju/work/NK_receptors/results/association_LD_newOnly/train/finns_train/association_finns_train" "./results/ID_lists/donors_train_relapse.txt"

# kat :
association_binary "/home/nihtiju/work/NK_receptors/results/extract_range/katalonia_plink_pgen_common_variants_filtered_ref_changed_newSNPsLDfiltered" "./results/ID_lists/train_sample_ids_donors.txt" "/home/nihtiju/work/NK_receptors/results/association_LD_newOnly/train/katalonia_train/association_katalonia_train" "./results/ID_lists/donors_train_relapse.txt"

# nc:
association_newcastle_binary "/home/nihtiju/work/NK_receptors/results/extract_range/newcastle_plink_pgen_common_variants_filtered_ref_changed_newSNPsLDfiltered" "./results/ID_lists/train_sample_ids_donors.txt" "/home/nihtiju/work/NK_receptors/results/association_LD_newOnly/train/newcastle_train/association_newcastle_train" "./results/ID_lists/donors_train_relapse.txt"

# poland with phenotype file without cGvHD severe
association_poland_binary "/home/nihtiju/work/NK_receptors/results/extract_range/poland_plink_pgen_common_variants_filtered_ref_changed_newSNPsLDfiltered" "./results/ID_lists/train_sample_ids_donors.txt" "/home/nihtiju/work/NK_receptors/results/association_LD_newOnly/train/poland_train/association_poland_train" "./results/ID_lists/donors_train_relapse.txt"



##################################################################################################################################################################

# for test set

#-------------------------------------------------------------------------------------------------------------------------------------------------------

# all populations: 
association_binary "/home/nihtiju/work/NK_receptors/results/extract_range/all_datasets_plink_pgen_common_variants_filtered_ref_changed_newSNPsLDfiltered" "./results/ID_lists/test_sample_ids_donors.txt" "/home/nihtiju/work/NK_receptors/results/association_LD_newOnly/test/all_test/association_all_test" "./results/ID_lists/donors_test_relapse.txt"

# fin:
association_binary "/home/nihtiju/work/NK_receptors/results/extract_range/finns_plink_pgen_common_variants_filtered_ref_changed_newSNPsLDfiltered" "./results/ID_lists/test_sample_ids_donors.txt" "/home/nihtiju/work/NK_receptors/results/association_LD_newOnly/test/finns_test/association_finns_test" "./results/ID_lists/donors_test_relapse.txt"

# kat :
association_binary "/home/nihtiju/work/NK_receptors/results/extract_range/katalonia_plink_pgen_common_variants_filtered_ref_changed_newSNPsLDfiltered" "./results/ID_lists/test_sample_ids_donors.txt" "/home/nihtiju/work/NK_receptors/results/association_LD_newOnly/test/katalonia_test/association_katalonia_test" "./results/ID_lists/donors_test_relapse.txt"

# nc:
association_binary "/home/nihtiju/work/NK_receptors/results/extract_range/newcastle_plink_pgen_common_variants_filtered_ref_changed_newSNPsLDfiltered" "./results/ID_lists/test_sample_ids_donors.txt" "/home/nihtiju/work/NK_receptors/results/association_LD_newOnly/test/newcastle_test/association_newcastle_test" "./results/ID_lists/donors_test_relapse.txt"

# poland with phenotype file without cGvHD severe
association_poland_binary "/home/nihtiju/work/NK_receptors/results/extract_range/poland_plink_pgen_common_variants_filtered_ref_changed_newSNPsLDfiltered" "./results/ID_lists/test_sample_ids_donors.txt" "/home/nihtiju/work/NK_receptors/results/association_LD_newOnly/test/poland_test/association_poland_test" "./results/ID_lists/donors_test_relapse.txt"




