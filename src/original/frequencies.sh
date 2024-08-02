# calculate allele and genotype frequenies for the HSCT data

#---------------------------------------------------

# make plink1 files

plink2 \
	--pfile ./results/pgen_NK_SNPs/all_datasets_plink_pgen_common_variants_filtered_ref_changed \
	--make-bed \
	--out ./results/frequencies/plink1_all

plink2 \
	--pfile ./results/pgen_NK_SNPs/finns_plink_pgen_common_variants_filtered_ref_changed \
	--make-bed \
	--out ./results/frequencies/plink1_finland

plink2 \
	--pfile ./results/pgen_NK_SNPs/newcastle_plink_pgen_common_variants_filtered_ref_changed \
	--make-bed \
	--out ./results/frequencies/plink1_uk

plink2 \
	--pfile ./results/pgen_NK_SNPs/katalonia_plink_pgen_common_variants_filtered_ref_changed \
	--make-bed \
	--out ./results/frequencies/plink1_spain

plink2 \
	--pfile ./results/pgen_NK_SNPs/poland_plink_pgen_common_variants_filtered_ref_changed \
	--make-bed \
	--out ./results/frequencies/plink1_poland


#---------------------------------------------------
# discovery set

# combined populations
plink1.9 \
	--bfile ./results/frequencies/plink1_all \
	--keep ./results/ID_lists/train_sample_ids_donors.txt \
	--freq \
	--out ./results/frequencies/combined_allelefreq
plink1.9 \
	--bfile ./results/frequencies/plink1_all \
	--keep ./results/ID_lists/train_sample_ids_donors.txt \
	--freqx \
	--out ./results/frequencies/combined_genotypefreq

# finland
plink1.9 \
	--bfile ./results/frequencies/plink1_finland \
	--keep ./results/ID_lists/train_sample_ids_donors.txt \
	--freq \
	--out ./results/frequencies/finland_allelefreq
plink1.9 \
	--bfile ./results/frequencies/plink1_finland \
	--keep ./results/ID_lists/train_sample_ids_donors.txt \
	--freqx \
	--out ./results/frequencies/finland_genotypefreq

# uk
plink1.9 \
	--bfile ./results/frequencies/plink1_uk \
	--keep ./results/ID_lists/train_sample_ids_donors.txt \
	--freq \
	--out ./results/frequencies/uk_allelefreq
plink1.9 \
	--bfile ./results/frequencies/plink1_uk \
	--keep ./results/ID_lists/train_sample_ids_donors.txt \
	--freqx \
	--out ./results/frequencies/uk_genotypefreq

# spain
plink1.9 \
	--bfile ./results/frequencies/plink1_spain \
	--keep ./results/ID_lists/train_sample_ids_donors.txt \
	--freq \
	--out ./results/frequencies/spain_allelefreq
plink1.9 \
	--bfile ./results/frequencies/plink1_spain \
	--keep ./results/ID_lists/train_sample_ids_donors.txt \
	--freqx \
	--out ./results/frequencies/spain_genotypefreq


# poland
plink1.9 \
	--bfile ./results/frequencies/plink1_poland \
	--keep ./results/ID_lists/train_sample_ids_donors.txt \
	--freq \
	--out ./results/frequencies/poland_allelefreq
plink1.9 \
	--bfile ./results/frequencies/plink1_poland \
	--keep ./results/ID_lists/train_sample_ids_donors.txt \
	--freqx \
	--out ./results/frequencies/poland_genotypefreq

#---------------------------------------------------
# replication set

# combined populations
plink1.9 \
	--bfile ./results/frequencies/plink1_all \
	--keep ./results/ID_lists/test_sample_ids_donors.txt \
	--freq \
	--out ./results/frequencies/combined_allelefreq_replication
plink1.9 \
	--bfile ./results/frequencies/plink1_all \
	--keep ./results/ID_lists/test_sample_ids_donors.txt \
	--freqx \
	--out ./results/frequencies/combined_genotypefreq_replication

# finland
plink1.9 \
	--bfile ./results/frequencies/plink1_finland \
	--keep ./results/ID_lists/test_sample_ids_donors.txt \
	--freq \
	--out ./results/frequencies/finland_allelefreq_replication
plink1.9 \
	--bfile ./results/frequencies/plink1_finland \
	--keep ./results/ID_lists/test_sample_ids_donors.txt \
	--freqx \
	--out ./results/frequencies/finland_genotypefreq_replication

# uk
plink1.9 \
	--bfile ./results/frequencies/plink1_uk \
	--keep ./results/ID_lists/test_sample_ids_donors.txt \
	--freq \
	--out ./results/frequencies/uk_allelefreq_replication
plink1.9 \
	--bfile ./results/frequencies/plink1_uk \
	--keep ./results/ID_lists/test_sample_ids_donors.txt \
	--freqx \
	--out ./results/frequencies/uk_genotypefreq_replication

# spain
plink1.9 \
	--bfile ./results/frequencies/plink1_spain \
	--keep ./results/ID_lists/test_sample_ids_donors.txt \
	--freq \
	--out ./results/frequencies/spain_allelefreq_replication
plink1.9 \
	--bfile ./results/frequencies/plink1_spain \
	--keep ./results/ID_lists/test_sample_ids_donors.txt \
	--freqx \
	--out ./results/frequencies/spain_genotypefreq_replication


# poland
plink1.9 \
	--bfile ./results/frequencies/plink1_poland \
	--keep ./results/ID_lists/test_sample_ids_donors.txt \
	--freq \
	--out ./results/frequencies/poland_allelefreq_replication
plink1.9 \
	--bfile ./results/frequencies/plink1_poland \
	--keep ./results/ID_lists/test_sample_ids_donors.txt \
	--freqx \
	--out ./results/frequencies/poland_genotypefreq_replication

############################################################################################################

# in vitro blood donors

plink1.9 \
	--bfile ./results/leena_donors/leena_donors_chr1_extracted \
	--keep ./results/frequencies/in_vitro_keep.txt \
	--freq \
	--out ./results/frequencies/blooddonors_chr_1_allelefreq
plink1.9 \
	--bfile ./results/leena_donors/leena_donors_chr1_extracted \
	--keep ./results/frequencies/in_vitro_keep.txt \
	--freqx \
	--out ./results/frequencies/blooddonors_chr_1_genotypefreq

plink1.9 \
	--bfile ./results/leena_donors/leena_donors_chr18_extracted \
	--keep ./results/frequencies/in_vitro_keep.txt \
	--freq \
	--out ./results/frequencies/blooddonors_chr_18_allelefreq
plink1.9 \
	--bfile ./results/leena_donors/leena_donors_chr18_extracted \
	--keep ./results/frequencies/in_vitro_keep.txt \
	--freqx \
	--out ./results/frequencies/blooddonors_chr_18_genotypefreq


############################################################################################################
############################################################################################################

# frequencies with plink2

#---------------------------------------------------
# discovery set

# combined populations
plink2 \
	--pfile ./results/pgen_NK_SNPs/all_datasets_plink_pgen_common_variants_filtered_ref_changed \
	--keep ./results/ID_lists/train_sample_ids_donors.txt \
	--freq \
	--out ./results/frequencies/combined_allelefreq_plink2
plink2 \
	--pfile ./results/pgen_NK_SNPs/all_datasets_plink_pgen_common_variants_filtered_ref_changed \
	--keep ./results/ID_lists/train_sample_ids_donors.txt \
	--geno-counts \
	--out ./results/frequencies/combined_genotypefreq_plink2

# finland
plink2 \
	--pfile ./results/pgen_NK_SNPs/finns_plink_pgen_common_variants_filtered_ref_changed \
	--keep ./results/ID_lists/train_sample_ids_donors.txt \
	--freq \
	--out ./results/frequencies/finland_allelefreq_plink2
plink2 \
	--pfile ./results/pgen_NK_SNPs/finns_plink_pgen_common_variants_filtered_ref_changed \
	--keep ./results/ID_lists/train_sample_ids_donors.txt \
	--geno-counts \
	--out ./results/frequencies/finland_genotypefreq_plink2

# uk
plink2 \
	--pfile ./results/pgen_NK_SNPs/newcastle_plink_pgen_common_variants_filtered_ref_changed \
	--keep ./results/ID_lists/train_sample_ids_donors.txt \
	--freq \
	--out ./results/frequencies/uk_allelefreq_plink2
plink2 \
	--pfile ./results/pgen_NK_SNPs/newcastle_plink_pgen_common_variants_filtered_ref_changed \
	--keep ./results/ID_lists/train_sample_ids_donors.txt \
	--geno-counts \
	--out ./results/frequencies/uk_genotypefreq_plink2

# spain
plink2 \
	--pfile ./results/pgen_NK_SNPs/katalonia_plink_pgen_common_variants_filtered_ref_changed \
	--keep ./results/ID_lists/train_sample_ids_donors.txt \
	--freq \
	--out ./results/frequencies/spain_allelefreq_plink2
plink2 \
	--pfile ./results/pgen_NK_SNPs/katalonia_plink_pgen_common_variants_filtered_ref_changed \
	--keep ./results/ID_lists/train_sample_ids_donors.txt \
	--geno-counts \
	--out ./results/frequencies/spain_genotypefreq_plink2


# poland
plink2 \
	--pfile ./results/pgen_NK_SNPs/poland_plink_pgen_common_variants_filtered_ref_changed \
	--keep ./results/ID_lists/train_sample_ids_donors.txt \
	--freq \
	--out ./results/frequencies/poland_allelefreq_plink2
plink2 \
	--pfile ./results/pgen_NK_SNPs/poland_plink_pgen_common_variants_filtered_ref_changed \
	--keep ./results/ID_lists/train_sample_ids_donors.txt \
	--geno-counts \
	--out ./results/frequencies/poland_genotypefreq_plink2

#---------------------------------------------------
# replication set

# combined populations
plink2 \
	--pfile ./results/pgen_NK_SNPs/all_datasets_plink_pgen_common_variants_filtered_ref_changed \
	--keep ./results/ID_lists/test_sample_ids_donors.txt \
	--freq \
	--out ./results/frequencies/combined_allelefreq_replication_plink2
plink2 \
	--pfile ./results/pgen_NK_SNPs/all_datasets_plink_pgen_common_variants_filtered_ref_changed \
	--keep ./results/ID_lists/test_sample_ids_donors.txt \
	--geno-counts \
	--out ./results/frequencies/combined_genotypefreq_replication_plink2

# finland
plink2 \
	--pfile ./results/pgen_NK_SNPs/finns_plink_pgen_common_variants_filtered_ref_changed \
	--keep ./results/ID_lists/test_sample_ids_donors.txt \
	--freq \
	--out ./results/frequencies/finland_allelefreq_replication_plink2
plink2 \
	--pfile ./results/pgen_NK_SNPs/finns_plink_pgen_common_variants_filtered_ref_changed \
	--keep ./results/ID_lists/test_sample_ids_donors.txt \
	--geno-counts \
	--out ./results/frequencies/finland_genotypefreq_replication_plink2

# uk
plink2 \
	--pfile ./results/pgen_NK_SNPs/newcastle_plink_pgen_common_variants_filtered_ref_changed \
	--keep ./results/ID_lists/test_sample_ids_donors.txt \
	--freq \
	--out ./results/frequencies/uk_allelefreq_replication_plink2
plink2 \
	--pfile ./results/pgen_NK_SNPs/newcastle_plink_pgen_common_variants_filtered_ref_changed \
	--keep ./results/ID_lists/test_sample_ids_donors.txt \
	--geno-counts \
	--out ./results/frequencies/uk_genotypefreq_replication_plink2

# spain
plink2 \
	--pfile ./results/pgen_NK_SNPs/katalonia_plink_pgen_common_variants_filtered_ref_changed \
	--keep ./results/ID_lists/test_sample_ids_donors.txt \
	--freq \
	--out ./results/frequencies/spain_allelefreq_replication_plink2
plink2 \
	--pfile ./results/pgen_NK_SNPs/katalonia_plink_pgen_common_variants_filtered_ref_changed \
	--keep ./results/ID_lists/test_sample_ids_donors.txt \
	--geno-counts \
	--out ./results/frequencies/spain_genotypefreq_replication_plink2


# poland
plink2 \
	--pfile ./results/pgen_NK_SNPs/poland_plink_pgen_common_variants_filtered_ref_changed \
	--keep ./results/ID_lists/test_sample_ids_donors.txt \
	--freq \
	--out ./results/frequencies/poland_allelefreq_replication_plink2
plink2 \
	--pfile ./results/pgen_NK_SNPs/poland_plink_pgen_common_variants_filtered_ref_changed \
	--keep ./results/ID_lists/test_sample_ids_donors.txt \
	--geno-counts \
	--out ./results/frequencies/poland_genotypefreq_replication_plink2




















