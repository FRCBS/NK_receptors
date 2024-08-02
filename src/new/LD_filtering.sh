# LD-prune variants

# find LD-based variants
plink2 \
	--pfile /home/nihtiju/work/NK_receptors/results/extract_range/all_datasets_plink_pgen_common_variants_filtered_ref_changed \
	--indep-pairwise 500kb 0.2 \
	--out /home/nihtiju/work/NK_receptors/results/extract_range/LD_pruning

# exracting LD variants
plink2 \
	--pfile /home/nihtiju/work/NK_receptors/results/extract_range/all_datasets_plink_pgen_common_variants_filtered_ref_changed \
	--extract /home/nihtiju/work/NK_receptors/results/extract_range/LD_pruning.prune.in \
	--make-pgen \
	--out /home/nihtiju/work/NK_receptors/results/extract_range/all_datasets_plink_pgen_common_variants_filtered_ref_changed_LD_pruned




# --indep-pairwise 100kb 0.8 --> 2390 variants remaining
# --indep-pairwise 200kb 0.5 --> 1510 variants remaining
# --indep-pairwise 500kb 0.2 --> 851 variants remaining

#-----------------------------------------------------------------------

# for finns
# exracting LD variants
plink2 \
	--pfile /home/nihtiju/work/NK_receptors/results/extract_range/finns_plink_pgen_common_variants_filtered_ref_changed \
	--extract /home/nihtiju/work/NK_receptors/results/extract_range/LD_pruning.prune.in \
	--make-pgen \
	--out /home/nihtiju/work/NK_receptors/results/extract_range/finns_plink_pgen_common_variants_filtered_ref_changed_LD_pruned

#-----------------------------------------------------------------------

# for uk
# exracting LD variants
plink2 \
	--pfile /home/nihtiju/work/NK_receptors/results/extract_range/newcastle_plink_pgen_common_variants_filtered_ref_changed \
	--extract /home/nihtiju/work/NK_receptors/results/extract_range/LD_pruning.prune.in \
	--make-pgen \
	--out /home/nihtiju/work/NK_receptors/results/extract_range/newcastle_plink_pgen_common_variants_filtered_ref_changed_LD_pruned

#-----------------------------------------------------------------------

# for spain
# exracting LD variants
plink2 \
	--pfile /home/nihtiju/work/NK_receptors/results/extract_range/katalonia_plink_pgen_common_variants_filtered_ref_changed \
	--extract /home/nihtiju/work/NK_receptors/results/extract_range/LD_pruning.prune.in \
	--make-pgen \
	--out /home/nihtiju/work/NK_receptors/results/extract_range/katalonia_plink_pgen_common_variants_filtered_ref_changed_LD_pruned

#-----------------------------------------------------------------------

# for poland
# exracting LD variants
plink2 \
	--pfile /home/nihtiju/work/NK_receptors/results/extract_range/poland_plink_pgen_common_variants_filtered_ref_changed \
	--extract /home/nihtiju/work/NK_receptors/results/extract_range/LD_pruning.prune.in \
	--make-pgen \
	--out /home/nihtiju/work/NK_receptors/results/extract_range/poland_plink_pgen_common_variants_filtered_ref_changed_LD_pruned



