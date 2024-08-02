# extract the NK receptor SNPs from all populations

# for eur samples
# leaving in only SNPs which are common to all datasets -> ok to extract only from one dataset here & then filter with the list of common SNPs made elsewhere

# form ic1: after imputation all chrs merged and turned into bed/bim/fam, not info-filtered

plink2 \
	--bfile /media/volume/datavolume/imputation/ic1/ic1_b38_imputed_info_all_chrs \
	--extract range /media/volume/datavolume/NK_cell_markers_and_relapse/ranges/range_file.txt \
	--make-bed \
	--out /media/volume/datavolume/NK_cell_markers_and_relapse/ranges/extracted_ranges_ic1











