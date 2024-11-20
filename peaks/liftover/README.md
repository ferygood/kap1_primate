These files are for the liftOver comparison

*_overlap.bed: the final results after the process of liftOver and bedtools intersect
*_sorted.bed: lift finish file and sort by bedtools intersect
*_unmapped.bed: content that cannot be lifted
*.bed: such as oran_to_human.bed, oran lift to human file
*_TE_only_liftoverInput.bed: Inputs for liftover in that species
human_TE_only.bed: peaks that overlap with TEs in human
human_TE_only_sorted.bed: peaks that overlap with TEs and sorted by `bedtools sort -i`
