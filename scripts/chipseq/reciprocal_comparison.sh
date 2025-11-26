#!/bin/bash

# ------------------------------
# Paths to BED files
# ------------------------------
HUMAN=human_TE_only_6col.bed
CHIMP=chimp_TE_only_6col.bed
ORAN=oran_TE_only_6col.bed

# ------------------------------
# Chain files
# ------------------------------
HG2CH=hg19ToPanTro4.over.chain
CH2HG=panTro4ToHg19.over.chain
HG2OR=hg19ToPonAbe2.over.chain
OR2HG=ponAbe2ToHg19.over.chain

# ------------------------------
# 1. Reciprocal liftover: human <-> chimp
# ------------------------------
liftOver $HUMAN $HG2CH human_to_chimp.bed unmapped_human_to_chimp.bed
liftOver human_to_chimp.bed $CH2HG human_recovered_from_chimp.bed unmapped_chimp_to_human.bed
bedtools intersect -u -a $HUMAN -b human_recovered_from_chimp.bed > human_chimp_reciprocal.bed

# ------------------------------
# 2. Reciprocal liftover: human <-> orangutan
# ------------------------------
liftOver $HUMAN $HG2OR human_to_oran.bed unmapped_human_to_oran.bed
liftOver human_to_oran.bed $OR2HG human_recovered_from_oran.bed unmapped_oran_to_human.bed
bedtools intersect -u -a $HUMAN -b human_recovered_from_oran.bed > human_oran_reciprocal.bed

# ------------------------------
# 3. Pairwise intersections
# ------------------------------
# human ∩ chimp
bedtools intersect -u -a human_chimp_reciprocal.bed -b $CHIMP > human_chimp_overlap.bed
# human ∩ orang
bedtools intersect -u -a human_oran_reciprocal.bed -b $ORAN > human_oran_overlap.bed
# chimp ∩ orang (reciprocal chimp ↔ orang)
liftOver $CHIMP hg19ToPonAbe2.over.chain.gz chimp_to_oran.bed unmapped_chimp_to_oran.bed
liftOver chimp_to_oran.bed ponAbe2ToHg19.over.chain.gz chimp_recovered_from_oran.bed unmapped_oran_to_chimp.bed
bedtools intersect -u -a $CHIMP -b chimp_recovered_from_oran.bed | bedtools intersect -u -a - -b $ORAN > chimp_oran_overlap.bed

# ------------------------------
# 4. Three-way conserved
# ------------------------------
bedtools intersect -u -a human_chimp_overlap.bed -b human_oran_overlap.bed > human_chimp_oran_conserved.bed

# ------------------------------
# 5. Species-specific peaks
# ------------------------------
bedtools subtract -a $HUMAN -b human_chimp_oran_conserved.bed > human_only.bed
bedtools subtract -a $CHIMP -b human_chimp_oran_conserved.bed > chimp_only.bed
bedtools subtract -a $ORAN -b human_chimp_oran_conserved.bed > oran_only.bed

# ------------------------------
# 6. Count numbers for table
# ------------------------------
H_ONLY=$(wc -l < human_only.bed)
C_ONLY=$(wc -l < chimp_only.bed)
O_ONLY=$(wc -l < oran_only.bed)

H_C=$(wc -l < human_chimp_overlap.bed)
H_O=$(wc -l < human_oran_overlap.bed)
C_O=$(wc -l < chimp_oran_overlap.bed)
H_C_O=$(wc -l < human_chimp_oran_conserved.bed)

# ------------------------------
# 7. Output table
# ------------------------------
echo -e "\tH\tC\tO"
echo -e "H\t$H_ONLY\t$H_C\t$H_O"
echo -e "C\t\t$C_ONLY\t$C_O"
echo -e "O\t\t\t$O_ONLY"
echo
echo "Three-way conserved (H∩C∩O): $H_C_O"

# ------------------------------
# 8. Move all generated BED files into reciprocal_bed folder
# ------------------------------
mkdir -p reciprocal_bed
mv human_to_chimp.bed human_recovered_from_chimp.bed human_chimp_reciprocal.bed \
   human_to_oran.bed human_recovered_from_oran.bed human_oran_reciprocal.bed \
   human_chimp_overlap.bed human_oran_overlap.bed chimp_to_oran.bed \
   chimp_recovered_from_oran.bed chimp_oran_overlap.bed human_chimp_oran_conserved.bed \
   human_only.bed chimp_only.bed oran_only.bed \
   unmapped_human_to_chimp.bed \
   unmapped_human_to_oran.bed \
   unmapped_chimp_to_human.bed \
   unmapped_chimp_to_oran.bed \
   unmapped_oran_to_human.bed \
   unmapped_oran_to_chimp.bed \
   reciprocal_bed/
