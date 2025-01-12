1.  `check_before_rmcorr.qmd` use PCA, and volcano plot to evaluate if we can use repeat measured correlation method.
2.  `rmcorr_human.qmd`, `rmcorr_chimp.qmd` and `rmcorr_oran.qmd`: normalize expression data to TPM and calculate repeated measured correlation in each species, human, chimp, oran. I additionally save the TPM files as `hm_expression_tpm.rds`, `chimp_expression_tpm.rds`, and `oran_expression_tpm.rds`
3.  `squire_to_correlation.qmd` this is the pilot calculation of if repeated measure correlation can use using squire output result.
4.  `squire.qmd` this file demonstrate how we use the SQUIRE package in command line.
5.  `tetranscripts_analysis.qmd` Do DE analysis in orthologous genes and transposable element across species. Output of DE genes and TEs tables are saved (DE_genes.csv, DE_TEs.csv).
6.  `tetranscripts_preprocess.qmd` The output of TEtranscript are raw counts. Thus, it is needed to first normalize and merge the results as an input for TEKRABber. The output is saved as `../../bcells_exp/tetranscriptDE.rds` including `hmchimpDE`, `hmoranDE`, `chimporanDE`, and `ortholog_hmchimp`.
7.  `filter_rmcorr.qmd` Use ChIP-seq detected results to filter if there are overlapped with rmcorr results across species.
8.  `check_repression.qmd` This script compare ChIP-seq and RNA-seq to see if binding TEs are not expressed in RNA-seq data in three species. There are 0, 26, and 3 overlapped in human, chimps, and orangutan respectively.
9.  `findOverlap_hmrmcorr_chipseq.qmd` This file indicate how to generate the 153 overlap pairs of TE:KRAB-ZNF network in human.
10. `tetranscripts_DE_visual.qmd` Visualize Differentially Expressed (DE) genes and TEs in heatmap. Output saved as jpeg images in "figures" folder.
11. `extract_motif_chipexo.qmd` There are 159 KRAB-ZNFs used in ChIP-exo, I want to extract their binding sequence.
12. `expression_profile_preparation.qmd` This file carry on the output from `rmcorr_human.qmd`, `rmcorr_chimp.qmd`, and `rmcorr_oran.qmd`. The output is `expression_merge.rds`.
13. `expression_profile_analysis.qmd` This file follows the previous file to keep the preparation step more clean and easy to maintain.
. 
