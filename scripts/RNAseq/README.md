1.  `check_before_rmcorr.qmd` use PCA, and volcano plot to evaluate if we can use repeat measured correlation method.
2.   `rmcorr_human.qmd`, `rmcorr_chimp.qmd` and `rmcorr_oran.qmd`: normalize expression data to TPM and calculate repeated measured correlation in each species, human, chimp, oran.
3.  `squire_to_correlation.qmd` this is the pilot calculation of if repeated measure correlation can use using squire output result.
4.  `squire.qmd` this file demonstrate how we use the SQUIRE package in command line.
5.  `tetranscripts_analysis.qmd` Do DE analysis in orthologous genes and transposable element across species. Output of DE genes and TEs tables are saved (DE_genes.csv, DE_TEs.csv).
6.  `tetranscripts_preprocess.qmd` The output of TEtranscript are raw counts. Thus, it is needed to first normalize and merge the results as an input for TEKRABber. The output is saved as `../../bcells_exp/tetranscriptDE.rds` including `hmchimpDE`, `hmoranDE`, `chimporanDE`, and `ortholog_hmchimp`.
7.  
