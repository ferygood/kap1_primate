1. `ChIPseq_findGene.qmd` The file use narrow peak results to find nearby genes.
2. `analysis.qmd` analyze peak count data, generating upset_TE.jpg, peak count barplot (fig2b.jpg), and save the interacted 589 TEs as csv file 
(tables/TE_589_intersect.csv).
3. `GO_genes_DE.qmd` Merge the information of GO term related to neuron development and differentially expressed information into a dataframe. Create a 
function `merge_de()` and save the output as a dataframe, `tables/GO_withDEinfo.csv`. Then the visualizations are created, including 
`upset_chimp_oran_GO_DE.jpg` and `GO_DE_barplot.jpg`.
