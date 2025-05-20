1. ReadCountTable_DESeq2_prefiltered_allSpecies: Filtered read count table (sum >= 10). 
   NHP were mapped to their corresponding ref genome and liftovered

2. GENCODE_incUnknown_genecounts_HSref: Very first read count table after feature count calling. 
   NHP were also mapped against the human ref genome fyi.

3. ids_names_extractedFromGTF: I called featurecount with gene id + version, 
   so I additionally made a list of gene id: gene symbol, in case you're also interested. 
   Gene symbols are extracted from the exact gtf file version used for counting
