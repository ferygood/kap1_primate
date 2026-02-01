library(dplyr)
library(data.table)
library(GenomicRanges)

#' Analyze Intersection between Genes, TEs, and Peaks
#' @param peak_data Dataframe containing peak calls (chr, start, end, peak_id)
#' @param gene_list Dataframe containing target Ensembl IDs and Symbols (ENSG, geneName)
#' @param gene_ref Dataframe of the genome annotation (GTF/TSV with #chrom, chromStart, chromEnd, etc.)
#' @param te_ref Dataframe of the RepeatMasker annotation (genoName, genoStart, genoEnd, etc.)
#' @return A filtered dataframe with gene-TE-peak intersections and strand relations
analyze_gene_te_peak_intersection <- function(peak_data, gene_list, gene_ref, te_ref) {
  
  # 1. Filter Gene Reference by target gene list
  # Using backticks for columns starting with # to avoid syntax errors
  message("Filtering gene reference...")
  gene_ref_filtered <- gene_ref %>%
    filter(ensembl_gene_id %in% gene_list$ENSG)
  
  # 2. Filter TE Reference by major classes
  message("Filtering TE reference...")
  te_ref_filtered <- te_ref %>%
    filter(repClass %in% c("SINE", "LINE", "LTR", "DNA"))
  
  # 3. Create GRanges objects for Genes and TEs
  gene_gr <- GRanges(
    seqnames = gene_ref_filtered$`#chrom`,
    ranges = IRanges(start = gene_ref_filtered$chromStart, end = gene_ref_filtered$chromEnd),
    strand = gene_ref_filtered$strand,
    ensembl_id = gene_ref_filtered$ensembl_gene_id
  )
  
  te_gr <- GRanges(
    seqnames = te_ref_filtered$genoName, 
    ranges = IRanges(start = te_ref_filtered$genoStart, end = te_ref_filtered$genoEnd),
    strand = te_ref_filtered$strand,
    repName = te_ref_filtered$repName
  )
  
  # 4. Find Gene-TE Overlaps (Ignoring strand for initial discovery)
  message("Finding Gene-TE overlaps...")
  hits_gt <- findOverlaps(gene_gr, te_gr, ignore.strand = TRUE)
  
  gt_output <- data.frame(
    chr_name    = as.character(seqnames(gene_gr)[queryHits(hits_gt)]),
    gene_id     = mcols(gene_gr)$ensembl_id[queryHits(hits_gt)],
    gene_start  = start(gene_gr)[queryHits(hits_gt)],
    gene_end    = end(gene_gr)[queryHits(hits_gt)],
    gene_strand = as.character(strand(gene_gr)[queryHits(hits_gt)]),
    te_name     = mcols(te_gr)$repName[subjectHits(hits_gt)],
    te_start    = start(te_gr)[subjectHits(hits_gt)],
    te_end      = end(te_gr)[subjectHits(hits_gt)],
    te_strand   = as.character(strand(te_gr)[subjectHits(hits_gt)])
  )
  
  # 5. Create GRanges for Peaks and intersect with the Gene-TE result
  message("Filtering by Peak regions...")
  peak_gr <- GRanges(
    seqnames = peak_data$chr,
    ranges   = IRanges(start = peak_data$start, end = peak_data$end),
    peak_id  = peak_data$peak_id
  )
  
  # Use TE coordinates from the Gene-TE result to check Peak overlap
  gt_res_gr <- GRanges(
    seqnames = gt_output$chr_name,
    ranges   = IRanges(start = gt_output$te_start, end = gt_output$te_end)
  )
  
  hits_peak <- findOverlaps(gt_res_gr, peak_gr, ignore.strand = TRUE)
  
  # Keep only Gene-TE pairs that overlap with a Peak
  final_res <- gt_output[unique(queryHits(hits_peak)), ]
  
  # 6. Merge with gene symbols and add Strand Relation metadata
  message("Merging symbols and calculating strand relations...")
  final_res <- merge(
    final_res, 
    gene_list[, c("ENSG", "geneName")], 
    by.x = "gene_id", by.y = "ENSG", all.x = TRUE
  )
  
  # Calculate Sense/Antisense relation
  final_res$strand_relation <- ifelse(
    final_res$gene_strand == final_res$te_strand, 
    "Sense", "Antisense"
  )
  
  # Reorganize columns for clarity
  final_res <- final_res %>%
    select(chr_name, ensembl_id = gene_id, gene_symbol = geneName, 
           gene_start, gene_end, gene_strand, 
           te_name, te_start, te_end, te_strand, strand_relation)
  
  message("Analysis Complete.")
  return(list(data = final_res, hits = hits_peak, peak_overlap = final_res))
}