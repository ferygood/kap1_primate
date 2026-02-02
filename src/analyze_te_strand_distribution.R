library(ggplot2)
library(scales)
library(tidyr)
library(dplyr)

#' Enhanced TE Strand Analysis for Class or Family
#' @param result_data The $data from intersection function
#' @param te_ref TE reference table
#' @param group_by_var Column to group by: "repClass" or "repFamily"
#' @param min_count Minimum total occurrences to include in the plot/test
analyze_te_strand_distribution <- function(result_data, te_ref, species_name = "Species", 
                                           group_by_var = "repClass", min_count = 1) {
  
  # Format display name for titles
  display_name <- if(group_by_var == "repClass") "TE Class" else "TE Family"
  
  message(paste("Processing", display_name, "statistics for", species_name, "..."))
  
  # 1. Prepare and Merge Metadata (Filtered to major TE classes)
  te_info <- unique(te_ref[, .(repName, repClass, repFamily)]) %>%
    filter(repClass %in% c("SINE", "LINE", "LTR", "DNA")) # Filter step added here
  
  # 2. Merge info back to result
  extended_data <- merge(result_data, te_info, by.x = "te_name", by.y = "repName", all.x = TRUE)
  
  # Remove rows that were filtered out from te_info
  extended_data <- extended_data %>% filter(!is.na(repClass))
  
  # 3. Calculate Distribution with dynamic grouping
  dist_stats <- extended_data %>%
    dplyr::rename(group_col = !!rlang::sym(group_by_var)) %>% 
    dplyr::group_by(group_col, strand_relation) %>%
    dplyr::summarise(count = n(), .groups = 'drop') %>%
    dplyr::group_by(group_col) %>%
    dplyr::mutate(total = sum(count),
                  percentage = (count / total)) %>%
    dplyr::filter(total >= min_count)
  
  if(nrow(dist_stats) == 0) {
    stop("No groups met the minimum count threshold after filtering.")
  }

  # 4. Statistical Test (Chi-squared on filtered data)
  test_data <- extended_data %>% filter(!!rlang::sym(group_by_var) %in% dist_stats$group_col)
  contingency_table <- table(test_data[[group_by_var]], test_data$strand_relation)
  chisq_test <- chisq.test(contingency_table)
  
  # 5. Visualization
  p <- ggplot(dist_stats, aes(x = reorder(group_col, -total), y = percentage, fill = strand_relation)) +
    geom_bar(stat = "identity", position = "fill", color = "white", width = 0.7) +
    scale_y_continuous(labels = percent_format()) +
    scale_fill_manual(values = c("Antisense" = "#E41A1C", "Sense" = "#377EB8")) +
    theme_minimal() +
    labs(title = paste(species_name, "-", display_name, "Strand Detection"),
         subtitle = paste("Chi-squared test p-value:", format.pval(chisq_test$p.value)),
         x = display_name, y = "Percentage", fill = "Relation") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(face = "bold"))
  
  return(list(extended_data = extended_data, stats = dist_stats, p_value = chisq_test$p.value, plot = p))
}

#' Updated Preference Test to support Class or Family
analyze_te_strand_preference <- function(extended_data, group_by_var = "repClass") {
  
  # Filter only the relevant classes before analysis
  filtered_data <- extended_data %>% 
    filter(repClass %in% c("SINE", "LINE", "LTR", "DNA"))

  summary_table <- filtered_data %>%
    group_by(!!rlang::sym(group_by_var), strand_relation) %>%
    summarise(count = n(), .groups = 'drop') %>%
    pivot_wider(names_from = strand_relation, values_from = count, values_fill = 0)
  
  # Ensure both columns exist for calculation
  if(!"Antisense" %in% colnames(summary_table)) summary_table$Antisense <- 0
  if(!"Sense" %in% colnames(summary_table)) summary_table$Sense <- 0

  summary_table <- summary_table %>%
    rowwise() %>%
    mutate(
      total = Sense + Antisense,
      antisense_prop = Antisense / total,
      p_value = binom.test(Antisense, total, p = 0.5)$p.value,
      sig = ifelse(p_value < 0.05, "*", "ns"),
      preference = ifelse(antisense_prop > 0.5, "Antisense Bias", 
                          ifelse(antisense_prop < 0.5, "Sense Bias", "No Bias"))
    ) %>%
    ungroup() %>%
    arrange(p_value)
  
  return(summary_table)
}