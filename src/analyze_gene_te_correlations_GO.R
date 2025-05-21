analyze_gene_te_correlations <- function(df, df_TE) {
    # Filter rows with non-empty TE_idx
    df <- df[!is.na(df$TE_idx) & df$TE_idx != "", ]

    # Extract TE expression matrix
    te_expr <- df_TE[, 2:10]
    rownames(te_expr) <- rownames(df_TE)

    # Initialize result list
    cor_results <- list()

    for (i in 1:nrow(df)) {
        gene_id <- df$ENSID[i]
        gene_exp <- as.numeric(df[i, 5:13])
        te_ids_str <- df$TE_idx[i]

        if (!is.na(te_ids_str) && te_ids_str != "") {
            te_ids <- strsplit(te_ids_str, ",")[[1]]

            # Get TE expression for these TEs, reshape to samples x TEs
            te_matrix <- t(te_expr[te_ids, , drop = FALSE])
            te_matrix <- as.data.frame(te_matrix)
            te_matrix$gene <- gene_exp

            # Fit multiple linear regression model
            fit <- lm(gene ~ ., data = te_matrix)
            summary_fit <- summary(fit)

            # Extract coefficients and p-values (excluding intercept)
            coefs <- summary_fit$coefficients[-1, "Estimate"]
            pvals <- summary_fit$coefficients[-1, "Pr(>|t|)"]
            fdrs  <- p.adjust(pvals, method = "fdr")

            # Save result
            cor_results[[length(cor_results) + 1]] <- data.frame(
                gene_id = gene_id,
                te_ids = paste(te_ids, collapse = ","),
                te_coef = paste(round(coefs, 4), collapse = ","),
                te_pval = paste(signif(pvals, 3), collapse = ","),
                te_fdr  = paste(signif(fdrs, 3), collapse = ","),
                adj_r_squared = summary_fit$adj.r.squared
            )
        }
    }

    # Combine results into a single data frame
    cor_df <- do.call(rbind, cor_results)
    return(cor_df)
}
