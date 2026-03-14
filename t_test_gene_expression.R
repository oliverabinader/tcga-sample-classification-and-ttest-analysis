# ============================================
# Sample classification + gene-wise t-tests
# ============================================


# -------------------------------------------------------
# Assign tumor/normal subtype from sample IDs (TCGA data)
# -------------------------------------------------------
get_subtype <- function(sample_id) {
  if (endsWith(sample_id, "01") || endsWith(sample_id, "02") || endsWith(sample_id, "06")) {
    return("cancer")
  } else if (endsWith(sample_id, "11")) {
    return("normal")
  } else {
    return("unknown")
  }
}

# Example usage:
# cancer_type$Subtype <- sapply(strsplit(cancer_type$sample, "-"), function(x) get_subtype(tail(x, n = 1)))
# cancer_type must hold metadat about the samples. The dataset should have a "sample" column stating all samples in the dataset
# and it will create the Subtype column and assign type for each sample.


# --------------------------------------------
# Perform gene-wise t-test (cancer vs normal)
# --------------------------------------------
run_gene_ttest <- function(genes_of_interest, cancer_data, normal_data) {

  # Pre-allocate vectors to store gene names and corresponding p-values
  gene_names <- vector("character", length(genes_of_interest)) 
  p_values <- numeric(length(genes_of_interest))

  # Loop through each gene in the list of genes of interest
  for (i in seq_along(genes_of_interest)) {
    gene <- genes_of_interest[i] # Extract the current gene name

    # Extract expression values for this gene from the cancer dataset. Make sure no gene_symbol is present in that dataset.
    cancer_expression <- as.numeric(cancer_data[cancer_data$gene_symbol == gene, ])

    # Do the same extarction for the normal dataset.
    normal_expression <- as.numeric(normal_data[normal_data$gene_symbol == gene, ])

    # Perform Welch two-sample t-test (does not assume equal variance)
    t_test_result <- t.test(cancer_expression, normal_expression, alternative = "two.sided", paired = FALSE, var.equal = FALSE)
    # tests for any difference in means
    # samples are independent
    # Welch correction for unequal variances

    # Store the gene name and the resulting p-value
    gene_names[i] <- gene
    p_values[i] <- t_test_result$p.value
  }

  # Combine results into a dataframe
  result_df <- data.frame(gene_name = gene_names, p_value = p_values, stringsAsFactors = FALSE) %>%
    arrange(gene_name)
  
  return(result_df)
}
