# Function: run diptest and filter matrix based on diptest results
suppressPackageStartupMessages({
  library(diptest)
  library(tidyverse)
})

perform_diptest <- function(count_matrix, pval_cutoff = 0.05){
  
  count_matrix_normalized <- count_matrix
  
  # define a dataframe to store the pval
  count_matrix_pval <- data.frame()
  count_matrix_normalized <- as.data.frame(count_matrix_normalized)
  for(i in 1:nrow(count_matrix_normalized)){
    # calculate dip.test for each gene
    diptest_each <- count_matrix_normalized[i,] %>% 
      as.numeric() %>%
      diptest::dip.test(simulate.p.value = FALSE, B = 2000)
    # gather the pvalue for each gene
    pval_each <- diptest_each$p.value
    # add another column to store pval
    count_matrix_pval_each <- count_matrix_normalized[i,] %>%
      as.data.frame() %>%
      dplyr::mutate(pval = pval_each)
    # cocgine each line back to the dataframe 
    count_matrix_pval <- bind_rows(count_matrix_pval,  count_matrix_pval_each)
  }
  
  # define rownames
  rownames(count_matrix_pval) <- rownames(count_matrix_normalized)
  
  # Filter to expression that only has pval<pval_cutoff if over min_n genes satisfy this or take top min_n
  # see how many genes have pval<pval_cutoff
  n <- count_matrix_pval %>% 
    dplyr::filter(pval < pval_cutoff) %>%
    nrow()
  
  print(n)
  
  count_matrix_pval <- count_matrix_pval %>% 
    dplyr::filter(pval < pval_cutoff) %>% 
    dplyr::select(-pval)
  return(count_matrix_pval)
}
