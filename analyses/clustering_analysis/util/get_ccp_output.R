# function to get ccp clustering output for a specific combination of distance + algorithm + % variable genes 
suppressPackageStartupMessages({
  library(ConsensusClusterPlus)
  library(DGCA)
  library(DESeq2)
})

get_ccp_output <- function(input_mat, data_type = c("raw_counts", "non_expr"), var_genes, cluster_algorithm, cluster_distance, prefix){
  
  # set seed for reproducibility
  set.seed(100)
  
  # remove rows where all values are 0
  input_mat <- input_mat[rowSums(input_mat[])>0,]
  
  # remove columns and rows with sd == 0 
  input_mat <- input_mat[apply(input_mat, MARGIN = 1,function(x) { sd(x) != 0} ),] 
  input_mat <- input_mat[,apply(input_mat, MARGIN = 2,function(x) { sd(x) != 0} )] 
  
  # filter features for expression data
  if(data_type == "raw_counts"){
    
    # filter by expression if var_prop > 0
    if(var_genes > 0){
      print("filter by expression")
      input_mat <- DGCA::filterGenes(inputMat = input_mat, 
                                     filterTypes = c("central", "dispersion"),
                                     filterDispersionType = "cv", 
                                     filterDispersionPercentile = 0.2,
                                     sequential= TRUE)
    }
    
    # transform raw counts using vst
    print("transform using VST")
    input_mat <- DESeq2::varianceStabilizingTransformation(round(as.matrix(input_mat)), blind = TRUE, fitType = "parametric")
    input_mat <- as.data.frame(input_mat)
  }
  
  # perform dip test if var_prop == 0 
  if(var_genes == 0){
    print("perform dip test")
    input_mat <- perform_diptest(count_matrix = input_mat)
  } else if(var_genes > 0){
    # filter for variability if var_prop > 0
    print("filter by % variability")
    var_genes <- as.numeric(var_genes)
    print(paste("subset to", var_genes, "%", "variable features"))
    n <- round(var_genes/100*nrow(input_mat))
    
    # reduce the dataset to the most variable features, measured by median absolute deviation
    mads <- apply(input_mat, 1, mad)
    input_mat <- input_mat[rev(order(mads))[1:n],] 
  }
  
  # output folder 
  ccp_dir <- file.path('output/ccp_output')
  dir.create(ccp_dir, showWarnings = F, recursive = T)
  
  # output prefix
  prefix <- paste(prefix, cluster_algorithm, cluster_distance, var_genes, sep = "_")
  
  # run ConsensusClusterPlus
  res.ccp <- ConsensusClusterPlus(d = as.matrix(input_mat),
                                  clusterAlg = cluster_algorithm,
                                  finalLinkage = "average",
                                  distance = cluster_distance,
                                  plot = "pdf",
                                  reps = 50, 
                                  maxK = 10, 
                                  pItem = 0.8,
                                  title = file.path(ccp_dir),
                                  seed = 100)
  
  # rename consensus pdf output
  old_name <- file.path(ccp_dir, 'consensus.pdf')
  new_name <- file.path(ccp_dir, paste0(prefix, '_consensus.pdf'))
  file.rename(from = old_name, to = new_name)
  
  # save entire consensus output 
  saveRDS(res.ccp, file = file.path(ccp_dir, paste0(prefix, '_ccp.rds')))
  
  # save input matrix (filtered/normalized) for downstream analyses
  saveRDS(input_mat, file = file.path(ccp_dir, paste0(prefix, '_matrix.rds')))
}