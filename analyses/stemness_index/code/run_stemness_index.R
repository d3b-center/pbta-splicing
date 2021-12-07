# Author: Komal S. Rathi
# Function: Run OCLR based Stemness profiling
# References: 
# Paper: https://www.cell.com/cell/pdf/S0092-8674(18)30358-1.pdf
# Tutorial: http://tcgabiolinks.fmrp.usp.br/PanCanStem/mRNAsi.html

# load libraries
# use synapser instead of synapseClient
suppressPackageStartupMessages({
  library(synapser) 
  library(gelnet)
  library(ggplot2)
})

# Synapse login required to get PCBC data
readRenviron('~/.Renviron')
SYNAPSE_ID <- Sys.getenv("SYNAPSE_ID")
SYNAPSE_PWD <- Sys.getenv("SYNAPSE_PWD")
synLogin(email = SYNAPSE_ID, password = SYNAPSE_PWD) # your id and password from https://www.synapse.org/

# root directory
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
analyses_dir <- file.path(root_dir, 'analyses', 'stemness_index')
output_dir <- file.path(analyses_dir, 'output')

# source scripts for training and prediction
source(file.path(analyses_dir, 'util', 'genes2hugo.R')) # convert entrez ids to hugo symbols
source(file.path(analyses_dir, 'util', 'main.train.R')) # train using PCBC dataset
source(file.path(analyses_dir, 'util', 'main.predict.R')) # predict using vector generated in main.train.R
source(file.path(analyses_dir, 'util', 'stemness_index.R')) # script to calculate and plot stemness scores

# run training set (the output of this should be vector of length 78 containing only 1s)
# stem/progenitor cells from the Progenitor Cell Biology Consortium
fname <- file.path(analyses_dir, 'input', 'pcbc-stemsig.tsv')
if(file.exists(fname)){
  print("Training set ready")
} else {
  print("Run training set")
  auc <- main.train(fnOut = fname, fnGenes = NULL)
}

# list (n = 8)
p <- stemness_index(expr = file.path(root_dir, 'data', 'tpm_matrix.rds'), 
                    meta = file.path(root_dir, 'data', 'metadata.rds'), 
                    fnSig = fname,
                    fnOut = file.path(output_dir, 'stemness_scores.tsv'), 
                    method = "t.test")
ggsave(filename = file.path(output_dir, 'stemness_index_plot.pdf'), device = "pdf", plot = p,width = 10, height = 7)
