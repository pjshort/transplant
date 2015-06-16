# script to correlate gene expression with clinical features via gene modules

# relies on the WGCNA package as well as the user pre-determining a soft thresholding
# parameter before running the analysis

# saves gene module to phenotype correlations and module assignments as RData file

### dependencies
library(optparse)
library(stringr)
library(oligo)
library(WGCNA)
source("../scripts/eset_tools.R")
source("../scripts/cluster_tools.R")g
source("../scripts/correlation_tools.R")
options(stringsAsFactors = FALSE)
allowWGCNAThreads()

library(yaml)
config = yaml.load_file("../config.yml")
creds <- list(username=config$credentials$username, password=config$credentials$pass)

### command line options
option_list <- list(
  make_option("--search_path", default="../data/TransplantCELs/discovery",
              help="Top level directory to start recursive search for *.CEL files."),
  make_option("--phenotype_continuous", default = NULL, help = "csv file with continuous 
              phenotypic outcomes for each patient."),
  make_option("--phenotype_discrete", default = NULL, help = "csv file with discrete 
              phenotypic outcomes for each patient."),
  make_option("--out", default="../results/discovery_results.RData",
              help="Location to save the ExpressionSet object. Defaults to ./eset.Rdata"),
  make_option("--soft_threshold_power", default = 1, help = "Set soft thresholding parameter 
              such that the resulting network approximates a scale free distribution."),
  make_option("--fdr", default = 0.05, help = "Set false discovery rate tolerated for 
              declaring a result significant (i.e. what proportion of significant results
              will be expected to be fall rejections of the null hypothesis)."),
  make_option("--merge", default = FALSE, help = "If TRUE, gene modules will be merged
              if eigengene expression levels are >80% correlated."),
  make_option("--verbose", action="store_true", default=FALSE,
              help="Print extra output advising the user of progression through the analysis.")
)
args <- parse_args(OptionParser(option_list=option_list))

args$search_path <- "../data/TransplantCELs/discovery"
args$phenotype_continuous <- "../data/transplant_continuous_phenotypes.csv"
args$phenotype_discrete <- "../data/transplant_discrete_phenotypes.csv"
args$soft_threshold_power <- 10

####     create expression set     ####

# retrieve the CEL files
CEL_files = list.files(path = args$search_path, pattern = ".*CEL", all.files = FALSE, full.names = TRUE, recursive = TRUE, ignore.case = TRUE)

# normalize the probe intensity
eset = rma_eset(CEL_list = CEL_files)

# process eset to keep only one probe per gene, only annotated genes
affy_annotation = "../data/HuGene-1_0-st-v1.na35.hg19.transcript.csv"
# (download this file from http://www.affymetrix.com/support/technical/annotationfilesmain.affx)

eset <- gene_filter_eset(eset, affy_annotation)

# extract the expression data and rotate so we have samples (rows) by genes (columns)
expr_data <- t(eset@assayData$exprs)

# change sample names to Tritan IDs
rownames(expr_data) <- str_extract(rownames(expr_data), "[6][0-9][0-9]") # matches 601 - 699
colnames(eset) <- rownames(expr_data)

n_genes = ncol(expr_data)
n_samples = nrow(expr_data)



####     build the network     ####

#TODO rewrite this section to calculate the tree from scratch

# soft_thresh_power = 10

# note - this takes a long time and has been pre-calculated and saved
#clust_out <- build_tree(expr_data, soft_thresh_power)
#gene_tree = clust_out[[1]]
#dissTOM = clust_out[[2]]

# module identification using dynamic tree cut
#modules = cutreeDynamic(dendro = gene_tree, distM = dissTOM,
#deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = 30)
#module_colors = labels2colors(modules)

#save(gene_tree, module_colors, file = "../data/geneTree.RData")
load("../data/geneTree.RData")

n_modules = length(unique(module_colors))

# define eigengenes with expression value as first PC
MEs = moduleEigengenes(expr_data, colors = module_colors)$eigengenes

if (args$merge) {
    
  # cut the tree at 0.20 (80% correlation)
  MEDissThres = 0.20
  
  merge = mergeCloseModules(expr_data, module_colors, cutHeight = MEDissThres, verbose = 3)
  module_colors = merge$colors
  MEs = merge$newMEs
}


####     sort phenotype data     ####

continuous_pheno <- read.csv(args$phenotype_continuous)
discrete_pheno <- read.csv(args$phenotype_discrete)

# sort on patient number then remove
continuous_pheno <- continuous_pheno[match(colnames(eset), continuous_pheno$patient_number), -1]
discrete_pheno <- discrete_pheno[match(colnames(eset), discrete_pheno$patient_number), -1]

# we will remove any phenotypes that have >75% missing data
continuous_pheno <- continuous_pheno[ , colSums(is.na(continuous_pheno)) < n_samples*0.75]
discrete_pheno <- discrete_pheno[ ,colSums(is.na(discrete_pheno)) < n_samples*0.75]

# reorder discrete and continuous pheno by correlating columns
cont_tree = hclust(as.dist(cor(data.matrix(continuous_pheno), use = "p", method = "pearson")), method="average")
disc_tree = hclust(as.dist(cor(data.matrix(discrete_pheno), use = "p", method = "spearman")), method = "average")

# reorder phenotypes by similarity
continuous_pheno <- continuous_pheno[ ,cont_tree$order]
discrete_pheno <- discrete_pheno[ ,disc_tree$order]


####     correlating gene modules with phenotypes     ####

# return any significant correlations as a data frame with unadjusted p value
continuous_significant <- significant_corrs(MEs, continuous_pheno, method = "pearson", fdr = args$fdr)
discrete_significant <- significant_corrs(MEs, discrete_pheno, method = "spearman", fdr = args$fdr)


####     save data for comparison with validation     #####

discovery_eset = eset
save(continuous_significant, discrete_significant, discovery_eset, 
     MEs, module_colors, file = args$out)

  