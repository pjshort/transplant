# prepare .GCT and .CLS files for Gene Set Enrichment Analysis

# load the expression set
load("../results/eset_batch23.RData")

# use array tools
library(ArrayTools)

phenoData(eset) = AnnotatedDataFrame(discrete_pheno)

createGSEAFiles(mydir = getwd(), eset, "rejection")
