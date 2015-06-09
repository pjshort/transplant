# script to load microarray data for transplant project

# what the code accomplishes:
# 1. search file system for any array corresponding to CD8 cells in transplant patients
# 2. annotate with phenotype information (requires X11 window)
# 3. normalize read intensities and calculate expression set with oligo::rma
# 4. filter down to one probe per gene (on largest IQR), only probes with gene annotations from hugene

# INPUT:
# --search_path -> location in file system to begin recursive search for *.CEL files (defaults to current directory)
# --out -> name for the resulting RData file with the annotated ExpressionSet object
#           this should be used for WGCNA analysis, QC, etc.
# optional: --cell_type -> cell type (defaults to "CD8")
# optional: --study -> study (defaults to "Transplant")
# optional: --time_point -> time_point (defaults to "0 months")
# optional: --platform -> defaults to Affy GeneST

# OUTPUT: RData file with the ExpressionSet object called eset with expression and phenotype information

# documentation notes:
# a triple hash (### description xyz) denotes a 'section header' in the code while (# comments..) denotes a more simple comment

### dependencies
creds <- list(username="***", password="***")
library(optparse)
library(ggplot2)
library(oligo)
library(Biobase)
library(genefilter)
library(arrayQualityMetrics)
library(SmithLabArray)
library(ClinStudyWeb)
library(WGCNA)
library(flashClust)
library(hugene10sttranscriptcluster.db)
library(genefilter)
options(stringsAsFactors = FALSE)
allowWGCNAThreads()

### command line options
option_list <- list(
  make_option("--search_path", default="../data/TransplantCELs/batch1_2009",
              help="Top level directory to start recursive search for *.CEL files."),
  make_option("--out", default="../data/transplant_eset.Rdata",
              help="Location to save the ExpressionSet object. Defaults to ./eset.Rdata"),
  make_option("--cell_type", default="CD8",
              help="Restrict by cell type (e.g. CD8, CD4, etc). Defaults to CD8."),
  make_option("--study", default="Transplant",
              help="Restrict by study. Defaults to Transplant."),
  make_option("--time_point", default="0 months",
              help="Restrict by time point Defaults to 0 months."),
  make_option("--platform", default="Affy GeneST",
              help="Microarray platform. Defaults to Affy GeneST"),
  make_option("--verbose", action="store_true", default=FALSE,
              help="Print extra output advising the user of progression through the analysis.")
)

args <- parse_args(OptionParser(option_list=option_list))

### LOAD data and filter for study type, cell type, time point, platform

# use ClinStudyWeb to locate all of the files in database with specified platform and cell type
CSW_files <- findAssayFiles(cell.type = args$cell_type, platform = args$platform, auth=creds)

# path on turing for file location: 
# /smith/smith/yet_to_be_sorted_data/

# grab path to any .CEL file beneath the search_path specified, remove duplicates, check for matches to files
files.long <- list.files(path = args$search_path, pattern = ".*CEL", all.files = FALSE, full.names = TRUE, recursive = TRUE, ignore.case = TRUE)
files.short <- sub(".*\\/([a-zA-Z0-9_ -\\[\\]]+.CEL)$", "\\1", files.long, perl=T)
all.files <- cbind(files.short, files.long)
all.files <- all.files[!duplicated(all.files[,1]),]
files.long <- all.files[,2]
files.short <- all.files[,1]
matched.files <- files.long[files.short %in% CSW_files]

# read .CEL files
cels <- read.celfiles(matched.files)

# annotate .CEL files with phenotype and metadata
cels <- reannoClinWeb(cels, auth=creds)

# keep only the cels that meet the specified criteria for cell type, study design, etc.
cels_pass = cels[,cels$studies == args$study & cels$cell_type == args$cell_type & cels$time_point == args$time_point]

# create ExpressionSet object using RMA from the CEL files
eset <- oligo::rma(cels_pass, target="core")

# get the phenotype data from the eset and filter out any with all NAs
pheno_data = eset@phenoData@data
pheno_data <- pheno_data[,colSums(is.na(pheno_data)) < nrow(pheno_data)]


### Annotate probes with entrezIDs to remove probes that will not serve in downstream analysis

# annotate the ExpressionSet object
fData <- data.frame(ID = featureNames(eset))
rownames(fData) <- fData$ID

# official gene name
ann <- mget(as.character(fData$ID), hugene10sttranscriptclusterSYMBOL, ifnotfound=NA)
ann <- lapply(ann, paste, collapse=";")
fData$Name <- unlist(ann)

# Entrez ID
ann <- mget(as.character(rownames(fData)),hugene10sttranscriptclusterENTREZID, ifnotfound=NA)
ann <- lapply(ann, paste, collapse=";")
fData$EntrezID <- unlist(ann)

# chromosome
ann <- mget(as.character(rownames(fData)), hugene10sttranscriptclusterCHR, ifnotfound=NA)
ann <- lapply(ann, paste, collapse=";")
fData$Chromosome <- unlist(ann)

# full gene name
ann <- mget(as.character(rownames(fData)), hugene10sttranscriptclusterGENENAME, ifnotfound=NA)
ann <- lapply(ann, paste, collapse=";")
fData$LongName <- unlist(ann)

# Affymetrix probe status annotation
# download this file from http://www.affymetrix.com/support/technical/annotationfilesmain.affx
affy.annot <- read.csv("../data/HuGene-1_0-st-v1.na35.hg19.transcript.csv", skip=21, header=T)

# use Affy annotation to get probes with annotation information
rownames(affy.annot) <- affy.annot$transcript_cluster_id
fData$ProbeStatus <- factor(as.character(affy.annot[rownames(fData),"category"]))
metadata <- data.frame(labelDescription = c("Manufacturers ID", "Official Symbol", "EntrezID", 
                                            "Chromosome", "Gene Name", "Affy Probe Status"), 
                       row.names=c("ID", "Name", "EntrezID", "Chromosome", "LongName", "ProbeStatus"))

colnames(fData) <- c("ID", "Name", "EntrezID", "Chromosome", "LongName", "ProbeStatus")
features <- new("AnnotatedDataFrame", data = fData, varMetadata = metadata)
featureData(eset) <- features

### remove any probes without gene annotations
entrezIds <- mget(featureNames(eset), envir = hugene10sttranscriptclusterENTREZID, ifnotfound=NA)
haveEntrezId <- names(entrezIds)[sapply(entrezIds, function(x) !is.na(x))]
numNoEntrezId <- length(featureNames(eset)) - length(haveEntrezId) 
eset <- eset[haveEntrezId, ]

# Then, make sure each probe only maps to 1 entrezID
esIqr <- apply(exprs(eset), 1, IQR)
uniqGenes <- findLargest(featureNames(eset), esIqr, "hugene10sttranscriptcluster")
eset <- eset[uniqGenes, ]
numSelected <- length(featureNames(eset))

# make gene symbol the featureName - makes more sense to work with
featureNames(eset) <- fData(eset)$Name
eset <- eset[fData(eset)$ProbeStatus == "main", ]

save(eset, file = args$out)
