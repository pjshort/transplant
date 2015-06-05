# search the file system for CD8 transplant arrays


# Initialize the environment:
rm(list=ls())
library(devtools)
library(R.utils)
library(oligo)
library(Biobase)
library(genefilter)
library(arrayQualityMetrics)
library(SmithLabArray)
library(ClinStudyWeb)
creds <- list(username="pjs90", password="pjs2810")

# find all of the CD8 affy gene ST file names
files <- findAssayFiles(cell.type = "CD8", platform = "Affy GeneST", auth=creds)

# path on turing for file location: 
# /smith/smith/yet_to_be_sorted_data/

files.long <- list.files(path = "/smith/smith/yet_to_be_sorted_data/", pattern = ".*CEL", all.files = FALSE, full.names = TRUE, recursive = TRUE, ignore.case = TRUE)
files.short <- sub(".*\\/([a-zA-Z0-9_ -\\[\\]]+.CEL)$", "\\1", files.long, perl=T)
all.files <- cbind(files.short, files.long)
all.files <- all.files[!duplicated(all.files[,1]),]
files.long <- all.files[,2]
files.short <- all.files[,1]
matched_files <- c(matched_files, files.long[files.short %in% files])

## Read in files to celfile
library(oligo)
cels <- read.celfiles(matched.files)
cels <- reannoClinWeb(cels, auth=creds)
cels.bk <- cels

## TA for these cell_types only
#cels <- cels[, grepl("(CD14|PBMC|CD19|CD8)", cels$cell_type)]

## assert that all CELs correspond to transplant and CD8
if (!all(cels$cell_type == "CD8")){ stop("CELs are not all CD8.")}
if (!all(cels$studies == "Transplant")){ stop("CELs are not all from transplant study.")}

eset <- oligo::rma(cels, target="core")

## create report QC report for all microarrays
arrayQualityMetrics(eset, outdir = "../data/normalized_transplant_eset_QC")



