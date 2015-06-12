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

library(yaml)
config = yaml.load_file("../config.yml")
creds <- list(username=config$credentials$username, password=config$credentials$pass)

# find all of the CD8 affy gene ST file names
files <- findAssayFiles(cell.type = "CD8", platform = "Affy GeneST", auth=creds)

# path on turing for file location: 
# /smith/smith/yet_to_be_sorted_data/

# make a wide search through all of /smithfor files with cd8 and tranplant tags
files.long <- list.files(path = "/smith", pattern = ".*CEL", all.files = FALSE, full.names = TRUE, recursive = TRUE, ignore.case = TRUE)
files.short <- sub(".*\\/([a-zA-Z0-9_ -\\[\\]]+.CEL)$", "\\1", files.long, perl=T)
all.files <- cbind(files.short, files.long)
all.files <- all.files[!duplicated(all.files[,1]),]
files.long <- all.files[,2]
files.short <- all.files[,1]
matched.files <- files.long[files.short %in% files]

## Read in files to celfile
library(oligo)
cels <- read.celfiles(matched.files)
cels <- reannoClinWeb(cels, auth=creds)
cels.bk <- cels

cd8_transplant = cels[,cels$studies == "Transplant" & cels$cell_type == "CD8"]

eset <- oligo::rma(cd8_transplant, target="core")

## create report QC report for all microarrays
arrayQualityMetrics(eset, outdir = "../data/normalized_transplant_eset_QC")



