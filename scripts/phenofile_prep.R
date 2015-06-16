# turn the csv phenotype files into a usable state

pheno_data <- read.csv("../data/full_phenotype_data.csv", header=TRUE)

# remove columns with all NAs
c_na = which(colSums(is.na(pheno_data)) == nrow(pheno_data))

# remove columns with all NULL
c_null = which(colSums(pheno_data == "NULL") == nrow(pheno_data))

pheno_data <- pheno_data[, -c(c_na, c_null)]

# get continuous data of interest
continuous_colnames <- c("Patient.number", "Length.of.RRT..months.", "CIT..hours.", "Number.of.MM", "Number.of.DR.MM", "IgG.level", "CMV_IgG", "K...pre.tx.", "eGFR.pre.tx")
tcell_ip_colnames <- names(pheno_data)[grep("^CD",names(pheno_data))]
continuous_pheno <- pheno_data[, c(continuous_colnames, tcell_ip_colnames)]

# get discrete data of interest
discrete_pheno <- pheno_data[, c("Patient.number", "Final.DGF.decision", "Immunological.risk", "Sensitised..Y.N.", "Final.biopsy.decision", "Infectious.complication", "CMV.donor", "CMV.recipient", "Rejection.Days")]

### transforming discrete data to ordered/binary ###
# 0 is none, 1 is DGF (Delayed Graft Function) - TODO: ask difference between DGF, functional DGF
discrete_pheno$Final.DGF.decision[discrete_pheno$Final.DGF.decision == "functional DGF"] = "DGF"
discrete_pheno$Final.DGF.decision <- factor(discrete_pheno$Final.DGF.decision, levels = c("none", "DGF"), ordered = TRUE)

# 0 is low, 1 is I/H (intermediate/high?)
discrete_pheno$Immunological.risk <- factor(discrete_pheno$Immunological.risk, levels = c("Low", "I/H"), ordered = TRUE )

# 0 is N, 1 is Y - TODO: ask for raw score on this (percentage of autoreactivity to panel?)
discrete_pheno$Sensitised..Y.N.[discrete_pheno$Sensitised..Y.N. %in% c("Y (highly)", "Y(highly)")] = "Y"
discrete_pheno$sensitised <- factor(discrete_pheno$Sensitised..Y.N., levels = c("N", "Y"), ordered = TRUE)

# true/false rejection at all
discrete_pheno$rejection = grepl("Rejection", discrete_pheno$Final.biopsy.decision)

# count the number of rejection events in each bin 0 to 1wk, 0 to 2wk, 1mo, 3mo, 1yr
discrete_pheno$rejection_1wk <- sapply(pheno_data$Rejection.Days, function(s) sum(as.integer(strsplit(as.character(s),",")[[1]]) <= 7))
discrete_pheno$rejection_2wk <- sapply(pheno_data$Rejection.Days, function(s) sum(as.integer(strsplit(as.character(s),",")[[1]]) <= 14))
discrete_pheno$rejection_1month <- sapply(pheno_data$Rejection.Days, function(s) sum(as.integer(strsplit(as.character(s),",")[[1]]) <= 31))
discrete_pheno$rejection_3month <- sapply(pheno_data$Rejection.Days, function(s) sum(as.integer(strsplit(as.character(s),",")[[1]]) <= 93))
discrete_pheno$rejection_1yr <- sapply(pheno_data$Rejection.Days, function(s) sum(as.integer(strsplit(as.character(s),",")[[1]]) <= 365))

# split Infectious.complication to viral_complication, bacterial_complication, fungal_complication
discrete_pheno$viral_complication = grepl("Viral", discrete_pheno$Infectious.complication)
discrete_pheno$bacterial_complication = grepl("Bacterial", discrete_pheno$Infectious.complication)
discrete_pheno$fungal_complication = grepl("Fungal", discrete_pheno$Infectious.complication)

# combine CMV donor/recipient into negative/negative (best), negative/positive (second best), positive/positive (third best), positive, negative (worst)
# if else is chosen such that above becomes 1, 2, 3, 4
discrete_pheno$CMV = ifelse(discrete_pheno$CMV.donor == "positive", 
                            ifelse(discrete_pheno$CMV.recipient == "positive", 3, 4),
                            ifelse(discrete_pheno$CMV.recipient == "negative", 1, 2))

discrete_pheno <- discrete_pheno[,c("Patient.number", "Final.DGF.decision", "Immunological.risk", 
                                    "sensitised", "rejection", "rejection_1wk", "rejection_2wk", 
                                    "rejection_1month", "rejection_3month", "rejection_1yr", 
                                    "CMV", "viral_complication", "fungal_complication", 
                                    "bacterial_complication")]

write.csv(continuous_pheno, "../data/transplant_continuous_phenotypes.csv", quote = FALSE, row.names = FALSE)
write.csv(discrete_pheno, "../data/transplant_discrete_phenotypes.csv", quote = FALSE, row.names = FALSE)


