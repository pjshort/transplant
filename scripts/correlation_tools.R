# tools for assessing correlations between gene modules and phenotypes

significant_corrs <- function(MEs, phenotypes, method = "pearson", fdr = 0.05){
  
  # takes module eigengene expression levels and phenotype data matrix
  # returns dataframe of module to phenotype correlations with correlation strength and p value
  
  corr = WGCNA::cor(MEs, phenotypes, use = "p", method = method)
  
  # get p value from correlation and degrees of freedom
  p = corPvalueStudent(corr, n_samples)
  
  # benjamini-hochberg error rate correction
  p_adjust = apply(p, MARGIN = 2, FUN = function(p) p.adjust(p, length(corr), method = "fdr"))
  
  # get row and column with adjusted p < false discovery rate
  l = which(p_adjust < fdr, arr.ind = TRUE)
  modules = names(MEs)[l[,1]]
  traits = colnames(phenotypes)[l[,2]]
  
  associations <- data.frame(cbind(modules, traits, corr, p))
  return(associations)
  
}