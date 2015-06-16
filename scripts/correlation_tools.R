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
  l_mat = which(p_adjust < fdr, arr.ind = TRUE)
  l_flat = which(p_adjust < fdr)
  modules = names(MEs)[l_mat[,1]]
  traits = colnames(phenotypes)[l_mat[,2]]
  cs = corr[l_flat]
  ps = p[l_flat]
  
  associations <- data.frame(cbind(modules, traits, cs, ps))
  return(associations)
  
}