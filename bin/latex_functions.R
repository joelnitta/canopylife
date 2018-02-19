# latex_functions.R
# various functions to assist with latex formatting when compiling 
# manuscript.tex and manuscript.pdf from manuscript.Rnw

# custom_cols
# add latex formatting / rename column names used in this manuscript
# this function is to be used as input to sanitize.colnames.function of print.xtable
custom_cols <- function (cols) {
  
  # phylo signal columns
  cols <- gsub ("^lambda$", "$\\\\lambda$", cols)
  cols <- gsub ("^lambda\\.pval$", "$\\\\pval_\\\\lambda$", cols)
  cols <- gsub ("^K$", "\\\\K", cols)
  cols <- gsub ("^k\\.pval$", "$\\\\pval_\\\\K$", cols)
  cols <- gsub ("^D$", "\\\\D", cols)
  cols <- gsub ("^prob_random$", "\\\\pval\\\\textsubscript{BM}", cols)
  cols <- gsub ("^prob_brownian$", "\\\\pval\\\\textsubscript{rnd}", cols)
  
  # binary pglmm columns
  cols <- gsub ("^sigmasq$", "$\\\\sigma^2$", cols)
  cols <- gsub ("^sigmap$", "$\\\\pval(\\\\sigma^2=0)$", cols)
  cols <- gsub ("^coeff$", "Coefficient", cols)
  cols <- gsub ("^se$", "\\\\stderr", cols)
  cols <- gsub ("^zscore$", "\\\\zval-score", cols)
  cols <- gsub ("^pvalue$", "\\\\pval", cols)
  
  # microclimate ancova columns
  cols <- gsub ("df\\.intercept", "\\\\df(intercept)", cols)
  cols <- gsub ("f\\.intercept", "\\\\fval(intercept)", cols)
  cols <- gsub ("pval\\.intercept", "\\\\pval(intercept)", cols)
  cols <- gsub ("df\\.slope", "\\\\df(slope)", cols)
  cols <- gsub ("f\\.slope", "\\\\fval(slope)", cols)
  cols <- gsub ("pval\\.slope", "\\\\pval(slope)", cols)
  
  # microclimate slopes
  cols <- gsub ("slope\\.epi", "Slope (epiphytic)", cols)
  cols <- gsub ("inter\\.epi", "Intercept (epiphytic)", cols)
  cols <- gsub ("slope\\.ter", "Slope (terrestrial)", cols)
  cols <- gsub ("inter\\.ter", "Intercept (terrestrial)", cols)
  
  return(cols)
}

# custom_rows
# add latex formatting / rename row names used in this manuscript
# this function is to be used as input to sanitize.rownames.function of print.xtable
custom_rows <- function (rows) {
  
  # traits
  rows <- gsub ("habit", "Growth habit", rows)
  rows <- gsub ("stipe", "Stipe length", rows)
  rows <- gsub ("^length$", "Frond length", rows)
  rows <- gsub ("width", "Frond width", rows)
  rows <- gsub ("dissection", "Frond dissection", rows)
  rows <- gsub ("pinna", "Pinna number", rows)
  rows <- gsub ("sla", "Specific leaf area", rows)
  rows <- gsub ("rhizome", "Rhizome \\\\diameter", rows)
  rows <- gsub ("gemmae", "Gemmae", rows)
  rows <- gsub ("glands", "Glands", rows)
  rows <- gsub ("hairs", "Hairs", rows)
  rows <- gsub ("morph_binary", "Morphotype", rows)
  
  # microclimate
  rows <- gsub ("max_temp", "Max. temperature", rows)
  rows <- gsub ("mean_temp", "Mean temperature", rows)
  rows <- gsub ("min_temp", "Min. temperature", rows)
  rows <- gsub ("sd_temp", "\\\\stdev temperature", rows)
  rows <- gsub ("max_RH", "Max. Rel. Hum.", rows)
  rows <- gsub ("mean_RH", "Mean Rel. Hum.", rows)
  rows <- gsub ("min_RH", "Min. Rel. Hum.", rows)
  rows <- gsub ("sd_RH", "\\\\stdev Rel. Hum.", rows)
  
  return(rows)
}


# embolden_p
# makes p-values less than chosen significance level bold in results table
# also replaces values of 0 with < 0.001 etc where 1 is placed at number of digits (3 by default)
# df is dataframe to be printed using xtable
# col.select are names of columns to be made bold, looks for pvalues by default
embolden_p <- function (df, col.select="pvalue", num.digits=3, sigval = 0.05) {
  for (i in 1:nrow(df)) {
    for (j in 1:ncol(df)) {
      if(!is.na(df[i,j])) {
        if (df[i,j] == 0 && colnames(df)[j] %in% col.select) {
          df[i,j] <- paste0("\\textbf{\\textless", 0.1^num.digits, "}")
        } else if (df[i,j] < sigval && colnames(df)[j] %in% col.select) {
          df[i,j] <- paste0("\\textbf{", df[i,j], "}")
        }
      }
    }
  }
  return(df)
}
