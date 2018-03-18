# latex_functions.R
# various functions to assist with latex formatting when compiling 
# manuscript.tex and manuscript.pdf from manuscript.Rnw

# add_part_label
# add a table part label on top of a table (as an extra row)
# df is dataframe to be printed using xtable
# labeltext is text to be printed
add_part_label <- function (df, labeltext) {

  part_label <- list ()
  part_label$pos <- list()
  part_label$pos[[1]] <- c(-1) # position label at top of table
  part_label$command <- c(paste0("\\multicolumn{", ncol(df)+1, "}{l}{", labeltext, "} \\\\ \n"))
  
  return(part_label)
  
}

# round_df
# round the numbers in a dataframe to the specified number of digits
round_df <- function (df, digits=3, col.select = NULL) {

  if (length(col.select) == 0) {
    df[sapply(df, class) == "numeric" ] <- sapply(df[sapply(df, class) == "numeric" ], round, digits)
  }
  
  df[,col.select] <- sapply(df[,col.select], round, digits)
    
  return(df)
}

# custom_cols
# add latex formatting / rename column names used in this manuscript
# this function is to be used as input to sanitize.colnames.function of print.xtable
custom_cols <- function (cols) {
  
  # Table 2
  # continuous phylogenetic signal
  cols <- gsub ("^lambda$", "$\\\\lambda$", cols)
  cols <- gsub ("^lambda\\.pval$", "$\\\\pval_\\\\lambda$", cols)
  cols <- gsub ("^K$", "\\\\K", cols)
  cols <- gsub ("^k\\.pval$", "$\\\\pval_\\\\K$", cols)
  
  # binary phylogenetic signal
  cols <- gsub ("^num_present$", "Number of presences", cols)
  cols <- gsub ("^num_absent$", "Number of absences", cols)
  cols <- gsub ("^D$", "\\\\D", cols)
  cols <- gsub ("^prob_random$", "\\\\pval\\\\textsubscript{rnd}", cols)
  cols <- gsub ("^prob_brownian$", "\\\\pval\\\\textsubscript{BM}", cols)
  
  # Table 3
  # quantitative PICs
  cols <- gsub ("^num_contrasts$", "Number of contrasts", cols)
  cols <- gsub ("^num_pos_con$", "Positive contrasts", cols)
  cols <- gsub ("^tval$", "\\\\tval", cols)
  cols <- gsub ("^pval$", "\\\\pval", cols)
  
  # Table 4
  # Pagel's correlated evolution
  cols <- gsub ("^logL_indep$", "LL (independent model)", cols)
  cols <- gsub ("^logL_dep$", "LL (dependent model)", cols)
  cols <- gsub ("^likelihood_ratio$", "Likelihood ratio", cols)
  
  # Table S1
  # continuous pglmmm columns
  cols <- gsub ("^Parameter.estimate$", "Estimate", cols)
  cols <- gsub ("^Lower.95..CI$", "Lower 95\\\\% CI", cols)
  cols <- gsub ("^Upper.95..CI$", "Upper 95\\\\% CI", cols)
  cols <- gsub ("^Effective.sample.size$", "Effective sample size", cols)
  cols <- gsub ("^P.value$", "\\\\pval", cols)
  
  # Table S2
  # binary pglmm columns
  cols <- gsub ("^sigmasq$", "$\\\\sigma^2$", cols)
  cols <- gsub ("^sigmap$", "\\\\pval($\\\\sigma^2=0$)", cols)
  cols <- gsub ("^coeff$", "Coefficient", cols)
  cols <- gsub ("^se$", "\\\\stderr", cols)
  cols <- gsub ("^zscore$", "\\\\zval-score", cols)
  cols <- gsub ("^pvalue$", "\\\\pval", cols)
  
  # Table S3
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
  cols <- gsub ("r\\.squared", "$\\\\rval^2$", cols)
  
  # Table S4
  # trait PCA
  cols <- gsub ("Std_PCA_Dim\\.1", "Standard PC1", cols)
  cols <- gsub ("Std_PCA_Dim\\.2", "Standard PC2", cols)
  cols <- gsub ("Phy_PCA_Dim\\.1", "Phylogenetic PC1", cols)
  cols <- gsub ("Phy_PCA_Dim\\.2", "Phylogenetic PC2", cols)
  
  return(cols)
}

# custom_rows
# add latex formatting / rename row names used in this manuscript
# this function is to be used as input to sanitize.rownames.function of print.xtable
custom_rows <- function (rows) {
  
  # continuous traits
  rows <- gsub ("stipe", "Stipe length", rows)
  rows <- gsub ("^length$", "Frond length", rows)
  rows <- gsub ("width", "Frond width", rows)
  rows <- gsub ("dissection", "Frond dissection", rows)
  rows <- gsub ("pinna", "Pinna number", rows)
  rows <- gsub ("sla", "Specific leaf area", rows)
  rows <- gsub ("rhizome", "Rhizome \\\\diameter", rows)
  
  # binary traits
  rows <- gsub ("habit", "Growth habit", rows)
  rows <- gsub ("epiphytic", "Epiphytic growth", rows)
  rows <- gsub ("gemmae", "Gemmae", rows)
  rows <- gsub ("glands", "Glands", rows)
  rows <- gsub ("hairs", "Hairs", rows)
  rows <- gsub ("cordate_morph", "Cordate morphotype", rows)
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
  
  # trait PCA
  rows <- gsub ("total_variance", "Total Variance (\\\\%)", rows)
  rows <- gsub ("cumulative_variance", "Cumulative Variance (\\\\%)", rows)
  
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

# embolden_top_values
# makes biggest (absolute) values in a column bold
# col.select are names of columns to be made bold, looks for all numeric columns by default
# row.select are names of columns to be made bold, uses all rows by default
# keep_top is number of top values to keep
# not run:
# bold_df <- embolden_top_values(PCA.result, row.select = c("dissection", "length", "pinna", "rhizome", "sla", "stipe", "width"))

embolden_top_values <- function (df, col.select = NULL, row.select = NULL, keep_top = 3) {
  # choose all columns of class numeric by default
  if (length(col.select) == 0) {
    col.select <- colnames(df)[which(apply(df, 2, class) == "numeric")]
  } 
  
  # choose all rows if not specified
  if (length(row.select) == 0) {
    row.select <- rownames(df)
  } 
  
  # first make list equal to length of selected cols, each with a vector indexing the top values
  col_list <- list()
  for (i in 1:length(col.select)) {
    col_list[[i]] <- order(abs(df[rownames(df) %in% row.select, col.select[i]]), decreasing=TRUE)[1:keep_top]
  }
  names(col_list) <- col.select
  
  # now loop through df and make these values bold
  for (i in 1:nrow(df)) {
    for (j in 1:ncol(df)) {
      if (colnames(df)[j] %in% col.select) {
        if (!is.na(df[i,j])) {
          if (i %in% col_list[[colnames(df)[j]]] ) {
            df[i,j] <- paste0("\\textbf{", df[i,j], "}")
          }
        }
      }
    }
  }
  
  return(df)
}
