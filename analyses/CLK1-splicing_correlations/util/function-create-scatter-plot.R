###--------------------------------------------------------------
#' Create scatter plot for CLK1 Exon 4 splicing and expression
#'
#' @param expr_psi_df 
#'
#' @return
#' @export
#'
#' @examples
create_scatterplot <- function(expr_psi_df) {

  ## first log2 expression
  expr_psi_df_log <- expr_psi_df %>%
    mutate(logExp = log(Expr, 2))
  
  goi <- unique(expr_psi_df_log$geneSymbol)
 
   ## create plot
  scatterplot <- ggscatter(expr_psi_df_log, 
                           x="IncLevel1", 
                           y="logExp", 
                           add = "reg.line", 
                           conf.int = TRUE, 
                           cor.coef = TRUE, 
                           cor.method = "pearson",
                           add.params = list(color = "red",
                                            fill = "pink"),
                           ticks = TRUE) + 
                           xlab(expression(bold(bolditalic("CLK1")~"Exon 4 Inclusion (PSI)"))) +
                           ylab(substitute(bold(bolditalic(var_name)~"Expression (log2 TPM)"), list(var_name = goi))) +
                          theme_Publication()
  
    return(scatterplot)
}


