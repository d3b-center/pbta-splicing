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

  ## create plot
  scatterplot <- ggscatter(expr_psi_df, 
                           x="IncLevel2", y="Expr", 
                           add = "reg.line", 
                           conf.int = TRUE, 
                           cor.coef = TRUE, 
                           cor.method = "pearson",
                           add.params = list(color = "red",
                                            fill = "pink"),
                                            ticks = TRUE,
                                            xlab = "Exon 4 Inclusion (PSI)", ylab = "CLK1 Expr (TPM)") + 
                          theme_Publication()
  
    return(scatterplot)
}