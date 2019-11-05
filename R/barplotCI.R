#' @title barplotCI
#' @description Visualising the CI using barplot
#'
#' @author Yingxin Lin
#' @param res results from \code{scDC_noClustering}, \code{scDC_clustering} function
#' @param condition a vector indicate the condition associate the results
#' @return Returns a data frame.
#'
#' @import ggplot2
#' @export
#' @examples
#' ## Loading example data
#' library(scDC)
#' data("sim")
#'
#' cellTypes = sim$sim_cellTypes
#' subject = sim$sim_subject
#' \dontrun{
#' res_scDC = scDC_noClustering(cellTypes, subject,
#' calCI = TRUE, calCI_method = c("BCa", "percentile"))
#' barplotCI(res_scDC, c("cond1","cond1","cond2","cond2"))
#'
#' }
#'


barplotCI <- function(res, condition){
  df_toPlot <- res$results
  df_toPlot$median <- apply(res$thetastar, 1, median)
  df_toPlot$cond <- condition
  df_toPlot$method <- factor(df_toPlot$method, levels = c("BCa", "percentile", "multinom"))
  n_method <- length(unique(df_toPlot$method))
  
  n_celltype = length(unique(df_toPlot$cellTypes))
  
  g_bar <- ggplot2::ggplot(df_toPlot, aes(x = subject, y = median, fill = cond)) +
    ggplot2::geom_bar(stat="identity", position = "dodge", alpha = 0.8) +
    ggplot2::theme_bw() +
    ggplot2::scale_fill_brewer(palette = "Set2") +
    ggplot2::ylab("Proportion") +
    ggplot2::geom_errorbar(aes(ymin=conf_low, ymax=conf_high, color = method), width=.3,lwd = 1,
                  position=position_dodge(width = 0.5)) +
    ggplot2::theme(axis.text.x = element_text(angle = 90), text = element_text(size = 12)) +
    ggplot2::scale_color_manual(values = .CIbarColor(n_method)) +
    ggplot2::facet_wrap(~cellTypes, ncol = n_celltype，
                        labeller = labeller(cellTypes  = label_wrap_gen(width = 10,  multi_line = TRUE))) +
    ggplot2::coord_flip()+
    ggplot2::ylim(c(0,1))+
    NULL

  g_bar
  return(g_bar)


}
