#' @title densityCI
#' @description Visualising the CI using density plots of proportions
#'
#' @author Yingxin Lin
#' @param res results from \code{scDC_noClustering}, \code{scDC_clustering} function
#' @param condition a vector indicate the condition associate the results
#' @return Returns a data frame.
#'
#' @import ggplot2
#' @import ggridges
#' @export
#' @examples
#' ## Loading example data
#'
#' data("sim")
#'
#' cellTypes = sim$sim_cellTypes
#' subject = sim$sim_subject
#' \dontrun{
#' res_scDC = scDC_noClustering(cellTypes, subject,
#' calCI = TRUE, calCI_method = c("BCa", "percentile"))
#' densityCI(res_scDC, c("cond1","cond1","cond2","cond2"))
#'
#' }
#'


densityCI <- function(res, condition){
  df_toPlot <- reshape2::melt(res$thetastar)
  df_toPlot$cellTypes <- res$info[,"cellTypes"][df_toPlot$Var1]
  df_toPlot$subject <- res$info[,"subject"][df_toPlot$Var1]
  
  df_toPlot$cond <- condition
  
  conf_line <- res$results
  conf_line$cond <- condition
  
  conf_line$method <- factor(conf_line$method, levels = c("BCa", "percentile", "multinom"))
  n_method <- length(unique(conf_line$method))
  n_celltype = length(unique(df_toPlot$cellTypes))
  
  if (unique(df_toPlot$subject) == 1){
    g_density <- ggplot2::ggplot(df_toPlot, aes(x = value, y = subject, fill= cond)) +
      ggridges::stat_density_ridges(alpha = 0.5) +
      ggplot2::geom_vline(data = conf_line, aes(xintercept = conf_low, 
                                                color = method, linetype = method),
                          lwd = 1, alpha = 0.8) +
      ggplot2::geom_vline(data = conf_line, aes(xintercept = conf_high,
                                                color = method, linetype = method),
                          lwd = 1, alpha = 0.8) +
      ggplot2::theme_bw() +
      ggplot2::theme(axis.text.x = element_text(angle = 90), text = element_text(size = 12)) +
      ggplot2::scale_color_manual(values = .CIbarColor(n_method)) +
      ggplot2::scale_fill_brewer(palette = "Set2") +
      ggplot2::xlab("Proportion") +
      ggplot2::facet_wrap(~cellTypes, ncol = n_celltype , scales = "free_x")
    
  } else {
    
    g_density <- ggplot2::ggplot(df_toPlot, aes(x = value, y = subject, fill= cond)) +
      ggridges::stat_density_ridges(alpha = 0.5) +
      ggplot2::geom_segment(data = conf_line, aes(x = conf_low, xend = conf_low, y = as.numeric(subject),
                                                  yend = as.numeric(subject) + .9,
                                                  color = method, linetype = method),
                            lwd = 1, alpha = 0.8) +
      ggplot2::geom_segment(data = conf_line, aes(x = conf_high, xend = conf_high, y = as.numeric(subject),
                                                  yend = as.numeric(subject) + .9,
                                                  color = method, linetype = method),
                            lwd = 1, alpha = 0.8) +
      ggplot2::theme_bw() +
      ggplot2::theme(axis.text.x = element_text(angle = 90), text = element_text(size = 12)) +
      ggplot2::scale_color_manual(values = .CIbarColor(n_method)) +
      ggplot2::scale_fill_brewer(palette = "Set2") +
      ggplot2::xlab("Proportion") +
      ggplot2::facet_wrap(~cellTypes, ncol = n_celltype , scales = "free_x")
    
  }
  g_density
  return(g_density)
}
