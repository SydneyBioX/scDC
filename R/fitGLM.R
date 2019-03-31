#' @title fitGLM
#' @description Visualising the CI using barplot
#'
#' @author Yingxin Lin
#' @param res results from \code{scDC_noClustering}, \code{scDC_clustering} function
#' @param condition a vector indicate the condition associate the results
#' @param subject_effect A logical input for whether fit the subject effect
#' @param pairwise A logical input for whether the subject in different condition are paired
#' @param fixed_only A logical input for whether only fiftting the fixed effect GLM model.
#' @param verbose A logical input for whether print the progress.
#' @return return GLM results
#'
#' @import mice
#' @import lme4
#' @export
#' @examples
#' ## Loading example data
#'

#' \dontrun{
#' library(scDC)
#' data("sim")
#'
#' cellTypes = sim$sim_cellTypes
#' subject = sim$sim_subject
#' res_scDC = scDC_noClustering(cellTypes, subject,
#' calCI = TRUE, calCI_method = c("BCa", "percentile"))
#' barplotCI(res_scDC, c("cond1","cond1","cond2","cond2"))
#'
#' }
#'

fitGLM <- function(res, condition, subject_effect = TRUE, pairwise = TRUE, fixed_only = FALSE, verbose = TRUE){

  fit_random <- list()
  fit_fixed <- list()
  for(i in 1:ncol(res$nstar)){
    # idx <- indexes_list[[i]]
    if(verbose){
      if(i%%10==0){
        print(paste("fitting GLM...", i))
      }
    }


    glm_df <-  cbind(res$info[,1:2], res$nstar[,i])

    # glm_df <- melt(glm_df)
    colnames(glm_df) <- c("cellTypes", "subject", "cell_count")
    glm_df$cond <- condition

    if(subject_effect){
      if(pairwise){

        fit_fixed[[i]] <- stats::glm(cell_count ~ cellTypes + cond +  cellTypes:cond + subject,
                                     data = glm_df,
                                     family = poisson(link=log))
        if(!fixed_only){
          fit_random[[i]] <- lme4::glmer(cell_count ~ cellTypes + cond +  cellTypes:cond + (1 | subject ),
                                         data = glm_df, family = poisson(link=log),
                                         control = glmerControl(nAGQ = 0L))
        }
      }else{
        fit_fixed[[i]] <- stats::glm(cell_count ~ cellTypes + cond +  cellTypes:cond + subject,
                                     data = glm_df, family = poisson(link=log))
        if(!fixed_only){

          fit_random[[i]] <- lme4::glmer(cell_count ~ cellTypes + cond +
                                           cellTypes:cond + (1 | subject ), data = glm_df,
                                         family = poisson(link=log),
                                         control = glmerControl(nAGQ = 0L))
        }
      }
    }else{
      fit_fixed[[i]] <- stats::glm(cell_count ~ cellTypes + cond +  cellTypes:cond, data = glm_df,
                                   family = poisson(link=log))
    }

  }

  if(!subject_effect){
    fixed_only = TRUE
  }

  if(!fixed_only){
    pool_res_random = mice::pool(fit_random)
    pool_res_fixed = mice::pool(fit_fixed)
    return(list(pool_res_random = pool_res_random,
                pool_res_fixed = pool_res_fixed,
                fit_random = fit_random,
                fit_fixed = fit_fixed))
  }else{
    pool_res_fixed = mice::pool(fit_fixed)
    return(list(pool_res_fixed = pool_res_fixed,
                fit_fixed = fit_fixed))
  }

}
