#' @title scDC_clustering
#' @description Single-cell Differential Composition Analysis with performing clustering
#'
#' @author Yingxin Lin
#' @param exprsMat logcounts expression matrix with each row represents gene, and each column represents cell
#' @param cellTypes A vector indicates the cell type info of the data
#' @param subject A vector indicates the subject info of the data
#' @param calCI A logical input for whether calculating the confidence interval for proportion
#' @param calCI_method A string indicates the method that is used to calculate confidence interval. Options include \code{BCa}, \code{percentile}, and \code{multinom}.
#' @param nboot Number of bootstrap. If \code{calCI = TRUE}, \code{nboot = 10000} by default. Otherwise, \code{nboot = 500}.
#' @param conf_level confidence level, with default 0.95
#' @param ncores Number of cores that are used.
#' @param verbose A logical input for whether print the progress.
#' @return Returns a data frame.
#'
#' @import parallel
#' @import reshape2
#' @import DescTools
#' @import scNet
#' @import stats
#' @export
#' @examples
#' ## Loading example data
#'
#' data("sim")
#'
#' cellTypes = sim$sim_cellTypes
#' subject = sim$sim_subject
#' \dontrun{
#' res_noCALCI = scDC_clustering(cellTypes, subject, calCI = FALSE)
#' res_BCa = scDC_clustering(cellTypes, subject, calCI = TRUE, calCI_method = "BCa")
#' res_percentile = scDC_clustering(cellTypes, subject, calCI = TRUE, calCI_method = "percentile")
#' res_multinom = scDC_clustering(cellTypes, subject, calCI = TRUE, calCI_method = "multinom")
#'
#' }
#'




scDC_clustering <- function (exprsMat = NULL,
                             cellTypes = NULL,
                             subject = NULL,
                             calCI = TRUE,
                             calCI_method = c("BCa", "multinom", "percentile"),
                             nboot = NULL,
                             conf_level = 0.95,
                             ncores = 1,
                             verbose = TRUE)
{

  x <- 1:ncol(exprsMat)
  n <- length(x)

  if(length(cellTypes)!=length(subject)){
    stop("the vector length of cell type info and subject info don't match!")
  }

  if(calCI){
    calCI_method <- match.arg(calCI_method, choices = c("BCa", "multinom",
                                                        "percentile"), several.ok = TRUE)
    if(is.null(nboot)){
      warnings("number of bootstrap is set as 10000")
      nboot = 10000
    }

  }

  if(!calCI&is.null(nboot)){
    warnings("number of bootstrap is set as 100")
    nboot = 100
  }

  cellTypes <- as.character(cellTypes)
  subject <- as.character(subject)

  tab <- table(cellTypes, subject)
  info <-  reshape2::melt(tab)

  if(verbose){
    print("Calculating sample proportion...")
  }
  df <- data.frame(cellTypes = cellTypes, subject = subject)
  # thetahat <- .bootstrap_clustering(x, exprsMat, cellTypes, subject, verbose = verbose)
  calProp_hat <- .bootstrap_clustering(x, exprsMat, cellTypes, subject, verbose = verbose)
  thetahat <- calProp_hat$prop
  nhat <- calProp_hat$count


  if(verbose){
    print("Calculating bootstrap proportion...")
  }
  bootsam <- do.call(rbind, .stratifiedBootstrap(1:length(cellTypes),
                                                strata = subject, times = nboot))

  # thetastar <-  do.call(cbind, parallel::mclapply(1:nboot, function(i){
  #   .bootstrap_clustering(bootsam[i,], exprsMat, cellTypes, subject, verbose = verbose)
  # }, mc.cores = ncores))
  #
  calProp_star <- parallel::mclapply(1:nboot, function(i){
    .bootstrap_clustering(bootsam[i,], exprsMat, cellTypes, subject, verbose = verbose)
  }, mc.cores = ncores)

  thetastar <-  do.call(cbind, lapply(calProp_star, "[[", "prop"))
  nstar <- do.call(cbind, lapply(calProp_star, "[[", "count"))


  res_multinom <- NULL
  res_BCa <- NULL
  res_percentile <- NULL

  if(calCI){
    alpha <-  c((1-conf_level)/2, 1-(1-conf_level)/2)
    if (!all(alpha < 1) || !all(alpha > 0))
      stop("All elements of alpha must be in (0,1)")
    alpha_sorted <- sort(alpha)
    if (nboot <= 1/min(alpha_sorted[1], 1 - alpha_sorted[length(alpha_sorted)]))
      warning("nboot is not large enough to estimate your chosen alpha.")

    if(verbose){
      print(paste("Calculating", calCI_method, "..."))
    }
    if("multinom" %in% calCI_method){


      ### Calculating the CI based on the multinomial


      confpoints = lapply(1:ncol(tab), function(i){
        multi_res <- DescTools::MultinomCI(tab[,i], conf.level = conf_level)
        multi_res <- multi_res[,-1]
        colnames(multi_res) <- c("conf_low", "conf_high")
        multi_res
      })
      res_multinom <- info[,1:2]
      res_multinom <- cbind(res_multinom, do.call(rbind, confpoints))
      rownames(res_multinom) <- 1:nrow(res_multinom)
      res_multinom$method <- "multinom"

      # return(list(results = res,
      #             confpoints = confpoints,
      #             thetastar = thetastar,
      #             thetahat = thetahat,
      #             nstar = nstar,
      #             nhat = nhat,
      #             info = info))

    }

    if("BCa" %in% calCI_method){

      if(verbose){
        print("Calculating z0 ...")
      }

      z0 <- sapply(1:nrow(thetastar), function(i) qnorm(sum(thetastar[i,] < thetahat[i])/nboot))
      if(verbose){
        print("Calculating acc ...")
      }

      u <- list()
      u <- parallel::mclapply(1:n, function(i){
        .bootstrap_clustering(x[-i], exprsMat, cellTypes, subject, verbose = verbose)$prop
      }, mc.cores = ncores)

      u <- do.call(cbind, u)
      uu <- matrix(rep(rowMeans(u), ncol(u)), ncol = ncol(u)) - u
      acc <- apply(uu, 1, function(x) sum(x * x * x)/(6 * (sum(x * x))^1.5))

      zalpha <- qnorm(alpha)
      tt <- list()
      for(i in 1:length(z0)){
        tt[[i]] <- pnorm(z0[i] + (z0[i] + zalpha)/(1 - acc[i] * (z0[i] + zalpha)))
      }

      confpoints <- lapply(1:length(tt), function(i) quantile(x = thetastar[i,], probs = tt[[i]], type = 1))
      confpoints_original <- confpoints
      confpoints <- lapply(confpoints, function(x) {
        names(x) <- NULL
        x})
      confpoints <- lapply(confpoints, function(x) cbind(alpha, x))
      confpoints <- lapply(confpoints, function(x){
        dimnames(x)[[2]] <- c("alpha", "bca point")
        x
      })

      res_BCa <- info[,1:2]
      res_BCa$conf_low <- do.call(rbind, lapply(confpoints, function(x)x[1,2]))
      res_BCa$conf_high <- do.call(rbind, lapply(confpoints, function(x)x[2,2]))
      colnames(res_BCa) <- c("cellTypes", "subject", "conf_low", "conf_high")
      res_BCa$method <- "BCa"
      # return(list(results = res,
      #             confpoints = confpoints,
      #             z0 = z0,
      #             acc = acc,
      #             u = u,
      #             thetastar = thetastar,
      #             thetahat = thetahat,
      #             nstar = nstar,
      #             nhat = nhat,
      #             info = info,
      #             confpoints_original = confpoints_original,
      #             tt = tt))
    }

    if("percentile" %in% calCI_method){

      confpoints <- apply(thetastar, 1, function(x)quantile(x, probs = alpha))

      res_percentile <- info[,1:2]
      res_percentile$conf_low <- confpoints[1,]
      res_percentile$conf_high <- confpoints[2,]
      res_percentile$method <- "percentile"
      #
      #         return(list(results = res,
      #                     confpoints = confpoints,
      #                     thetastar = thetastar,
      #                     thetahat = thetahat,
      #                     nstar = nstar,
      #                     nhat = nhat,
      #                     info = info))

    }



  }

  res <- rbind(res_BCa, res_percentile, res_multinom)


  return(list(results = res,
              thetastar = thetastar,
              thetahat = thetahat,
              nstar = nstar,
              nhat = nhat,
              info = info))


}
