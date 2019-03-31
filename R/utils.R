

######################################################################
# Helper function
######################################################################

.stratifiedBootstrap <- function(idx, strata, times = 1000){
  strata_list <- unique(strata)
  index_list <- lapply(1:times, function(x){
    unlist(lapply(strata_list, function(s){
      idx_strata = idx[strata == s]
      sample(idx_strata, length(idx_strata), replace = T)
    }))
  })
  return(index_list)
}



.calculateProp <- function(x, xdata){
  tab <- table(xdata[x,]$cellTypes, xdata[x,]$subject)
  tab_prop <- c(tab/matrix(rep(colSums(tab), nrow(tab)), nrow = nrow(tab), byrow = T))
  tab_count <- c(tab)
  return(list(count = tab_count, prop = tab_prop))
}

.scran_high_var <- function(exprsMat,topn=1000){
  topn <- min(topn, nrow(exprsMat))
  var.fit <- scran::trendVar(exprsMat, method="loess")
  var.out <- scran::decomposeVar(exprsMat, var.fit)
  hvg.out <- var.out[order(var.out$bio, decreasing=TRUE)[1:topn], ]
  return(rownames(hvg.out))
}


.bootstrap_clustering <- function(x, exprsMat, cellTypes, subject, verbose = FALSE){
  #get the expression matrix
  exprsMat <- exprsMat[,x]

  # filtering the genes with zero variance
  nonZeroVar <- apply(exprsMat, 1, stats::var)!=0


  exprsMat <- exprsMat[nonZeroVar,]

  hvg <- .scran_high_var(exprsMat)
  pca <- stats::prcomp(t(exprsMat[hvg,]), scale. = TRUE)$x[,1:10]

  num_G = length(unique(cellTypes))
  label = NULL

  while (is.null(label) | length(unique(label)) != length(unique(cellTypes))){
    label = NULL

    kmeans.result <- scClust::scClust(t(pca), num_G, similarity = "pearson", method = "kmeans", seed = 1, nstart = 100, iter.max = 1000)


    confusion <- table(cellTypes[x], kmeans.result$cluster)

    label <- rownames(confusion)[apply(confusion, 2, which.max)]
    # if(verbose){
    #   print(label)
    # }
    num_G <- num_G + 1
  }

  clusterRes <- label[kmeans.result$cluster]
  tab <- table(clusterRes, subject[x])
  # if(verbose){
  #   print(tab)
  # }
  tab_count <- c(tab)
  tab_prop <- c(tab/matrix(rep(colSums(tab), nrow(tab)), nrow = nrow(tab), byrow = T))
  return(list(count = tab_count, prop = tab_prop))
}



.CIbarColor <- function(n){
  if(n==1){
    colour <- "black"
  }else{
    colour <- c("red", "blue", "purple")
  }
  return(colour)
}
