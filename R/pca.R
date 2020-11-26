
#' Title
#'
#' @param dt  data points
#'
#' @return {PCA=principal, cumVar=cumVar }
#' @import stats
#' @export

pca <- function(dt){
  #center each column
  sca.dt <- scale(dt, center=TRUE, scale=FALSE)
  #calculate the Covariance matrix
  cov.sca.dt <- cov(sca.dt)
  #Singular value decomposition
  svd.dt <- svd(cov.sca.dt)

  svd.value <- svd.dt$d
  svd.vector <-svd.dt$u

  order_value <- order(svd.value,decreasing = T)
  values <- svd.value[order_value]
  valueSum <- sum(values)
  cumVar <- cumsum(values)/valueSum * 100

  order_vector <- svd.vector[,order_value]

  principal <- sca.dt %*% order_vector
  return(list(PCA=principal, cumVar=cumVar ))
}
