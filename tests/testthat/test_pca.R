context("test my own function with the R package function")

test_that("Whether PCA function give us the reasonable output",{

  i.df <- iris[, 2:3]

  pc.fit <- prcomp(i.df)
  PC1 <- pc.fit[["rotation"]][,1]
  PC1.mat <- matrix(PC1, nrow=nrow(iris), ncol=2, byrow=TRUE)
  mean.vec <- colMeans(i.df)
  mean.mat <- matrix(mean.vec, nrow=nrow(iris), ncol=2, byrow=TRUE)
  pred.mat <- mean.mat + PC1.mat * pc.fit[["x"]][, 1]
  pc.error = sum( (pred.mat - i.df)^2)/nrow(iris)
  pc.error
  my.pc.fit <-pca(i.df)
  my.error = sum(my.pc.fit$cumVar)/nrow(iris)
  my.error
  expect_less_than(abs(pc.error - my.error),1.5)
})
