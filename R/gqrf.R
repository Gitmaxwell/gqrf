#' @title Geographic quantile regression forest
#' @description Predict values at unknown locations using geographic quantile regression forest
#' @param formula Formula in the format y~X where y is the target variable and X are predictor/s
#' @param data Training data of class SpatialPointsDataFrame (sp)
#' @param newdata Data of class SpatialPointsDataFrame (sp)comprising of new values to be estimated
#' @param k Number of nearest neighbors to use. The default is the number of rows in data minus 1
#' @param ntree Number of trees.
#' @param mtry Number of variables to possibly split at in each node. Default is the (rounded down) square root of the number variables.
#' @param quantiles Vector of the upper, lower and middle quantiles. Default is c(0.05,0.5,0.95)
#' @details
#' Uses geographic quantile regression forest (Maxwell, 2020) for spatially interpolating values at new, unknown locations.
#' @export
#' @import data.table
#' @return A data frame consisting of x (x coordinate), y (y coordinate), target (the target variable), covars (covariates/auxilary data), pred (the predicted value) and var (the variance at the chosen quantiles)
#' @references
#' * Maxwell, K., Rajabi, M., Esterle, J. (2020). Spatial interpolation of coal geochemical properties using geographic quantile regression forest. Manuscript submitted for publication.
#'
#' * Georganos, S., Grippa, T., Niang Gadiaga, A., Linard, C., Lennert, M., Vanhuysse, S., ... & Kalogirou, S. (2019). Geographical random forests: a spatial extension of the random forest algorithm to address spatial heterogeneity in remote sensing and population modelling. Geocarto International, 1-16
#' @examples
#' #example

gqrf = function(formula,data,newdata,k=NULL,ntree=500,mtry=NULL,quantiles = c(0.05,0.500,0.95))

{
  f <- formula(formula)
  covars <- attr(stats::terms(f),"term.labels")
  target = all.vars(f)[1]
  ntree <- ntree
  if (is.null(k)) {k =nrow(data)-1}
  if (is.null(mtry)) {mtry= max(floor(length(covars)/3), 1)}
  get_nnn_rf = function(x) {
    train = data
    test = newdata[x,]
    nn1 <- nabor::knn(sp::coordinates(train),sp::coordinates(test),k)
    traindt = data.table::setDT(as.data.frame(train))
    trainnn = traindt[nn1$nn.idx[1,]]
    trainnn$D = nn1$nn.dists[1,]
    range02 <- function(x){(x-max(x))/(min(x)-max(x))} #scales 1 to 0
    trainnn$Wts = range02(trainnn$D)
    quantiles = quantiles
    rl <- ranger::ranger(formula=f, data = trainnn, importance = 'none', case.weights = trainnn$Wts, oob.error = FALSE,quantreg=TRUE, num.trees = ntree)
    test = data.table::setDT(as.data.frame(test))
    pred = stats::predict(rl, data = test,type="quantiles",quantiles=quantiles)$predictions
    pred_df = data.frame(Predicted=pred[,2], var = (pred[,3]-pred[,1])/2)
    test$pred = pred_df$Predicted
    test$var = pred_df$var
    test$target = target
    test$covars = paste(covars,collapse = " , ")
    test = test[,c("x", "y","target","covars","pred","var")]
  }
  l=  lapply(1:nrow(newdata),  get_nnn_rf)
  r =data.table::rbindlist(l)
  return(r)
}
