library(xgboost)
library(SHAPforxgboost)
library(caret)
library(tidyverse)
library(ggplot2)
library(stringr)
library(dplyr)
library(gridExtra)
library(grid)
library(ggplotify)
library(data.table)
library(reshape2)

source("Source/MiDataProc.ML.Models.R")

# Extreme Gradient Boosting ----------------------

xgb.cla <- function(data, sam.dat.na, y.name, eta = 0.001, nrounds = 500, nfold = c(5, 10), alpha = 0, lambda = 1, stratified = TRUE,
                    loss.func = c("error", "auc", "logloss"), name){
  xgb.list <- list()  
  X = as.matrix(data[[name]])
  y = sam.dat.na[[y.name]]
  cat.names <- category.names(sam.dat.na, y.name)
  data.dmatrix <- xgb.DMatrix(data = X, label = y)
  
  ind <- 1
  xgb.cv.list <- list()
  output <- data.frame()
  
  for(currentMaxDepth in seq(4, 10, 2)){
    for(currentMCW in seq(2, 6, 2)){
      for(currentCbT in seq(0.5, 0.75, 1)){
        for(currentGamma in seq(0, 0.6, 0.2)){
          set.seed(578)
          xgboostModelCV <- xgb.cv(data =  data.dmatrix, nrounds = nrounds, nfold = nfold, showsd = TRUE, stratified = stratified, tree_method = "exact",
                                   metrics = loss.func, verbose = TRUE, objective = "binary:logistic", max_depth = currentMaxDepth, colsample_bytree = currentCbT,
                                   min_child_weight = currentMCW ,eta = eta, print_every_n = 10, booster = "gbtree", gamma = currentGamma, alpha = alpha, lambda = lambda,
                                   early_stopping_rounds = 150, nthread = 1)
          xvalidationScores <- as.data.frame(xgboostModelCV$evaluation_log)
          niter <- xgboostModelCV$best_iteration
          test <- as.numeric(tail(xvalidationScores[,4], 1))
          train <- as.numeric(tail(xvalidationScores[,2],1))

          xgb.cv.list[[ind]] <- xgboostModelCV
          ind <- ind + 1
          output <- rbind(output, c(niter, train, test, currentMaxDepth, currentMCW, currentCbT, currentGamma))
        }
      }
    }
  }
  varnames <- c("Optimal_trees", "Train", "Test", "max_depth", "min_child_weight", "colsample_bytree", "gamma")
  names(output) <- varnames
  
  if(loss.func == "auc"){
    cv.ind <- which.max(output$Test)
    cv.model <- xgb.cv.list[[cv.ind]]
    opt.tree <- output[cv.ind,][[1]]
    max_depth <- output[cv.ind,][[4]]
    min_child_weight <- output[cv.ind,][[5]]
    colsample_bytree <- output[cv.ind,][[6]]
    gamma <- output[cv.ind,][[7]]
  }
  else{
    cv.ind <- which.min(output$Test)
    cv.model <- xgb.cv.list[[cv.ind]]
    opt.tree <- output[cv.ind,][[1]]
    max_depth <- output[cv.ind,][[4]]
    min_child_weight <- output[cv.ind,][[5]]
    colsample_bytree <- output[cv.ind,][[6]]
    gamma <- output[cv.ind,][[7]]
  }
  
  set.seed(578)
  xgboostFit <- xgboost(data = X, label = y, nrounds = opt.tree, objective = "binary:logistic", tree_method = "exact", verbose = TRUE, eval_metric = loss.func,
                          max_depth = max_depth, colsample_bytree = colsample_bytree, eta = eta, min_child_weight = min_child_weight, print_every_n = 10,
                          booster = "gbtree", gamma = gamma, alpha = alpha, lambda = lambda, nthread = 1)
  
  xgb.list[["model"]] <- xgboostFit
  xgb.list[["cv.model"]] <- cv.model
  xgb.list[["loss"]] <- loss.func
  xgb.list[["tuning.result"]] <- output
  xgb.list[["data"]] <- list(x = X, y = y)
  return(xgb.list)
}

xgb.reg <- function(data, sam.dat.na, y.name, eta = 0.001, nrounds = 500, nfold = c(5, 10), alpha = 0, lambda = 1,
                    loss.func = c("huber", "rss"), name){
  xgb.list <- list()
  X = as.matrix(data[[name]])
  y = sam.dat.na[[y.name]]
  data.dmatrix <- xgb.DMatrix(data = X, label = y)
  
  if(loss.func == "huber"){
    # loss <- "reg:pseudohubererror"
    loss <- "reg:squarederror"
    eval.metric <- "mphe"
  }
  else if(loss.func == "rss"){
    loss <- "reg:squarederror"
    eval.metric <- "rmse"
  }
  
  ind <- 1
  xgb.cv.list <- list()
  output <- data.frame()
  
  for(currentMaxDepth in seq(4, 10, 2)){
    for(currentMCW in seq(2, 6, 2)){
      for(currentCbT in seq(0.5, 0.75, 1)){
        for(currentGamma in seq(0, 0.6, 0.2)){
          set.seed(578)
          xgboostModelCV <- xgb.cv(data = data.dmatrix, nrounds = nrounds, nfold = nfold, showsd = TRUE, tree_method = "exact",
                                   metrics = eval.metric, verbose = TRUE, objective = loss, max_depth = currentMaxDepth, colsample_bytree = currentCbT,
                                   min_child_weight = currentMCW ,eta = eta, print_every_n = 10, booster = "gbtree", gamma = currentGamma, alpha = alpha, lambda = lambda,
                                   early_stopping_rounds = 150, nthread = 1)
          
          xvalidationScores <- as.data.frame(xgboostModelCV$evaluation_log)
          niter <- xgboostModelCV$best_iteration
          test <- as.numeric(tail(xvalidationScores[,4], 1))
          train <- as.numeric(tail(xvalidationScores[,2],1))
          
          xgb.cv.list[[ind]] <- xgboostModelCV
          ind <- ind + 1
          output <- rbind(output, c(niter, train, test, currentMaxDepth, currentMCW, currentCbT, currentGamma))
        }
      }
    }
  }
  varnames <- c("Optimal_trees", "Train", "Test", "max_depth", "min_child_weight", "colsample_bytree", "gamma")
  names(output) <- varnames
  
  cv.ind <- which.min(output$Test)
  cv.model <- xgb.cv.list[[cv.ind]]
  opt.tree <- output[cv.ind,][[1]]
  max_depth <- output[cv.ind,][[4]]
  min_child_weight <- output[cv.ind,][[5]]
  colsample_bytree <- output[cv.ind,][[6]]
  gamma <- output[cv.ind,][[7]]
  
  set.seed(578)
  xgboostFit <- xgb.train(data = data.dmatrix, nrounds = opt.tree, objective = loss, tree_method = "exact", verbose = TRUE, eval_metric = eval.metric,
                          max_depth = max_depth, colsample_bytree = colsample_bytree, eta = eta, min_child_weight = min_child_weight, print_every_n = 10,
                          booster = "gbtree", gamma = gamma, alpha = alpha, lambda = lambda, nthread = 1)
  
  xgb.list[["model"]] <- xgboostFit
  xgb.list[["cv.model"]] <- cv.model
  xgb.list[["loss"]] <- loss.func
  xgb.list[["tuning.result"]] <- output
  xgb.list[["data"]] <- list(x = X, y = y)
  return(xgb.list)
}

xgb.mult <- function(data, sam.dat.na, y.name, eta = 0.001, nfold = c(5, 10), nrounds = 500, alpha = 0, lambda = 1, stratified = TRUE, name){
  xgb.list <- list()
  X = as.matrix(data[[name]])
  y = sam.dat.na[[y.name]]
  cat.names <- category.names(sam.dat.na, y.name)
  data.dmatrix <- xgb.DMatrix(data = X, label = y)
  num_class <- length(cat.names)
  ind <- 1
  xgb.cv.list <- list()
  output <- data.frame()
  
  # Tuning Process
  for(currentMaxDepth in seq(4, 10, 2)){
    for(currentMCW in seq(2, 6, 2)){
      for(currentCbT in seq(0.5, 0.75, 1)){
        for(currentGamma in seq(0, 0.6, 0.2)){
          set.seed(578) # metrics = loss.func, 
          xgboostModelCV <- xgb.cv(data =  data.dmatrix, nrounds = nrounds, nfold = nfold, showsd = TRUE, stratified = stratified, tree_method = "exact",
                                   num_class = num_class, verbose = TRUE, objective = "multi:softprob", max_depth = currentMaxDepth, colsample_bytree = currentCbT,
                                   min_child_weight = currentMCW ,eta = eta, print_every_n = 10, booster = "gbtree", gamma = currentGamma, alpha = alpha, lambda = lambda,
                                   early_stopping_rounds = 150, nthread = 1)
          xvalidationScores <- as.data.frame(xgboostModelCV$evaluation_log)
          niter <- xgboostModelCV$best_iteration
          test <- as.numeric(tail(xvalidationScores[,4], 1))
          train <- as.numeric(tail(xvalidationScores[,2],1))
          
          xgb.cv.list[[ind]] <- xgboostModelCV
          ind <- ind + 1
          output <- rbind(output, c(niter, train, test, currentMaxDepth, currentMCW, currentCbT, currentGamma))
        }
      }
    }
  }
  
  varnames <- c("Optimal_trees", "Train", "Test", "max_depth", "min_child_weight", "colsample_bytree", "gamma")
  names(output) <- varnames
  
  cv.ind <- which.min(output$Test)
  cv.model <- xgb.cv.list[[cv.ind]]
  opt.tree <- output[cv.ind,][[1]]
  max_depth <- output[cv.ind,][[4]]
  min_child_weight <- output[cv.ind,][[5]]
  colsample_bytree <- output[cv.ind,][[6]]
  gamma <- output[cv.ind,][[7]]
  
  set.seed(578) # , eval_metric = loss.func
  xgboostFit <- xgboost(data = X, label = y, nrounds = opt.tree, objective = "multi:softprob", tree_method = "exact", verbose = TRUE, num_class = num_class, 
                          max_depth = max_depth, colsample_bytree = colsample_bytree, eta = eta, min_child_weight = min_child_weight, print_every_n = 10,
                          booster = "gbtree", gamma = gamma, alpha = alpha, lambda = lambda, nthread = 1)
  
  xgb.list[["model"]] <- xgboostFit
  xgb.list[["cv.model"]] <- cv.model
  xgb.list[["tuning.result"]] <- output
  xgb.list[["data"]] <- list(x = X, y = y)

  return(xgb.list)
}

xgb.error.plot.2 <- function(xgb.list, rank.name){
  options(scipen = 999)
  eval.log <- xgb.list[[rank.name]][["cv.model"]]$evaluation_log
  std <- names(eval.log[,2]) %>% gsub("train_", "",.) %>% gsub("_mean", "",.)
  if(std == "mlogloss"){
    std <- "Cross entropy"
  }
  else if(std == "error"){
    std <- "Error rate"
  }
  else if(std == "auc"){
    std <- "AUC"
  }
  else if(std == "rmse"){
    std <- "Mean Squared Error"
    eval.log[,2] <- eval.log[,2]^2
    eval.log[,4] <- eval.log[,4]^2
  }
  else if(std == "mphe"){
    std <- "Mean Pseudo-Huber error"
  }
  data.frame(error = c(unlist(eval.log[,2]), unlist(eval.log[,4])),
             class = c(rep("Train Error", nrow(eval.log)),
                       rep("CV Error", nrow(eval.log))),
             nround = rep(1:nrow(eval.log), 2)) %>%
    ggplot(aes(nround, error, col = class)) +
    # geom_point(alpha = 0.5) +
    geom_line(linetype = "solid", linewidth = 1.2) +
    # geom_smooth(alpha = 0.1, se = FALSE, method = "gam") +
    theme_bw() +
    ggtitle("XGBoost Final Model",
            subtitle = sprintf("Level : %s  Iterations : %i", str_to_title(rank.name), dim(eval.log)[1])) +
    xlab("# Iterations") +
    ylab(std) +
    theme(axis.text = element_text(size = 11),
          axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14, vjust = +2),
          plot.title = element_text(size = 17),
          plot.subtitle = element_text(size = 15),
          legend.title = element_blank())
}

xgb.imp.list <- function(xgb.list, level.names){
  imp.list <- list()
  for(name in level.names){
    imp.list[[name]] <- xgb.importance(model = xgb.list[[name]]$model)
  }
  return(imp.list)
}

xgb.imp.list.ma <- function(xgb.list, level.names, data){
  xgb.importance <- xgb.imp.list(xgb.list, level.names)
  
  d <- subset(as.data.frame(xgb.importance[[level.names]]), select = c("Feature", "Gain"))
  ma <- colMeans(data[[level.names]])
  d$MeanAbundance <- rep(NA, nrow(d))
  
  for(v in names(ma)){
    if(v %in% d$Feature){
      ind <- which(v == d$Feature)
      d$MeanAbundance[ind] <- ma[[v]]
    }
  }
  return(d)
}

xgb.imp.plot <- function(xgb.list, level.names, data, data.type, n = 30, is.cat = TRUE){
  imp.df <- xgb.imp.list.ma(xgb.list, level.names, data)
  
  if(data.type == "clr"){
    data.type <- "CLR"
  }
  else if(data.type == "prop"){
    data.type <- "Proportion"
  }
  else if(data.type == "rare.count"){
    data.type <- "Rarefied Count"
  }
  else if(data.type == "arcsin"){
    data.type <- "Arcsine-Root"
  }
  
  yl <- ifelse(is.cat, "Decrease in Gini Impurity", "Decrease in Mean Squared Error")
  
  imp.df <- imp.df %>% 
    arrange(desc(Gain)) %>%
    top_n(n, Gain)
  
  imp.df <- imp.df %>%
    ggplot(aes(x=reorder(Feature, Gain), y=Gain)) +
    geom_segment( aes(x=reorder(Feature, Gain), xend=reorder(Feature, Gain), y=0, yend=Gain), color="grey") +
    geom_point(aes(color = MeanAbundance), size=5, alpha=0.8) +
    scale_color_gradient(low = "blue", high = "red", name = sprintf("Mean\nAbundance\n(%s)", data.type)) +
    theme_light() +
    coord_flip() +
    theme(
      panel.grid.major.y = element_blank(),
      panel.border = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.y = element_text(size = 13),
      axis.title.x = element_text(size = 13)
    ) +
    xlab(element_blank()) +
    ylab(yl)

  imp.df
}

xgb.shap.summary <- function(xgb.list, rank.name, n = 20){
  X <- data.matrix(xgb.list[[rank.name]]$data$x)
  shap <- shap.prep(xgb.list[[rank.name]]$model, X_train = X, top_n = n)
  shap.plot.summary(shap, scientific = FALSE) +
    theme(
      axis.line.y = element_blank(), 
      axis.ticks.y = element_blank(), 
      legend.position = "bottom", 
      legend.title = element_text(size = 11), 
      legend.text = element_text(size = 9),
      axis.title.x = element_text(size = 13),
      axis.text.y = element_text(size = 13)
    )
}

xgb.shap.summary.2 <- function (data, top_n = 20, model) {
  xgb.ggplot.shap.summary(data = data, top_n = top_n, model = model) +
    theme_light() +
    xlab(element_blank()) +
    ylab("SHAP value (impact on model output)") +
    labs(colour = "Feature\nvalue") +
    theme(
      axis.line.y = element_blank(), 
      axis.ticks.y = element_blank(), 
      legend.title = element_text(size = 11), 
      legend.text = element_text(size = 9),
      axis.title.x = element_text(size = 13),
      axis.text.y = element_text(size = 13)
    )
}



xgb.shap.prep <- function(xgb.list, rank.name, n = 20){
  X <- data.matrix(xgb.list[[rank.name]]$data$x)
  shap <- shap.prep(xgb.list[[rank.name]]$model, X_train = X, top_n = n)
  return(shap)
}

xgb.shap.dependence <- function(xgb.list, rank.name, n = 20){
  shap <- xgb.shap.prep(xgb.list, rank.name, n = n)
  
  plot.list <- list()
  for (v in shap.importance(shap, names_only = TRUE)) {
    p <- shap.plot.dependence(shap, x = v, color_feature = v, alpha = 0.8) +
      labs(color = "") +
      theme(
        plot.title = element_text(size = 10),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.margin = grid::unit(c(0, 0.75, 0.2, 0.75), unit = "cm")
      ) +
      ggtitle(v)
    plot.list[[v]] <- p
  }
  grid.arrange(grobs = plot.list, nrow = 5, left = textGrob("SHAP", rot = 90), bottom = textGrob("Feature Value", rot = 0), as.table = FALSE)
}

xgb.imp.var <- function(xgb.importance, rank.name, n = 10){
  imp.var <- (xgb.importance[[rank.name]] %>% top_n(n, Gain))$Feature
  return(imp.var)
}

xgb.shap.data.ext <- function (data, shap_contrib = NULL, features = NULL, top_n = 1, 
                           model = NULL, trees = NULL, target_class = NULL, approxcontrib = FALSE, 
                           subsample = NULL, max_observations = 100000) {
  if (!is.matrix(data) && !inherits(data, "dgCMatrix")) 
    stop("data: must be either matrix or dgCMatrix")
  if (is.null(shap_contrib) && (is.null(model) || !inherits(model, 
                                                            "xgb.Booster"))) 
    stop("when shap_contrib is not provided, one must provide an xgb.Booster model")
  if (is.null(features) && (is.null(model) || !inherits(model, 
                                                        "xgb.Booster"))) 
    stop("when features are not provided, one must provide an xgb.Booster model to rank the features")
  if (!is.null(shap_contrib) && (!is.matrix(shap_contrib) || 
                                 nrow(shap_contrib) != nrow(data) || ncol(shap_contrib) != 
                                 ncol(data) + 1)) 
    stop("shap_contrib is not compatible with the provided data")
  if (is.character(features) && is.null(colnames(data))) 
    stop("either provide `data` with column names or provide `features` as column indices")
  if (is.null(model$feature_names) && model$nfeatures != ncol(data)) 
    stop("if model has no feature_names, columns in `data` must match features in model")
  if (!is.null(subsample)) {
    idx <- sample(x = seq_len(nrow(data)), size = as.integer(subsample * 
                                                               nrow(data)), replace = FALSE)
  }
  else {
    idx <- seq_len(min(nrow(data), max_observations))
  }
  data <- data[idx, ]
  if (is.null(colnames(data))) {
    colnames(data) <- paste0("X", seq_len(ncol(data)))
  }
  if (!is.null(shap_contrib)) {
    if (is.list(shap_contrib)) {
      shap_contrib <- if (!is.null(target_class)) 
        shap_contrib[[target_class + 1]]
      else Reduce("+", lapply(shap_contrib, abs))
    }
    shap_contrib <- shap_contrib[idx, ]
    if (is.null(colnames(shap_contrib))) {
      colnames(shap_contrib) <- paste0("X", seq_len(ncol(data)))
    }
  }
  else {
    shap_contrib <- predict(model, newdata = data, predcontrib = TRUE, 
                            approxcontrib = approxcontrib)
    if (is.list(shap_contrib)) {
      shap_contrib <- if (!is.null(target_class)) 
        shap_contrib[[target_class + 1]]
      else Reduce("+", lapply(shap_contrib, abs))
    }
  }
  if (is.null(features)) {
    if (!is.null(model$feature_names)) {
      imp <- xgb.importance(model = model, trees = trees)
    }
    else {
      imp <- xgb.importance(model = model, trees = trees, 
                            feature_names = colnames(data))
    }
    top_n <- top_n[1]
    if (top_n < 1 | top_n > 100) 
      stop("top_n: must be an integer within [1, 100]")
    features <- imp$Feature[1:min(top_n, NROW(imp))]
  }
  if (is.character(features)) {
    features <- match(features, colnames(data))
  }
  shap_contrib <- shap_contrib[, features, drop = FALSE]
  data <- data[, features, drop = FALSE]
  list(data = data, shap_contrib = shap_contrib)
}

xgb.shap.imp.var <- function(xgb.list, rank.name, colnames.list, n = 20){
  X <- xgb.list[[rank.name]]$data$x
  fit <- xgb.list[[rank.name]]$model
  xgb.importance <- xgb.shap.data.ext(data = xgb.list[[rank.name]]$data$x, top_n = n, model = xgb.list[[rank.name]]$model)$shap_contrib %>% colnames
  if(n > ncol(xgb.list[[rank.name]]$data$x)){
    n <- ncol(xgb.list[[rank.name]]$data$x)
  }
  var <- xgb.importance[1:n]
  nm <- colnames.df(colnames.list, rank.name)
  ind <- integer()
  for(v in var){
    ind <- c(ind, which(v == rownames(nm)))
  }
  # new.ind <- sort(ind)
  nm <- nm[ind,, drop = FALSE]
  return(nm)
}

xgb.pdp.mult <- function(xgb.list, n, name, data.type, level.names, label){
  X <- xgb.list[[name]]$data$x
  fit <- xgb.list[[name]]$model
  xgb.importance <- xgb.shap.data.ext(data = xgb.list[[name]]$data$x, top_n = n, model = xgb.list[[name]]$model)$shap_contrib %>% colnames
  n <- min(length(xgb.importance), n)
  feature <- xgb.importance[1:n]
  result <- data.frame()
  
  if(data.type == "clr"){
    type = "CLR"
  }
  else if(data.type == "prop"){
    type = "Proportion"
  }
  else if(data.type == "rare.count"){
    type = "Rarefied Count"
  }
  else if(data.type == "arcsin"){
    type = "Arcsine-Root"
  }
  
  for(taxon.name in feature){
    val <- numeric()
    prob <- list()
    label.len <- length(label)
    for(i in 1:label.len){
      prob[[i]] <- numeric()
    }
    # p1 <- numeric()
    # p2 <- numeric()
    # p3 <- numeric()
    taxon.val <- seq(min(X[,taxon.name]), max(X[,taxon.name]),len = 100)
    
    for(i in 1:length(taxon.val)){
      newX <- X
      newX[,taxon.name] <- rep(taxon.val[i], nrow(newX))
      y_pred_prob <- predict(fit, newX, reshape = TRUE)
      y_pred_prob_mean <- apply(y_pred_prob, 2, mean)
      val <- c(val, taxon.val[i])
      for(j in 1:label.len){
        prob[[j]] <- c(prob[[j]], y_pred_prob_mean[j])
      }
      # p1 <- c(p1, y_pred_prob_mean[1])
      # p2 <- c(p2, y_pred_prob_mean[2])
      # p3 <- c(p3, y_pred_prob_mean[3])
    }
    prob.result <- data.frame(prob)
    colnames(prob.result) <- paste0("p", 1:label.len)
    prob.result <- data.frame(val, prob.result) %>%
      reshape2::melt(id.vars = "val", variable.name = "Category", value.name = "Prob")
    prob.result$title <- rep(taxon.name, nrow(prob.result))
    result <- rbind(result, prob.result)
  }
  result$title <- factor(result$title, levels = feature)
  
  p <- ggplot(result, aes(val, Prob)) + 
    geom_line(aes(color = Category), size = 0.8) +
    theme_light() +
    xlab(type) + 
    ylab("Predicted Value") +
    scale_color_discrete(name = "Category", labels = label) +
    theme(
      axis.text.x = element_text(size = 7.5),
      axis.text.y = element_text(size = 7.5),
      strip.text = element_text(size=12),
      panel.spacing.x = unit(1, "lines"),
      legend.title = element_blank()
    ) +
    facet_wrap(~ title, scales = "free_x", nrow = 5, dir = "v")
  p
}

produce.new.levels <- function(level){
  new.levels <- try(as.numeric(level), silent = TRUE)
  if(NA %in% new.levels){
    return(1:length(level)-1)
  } 
  else{
    if(is.factor(level)){
      return(1:length(level)-1)
    }
    else{
      return(level)
    }
  }
}