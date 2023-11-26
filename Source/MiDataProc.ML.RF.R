library(randomForest)
library(caret)
library(tidyverse)
library(ggplot2)
library(stringr)
library(edarf)
library(data.table)
library(dplyr)
library(reshape2)

source("Source/MiDataProc.ML.Models.R")

# Random Forest ----------------------

get.opt.pred <- function(rfcv){
  error.cv <- rfcv$error.cv
  return(min(as.numeric(names(which(error.cv == min(error.cv))))))
}

rf.cla.rev <- function(data, sam.dat.na, y.name, nfold = c(5, 10), ntree, stratified = TRUE, name, p = 0.75){
  rf.list <- list()
  
  X = data[[name]]
  y = as.factor(sam.dat.na[[y.name]])
  
  tr.ind <- y %>% createDataPartition(p = p, list = FALSE)

  # Train / Test Split
  train_X = X[tr.ind,]
  train_Y = y[tr.ind]

  test_X = X[-tr.ind,]
  test_Y = y[-tr.ind]
  
  set.seed(578)
  rfcv <- rfcv(train_X, train_Y, cv.fold=nfold, scale="log", step=0.8, recursive=FALSE, ntree = ntree)
  (mtry <- get.opt.pred(rfcv))
  
  train_X <- remove.symb(train_X)
  new.dat <- cbind(train_X, train_Y)
  colnames(new.dat)[dim(new.dat)[2]] <- y.name

  f1 <<- as.formula(paste(y.name, "~", ".", sep=" "))

  if(stratified){
    set.seed(578)
    fit <- randomForest(f1, data = new.dat, mtry = mtry, ntree = ntree, importance = TRUE, strata = new.dat[[y.name]])
  }
  else{
    set.seed(578)
    fit <- randomForest(f1, data = new.dat, mtry = mtry, ntree = ntree, importance = TRUE)
  }
  
  rf.list[["rfcv"]] <- rfcv
  rf.list[["fit"]] <- fit
  if(p != 1) rf.list[["predicted"]] <- predict(fit, test_X)
  rf.list[["data"]] <- list(x = X, y = y)
  rf.list[["train"]] <- list(x = train_X, y = train_Y)
  if(p != 1) rf.list[["test"]] <- list(x = test_X, y = test_Y)
  return(rf.list)
}

rf.reg.rev <- function(data, sam.dat.na, y.name, nfold = c(5, 10), ntree, name, p = 0.75){
  rf.list <- list()
  
  X = data[[name]]
  y = sam.dat.na[[y.name]]
  
  set.seed(487)
  tr.ind <- y %>% createDataPartition(p = p, list = FALSE)

  # Train / Test Split
  train_X = X[tr.ind,]
  train_Y = y[tr.ind]

  test_X = X[-tr.ind,]
  test_Y = y[-tr.ind]
  
  set.seed(578)
  rfcv <- rfcv(train_X, train_Y, cv.fold=nfold, scale="log", step=0.8, recursive=FALSE, ntree = ntree)
  (mtry <- get.opt.pred(rfcv))
  
  train_X <- remove.symb(train_X)
  new.dat <- cbind(train_X, train_Y)
  colnames(new.dat)[dim(new.dat)[2]] <- y.name
  
  f1 <<- as.formula(paste(y.name, "~", ".", sep=" "))
  set.seed(578)
  fit <- randomForest(f1, data = new.dat, mtry = mtry, ntree = ntree, importance = TRUE)
  
  rf.list[["rfcv"]] <- rfcv
  rf.list[["fit"]] <- fit
  if(p != 1) rf.list[["predicted"]] <- predict(fit, test_X)
  rf.list[["data"]] <- list(x = X, y = y)
  rf.list[["train"]] <- list(x = train_X, y = train_Y)
  if(p != 1) rf.list[["test"]] <- list(x = test_X, y = test_Y)
  return(rf.list)
}

cv.mtry <- function(rf.list, rank.name, is.cat = TRUE){
  ylab <- ifelse(is.cat, "CV Error", "CV MSE")
  
  cvdf <- data.frame(rf.list[[rank.name]]$rfcv$error.cv)
  cvdf <- cbind(cvdf, as.numeric(rownames(cvdf)))
  colnames(cvdf) <- c("cv.error", "mtry")
  
  p <- ggplot(cvdf, aes(x = mtry, y = cv.error)) +
    geom_line(linetype = "solid", size = 1)+
    geom_point() +
    labs(x = "# Randomly Selected Taxa", y = ylab) +
    theme_bw()
  p + theme(
    axis.text.x = element_text(size = 13),
    axis.text.y = element_text(size = 13),
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15, vjust = +2))
}

error.plot <- function(rf.list, rank.name, is.cat = TRUE){
  if(is.cat){
    train.error <- rf.list[[rank.name]]$fit$err.rate[,1]
    err1 <- tail(rf.list[[rank.name]]$fit$err.rate[,1], 1)
    err2 <- 1 - mean(rf.list[[rank.name]]$predicted == rf.list[[rank.name]]$test$y)
    # sub <- sprintf("Train Error: %0.4f", err1)
  }
  else{
    train.error <- rf.list[[rank.name]]$fit$mse
    err1 <- tail(rf.list[[rank.name]]$fit$err.rate[,1], 1)
    err2 <- mean((rf.list[[rank.name]]$predicted - rf.list[[rank.name]]$test$y)^2)
    # sub <- sprintf("Train MSE: %0.4f", err1)
  }
  
  train.error.df <- data.frame(error = train.error, ntree = 1:length(train.error), Group = rep("Train", length(train.error)))
  # test.error.df <- data.frame(error = test.error, ntree = 1:length(test.error), Group = rep("Test", length(test.error)))
  
  # error.df <- rbind(train.error.df, test.error.df)
  
  p <- ggplot(train.error.df, aes(x =ntree, y = error, col = Group)) +
    geom_line(linetype = "solid", size = 0.8) +
    labs(x = "# Trees", y = "OOB Error") +
    theme_bw()
  p + theme(
    axis.text.x = element_text(size = 13),
    axis.text.y = element_text(size = 13),
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    legend.text = element_text(size = 12),
    legend.title = element_blank(),
    legend.position = "none")
}

rf.imp.df <- function(rf.list, rank.name, type){
  fit <- rf.list[[rank.name]]$fit
  if(type == 0){
    imp <- randomForest::importance(fit)
  }
  else if(type == 1){
    imp <- randomForest::importance(fit, type = 1)
  }
  else if(type == 2){
    imp <- randomForest::importance(fit, type = 2)
  }
  if("%IncMSE" %in% colnames(imp)){
    ind <- which(colnames(imp) == "%IncMSE")
    colnames(imp)[ind] <- "IncMSE"
  }
  return(as.data.frame(imp))
}

rf.imp.df.order <- function(rf.list, rank.name, n = 10, is.cat = TRUE, colnames.list){
  if(is.cat){
    imp.df <- rf.imp.df(rf.list, rank.name, type = 2)
    new.names <- colnames.list$origin[[rank.name]]
    rownames(imp.df) <- new.names
    
    imp.df <- imp.df %>% 
      mutate(names = rownames(imp.df)) %>% 
      arrange(desc(MeanDecreaseGini)) %>%
      top_n(n, MeanDecreaseGini) 
  }
  else{
    imp.df <- rf.imp.df(rf.list, rank.name, type = 1)
    new.names <- colnames.list$origin[[rank.name]]
    rownames(imp.df) <- new.names
    
    imp.df <- imp.df %>% 
      mutate(names = rownames(imp.df)) %>% 
      arrange(desc(IncMSE)) %>%
      top_n(n, IncMSE) 
  }
  
  return(imp.df)
}

rf.imp.plot <- function(rf.list, rank.name, type, n = 30, is.cat = TRUE, data, data.type){
  
  imp.df1 <- rf.imp.df(rf.list, rank.name, type = 1)
  imp.df2 <- rf.imp.df(rf.list, rank.name, type = 2)
  
  ma <- get.mean.abundance(data, rank.name)
  imp.df1 <- data.frame(imp.df1, ma)
  colnames(imp.df1)[2] <- "MeanAbundance"
  imp.df2 <- data.frame(imp.df2, ma)
  colnames(imp.df2)[2] <- "MeanAbundance"
  
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
  
  if(is.cat){
    imp.df1 <- imp.df1 %>% 
      mutate(names = rownames(imp.df1)) %>% 
      arrange(desc(MeanDecreaseAccuracy)) %>%
      top_n(n, MeanDecreaseAccuracy)
    
    imp.df1 <- imp.df1 %>% 
      ggplot(aes(x=reorder(names, MeanDecreaseAccuracy), y=MeanDecreaseAccuracy)) +
      geom_segment( aes(x=reorder(names, MeanDecreaseAccuracy), xend=reorder(names, MeanDecreaseAccuracy), y=0, yend=MeanDecreaseAccuracy), color="grey") +
      # geom_point( color=ifelse(imp.df1$MeanDecreaseAccuracy > 0, "#023020", "#8B0000"), size=5, alpha=0.8) +
      geom_point(aes(color = MeanAbundance), size=5, alpha=0.8) +
      scale_color_gradient(low = "blue", high = "red", name = sprintf("Mean\nAbundance\n(%s)", data.type)) +
      theme_light() +
      coord_flip() +
      theme(
        panel.grid.major.y = element_blank(),
        panel.border = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 12)
      ) +
      xlab(element_blank()) +
      ylab("Decrease in Accuracy")
    
    imp.df2 <- imp.df2 %>% 
      mutate(names = rownames(imp.df2)) %>% 
      arrange(desc(MeanDecreaseGini)) %>%
      top_n(n, MeanDecreaseGini) 
    
    imp.df2 <- imp.df2 %>% 
      ggplot(aes(x=reorder(names, MeanDecreaseGini), y=MeanDecreaseGini)) +
      geom_segment( aes(x=reorder(names, MeanDecreaseGini), xend=reorder(names, MeanDecreaseGini), y=0, yend=MeanDecreaseGini), color="grey") +
      # geom_point( color=ifelse(imp.df2$MeanDecreaseGini > 0, "#023020", "#8B0000"), size=5, alpha=0.8) +
      geom_point(aes(color = MeanAbundance), size=5, alpha=0.8) +
      scale_color_gradient(low = "blue", high = "red", name = sprintf("Mean\nAbundance\n(%s)", data.type)) +
      theme_light() +
      coord_flip() +
      theme(
        panel.grid.major.y = element_blank(),
        panel.border = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 12)
      ) +
      xlab(element_blank()) +
      ylab("Decrease in Gini Impurity")
  }
  else{
    names(imp.df1)[1] <- "IncMSE"
    imp.df1 <- imp.df1 %>% 
      mutate(names = rownames(imp.df1)) %>% 
      arrange(desc(IncMSE)) %>%
      top_n(n, IncMSE) 
    
    imp.df1 <- imp.df1 %>% 
      ggplot(aes(x=reorder(names, IncMSE), y=IncMSE)) +
      geom_segment( aes(x=reorder(names, IncMSE), xend=reorder(names, IncMSE), y=0, yend=IncMSE), color="grey") +
      # geom_point( color=ifelse(imp.df1$IncMSE > 0, "#023020", "#8B0000"), size=5, alpha=0.8) +
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
      ylab("Decrease in Mean Squared Error")
    
    imp.df2 <- imp.df2 %>% 
      mutate(names = rownames(imp.df2)) %>% 
      arrange(desc(IncNodePurity)) %>%
    top_n(n, IncNodePurity)
    
    imp.df2 <- imp.df2 %>% 
      ggplot(aes(x=reorder(names, IncNodePurity), y=IncNodePurity)) +
      geom_segment( aes(x=reorder(names, IncNodePurity), xend=reorder(names, IncNodePurity), y=0, yend=IncNodePurity), color="grey") +
      # geom_point( color=ifelse(imp.df2$IncNodePurity > 0, "#023020", "#8B0000"), size=5, alpha=0.8) +
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
      ylab("Decrease in Node Impurity")
  }
  if(type == 0){
    return(imp.df1 + imp.df2)
  }
  else if(type == 1){
    return(imp.df1)
  }
  else if(type == 2){
    return(imp.df2)
  }
}

pd.plot <- function(rf.list, rank.name){
  pd <- partial_dependence(fit = rf.list[[rank.name]]$fit, vars = colnames(rf.list[[rank.name]]$data$x), data = cbind(rf.list[[rank.name]]$data$x, rf.list[[rank.name]]$data$y))
  plot_pd(pd)
}

rf.imp.plot.order <- function(rf.list, rank.name, is.cat = TRUE){
  imp_df <- rf.imp.df(rf.list, rank.name, type = 0)
  if(is.cat){
    imp_df <- imp_df %>% 
      mutate(names = rownames(imp_df)) %>% 
      arrange(desc(MeanDecreaseGini))
  }
  else{
    imp_df <- imp_df %>% 
      mutate(names = rownames(imp_df)) %>% 
      arrange(desc(IncMSE))
  }
  nm <- imp_df$names
  return(list(imp_df = imp_df, nm = nm))
}

rf.pd.plot <- function (rf.list, rank.name, facet = NULL, is.cat = TRUE, n = 10, data.type){
  imp <- rf.imp.plot.order(rf.list, rank.name, is.cat = is.cat)
  imp_df <- imp$imp_df
  if(n > ncol(rf.list[[rank.name]]$train$x)){
    n <- ncol(rf.list[[rank.name]]$train$x)
  }
  nm <- imp$nm[1:n]
  
  type <- character()
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
  
  pd <- partial_dependence(fit = rf.list[[rank.name]]$fit, vars = nm, data = cbind(rf.list[[rank.name]]$train$x, rf.list[[rank.name]]$train$y))
  atts = attributes(pd)
  if (length(atts$vars) > 2 & atts$interaction) {
    stop("too many variables to plot.")
  }
  if (!is.null(facet)) {
    if (!any(facet %in% atts$vars)) 
      stop("facet must be one of the variables in the pd argument.")
    ordering = unique(sort(pd[[facet]]))
    if (is.factor(pd[[facet]])) 
      labels = paste0(facet, " = ", as.character(ordering))
    else labels = paste0(facet, " = ", as.character(signif(ordering, 
                                                           3)))
    pd[[facet]] = factor(pd[[facet]], levels = ordering, 
                         labels = labels)
    vars = atts$vars[atts$vars != facet]
  }
  else {
    vars = atts$vars
  }
  if (!is.null(facet) | !atts$interaction) {
    dat = data.table::melt(pd, id.vars = c(atts$target, facet), na.rm = TRUE)
    if (is.character(dat$value)) {
      dat$value = as.numeric(dat$value)
    }
    if (length(atts$target) > 1) 
      dat = data.table::melt(dat, id.vars = c("variable", "value", 
                                  facet), value.name = "prediction", variable.name = "class", 
                 na.rm = TRUE)
    if (length(atts$target) == 1) {
      p = ggplot(dat, aes_string("value", atts$target))
    }
    else {
      p = ggplot(dat, aes_string("value", "prediction", 
                                 colour = "class"))
    }
    p = p +
      geom_line() +
      geom_point() + 
      theme_light()
  }
  else {
    dat = pd
    dat = data.table::melt(dat, id.vars = vars)
    p = ggplot(dat, aes_string(vars[1], vars[2], fill = "value")) + 
      geom_raster()
    facet = "variable"
  }
  if (length(vars) == 1) 
    p = p + labs(x = vars)
  if (length(atts$vars) > 1) {
    if (!atts$interaction) {
      p = p +
        facet_wrap(~variable, scales = "free_x", nrow = 5, dir = "v") +
        xlab(type) + ylab("Predicted Value") +
        theme(
          axis.text.x = element_text(size = 7.5),
          axis.text.y = element_text(size = 7.5),
          strip.text = element_text(size=12),
          panel.spacing.x = unit(1, "lines"),
          legend.title = element_blank()
        )
        # theme_light()
    }
    else {
      p = p +
        facet_wrap(as.formula(paste0("~", facet)), nrow = 5, dir = "v") +
        xlab(type) + ylab("Predicted Value") +
        theme(
          axis.text.x = element_text(size = 7.5),
          axis.text.y = element_text(size = 7.5),
          strip.text = element_text(size=12),
          panel.spacing.x = unit(1, "lines"),
          legend.title = element_blank()
        )
    }
  }
  p
}

rf.pd.var.used <- function(rf.list, rank.name, colnames.list, n = 10, is.cat = TRUE){
  imp <- rf.imp.plot.order(rf.list, rank.name, is.cat = is.cat)
  imp_df <- imp$imp_df
  if(n > ncol(rf.list[[rank.name]]$data$x)){
    n <- ncol(rf.list[[rank.name]]$data$x)
  }
  var <- imp$nm[1:n]
  nm <- colnames.df(colnames.list, rank.name)
  ind <- integer()
  for(v in var){
    ind <- c(ind, which(v == rownames(nm)))
  }
  # new.ind <- sort(ind)
  nm <- nm[ind,, drop = FALSE]
  return(nm)
}

# rf.pdp.bin <- function(rf.list, n, name, data.type, level.names, label){
#   X <- rf.list[[name]]$data$x
#   fit <- rf.list[[name]]$fit
#   rf.importance <- rf.imp.df(rf.list, name, type = 2) %>%
#     mutate(taxon = rownames(rf.imp.df(rf.list, name, type = 2))) %>%
#     arrange(-MeanDecreaseGini) %>%
#     rownames
#   n <- min(length(rf.importance), n)
#   feature <- rf.importance[1:n]
#   result <- data.frame()
#   
#   if(data.type == "clr"){
#     type = "CLR"
#   }
#   else if(data.type == "prop"){
#     type = "Proportion"
#   }
#   else if(data.type == "rare.count"){
#     type = "Rarefied Count"
#   }
#   else if(data.type == "arcsin"){
#     type = "Arcsine-Root"
#   }
#   
#   for(taxon.name in feature){
#     val <- numeric()
#     p1 <- numeric()
#     p2 <- numeric()
#     taxon.val <- seq(min(X[,taxon.name]), max(X[,taxon.name]),len = 100)
#     
#     for(i in 1:length(taxon.val)){
#       newX <- X
#       newX[,taxon.name] <- rep(taxon.val[i], nrow(newX))
#       y_pred_prob <- predict(fit, newX, type = "prob")
#       y_pred_prob_mean <- apply(y_pred_prob, 2, mean)
#       val <- c(val, taxon.val[i])
#       p1 <- c(p1, y_pred_prob_mean[1])
#       p2 <- c(p2, y_pred_prob_mean[2])
#     }
#     prob.result <- data.frame(val, p1, p2) %>%
#       reshape2::melt(id.vars = "val", variable.name = "Category", value.name = "Prob")
#     prob.result$title <- rep(taxon.name, nrow(prob.result))
#     result <- rbind(result, prob.result)
#   }
#   result$title <- factor(result$title, levels = feature)
#   
#   p <- ggplot(result, aes(val, Prob)) + 
#     geom_line(aes(color = Category), size = 0.8) +
#     theme_light() +
#     xlab(type) + 
#     ylab("Predicted Value") +
#     scale_color_discrete(name = "Category", labels = c("0", "1", "2")) +
#     theme(
#       axis.text.x = element_text(size = 7.5),
#       axis.text.y = element_text(size = 7.5),
#       strip.text = element_text(size=12),
#       panel.spacing.x = unit(1, "lines"),
#       legend.title = element_blank()
#     ) +
#     facet_wrap(~ title, scales = "free_x", nrow = 5, dir = "v")
#   p
# }
# 
# rf.pdp.reg <- function(rf.list, n, name, data.type, level.names){
#   X <- rf.list[[name]]$data$x
#   fit <- rf.list[[name]]$fit
#   rf.importance <- rf.imp.df(rf.list, name, type = 1) %>%
#     mutate(taxon = rownames(rf.imp.df(rf.list, name, type = 1))) %>%
#     arrange(-IncMSE) %>%
#     rownames
#   n <- min(length(rf.importance), n)
#   feature <- rf.importance[1:n]
#   result <- data.frame()
#   
#   if(data.type == "clr"){
#     type = "CLR"
#   }
#   else if(data.type == "prop"){
#     type = "Proportion"
#   }
#   else if(data.type == "rare.count"){
#     type = "Rarefied Count"
#   }
#   else if(data.type == "arcsin"){
#     type = "Arcsine-Root"
#   }
#   
#   for(taxon.name in feature){
#     val <- numeric()
#     y_hat <- numeric()
#     taxon.val <- seq(min(X[,taxon.name]), max(X[,taxon.name]),len = 100)
#     
#     for(i in 1:length(taxon.val)){
#       newX <- X
#       newX[,taxon.name] <- rep(taxon.val[i], nrow(newX))
#       y_pred <- predict(fit, newX)
#       val <- c(val, taxon.val[i])
#       y_hat <- c(y_hat, mean(y_pred))
#     }
#     pred.result <- data.frame(val, y_hat) %>%
#       reshape2::melt(id.vars = "val", value.name = "Prediction")
#     pred.result$title <- rep(taxon.name, nrow(pred.result))
#     result <- rbind(result, pred.result)
#   }
#   result$title <- factor(result$title, levels = feature)
#   
#   p <- ggplot(result, aes(val, Prediction)) + 
#     geom_line(size = 0.8) +
#     theme_light() +
#     xlab(type) + 
#     ylab("Predicted Value") +
#     theme(
#       axis.text.x = element_text(size = 7.5),
#       axis.text.y = element_text(size = 7.5),
#       strip.text = element_text(size=12),
#       panel.spacing.x = unit(1, "lines"),
#       legend.title = element_blank()
#     ) +
#     facet_wrap(~ title, scales = "free_x", nrow = 5, dir = "v")
#   p
# }

rf.pdp.mult <- function(rf.list, n, name, data.type, level.names, label){
  X <- rf.list[[name]]$data$x
  fit <- rf.list[[name]]$fit
  rf.importance <- rf.imp.df(rf.list, name, type = 2) %>%
    mutate(taxon = rownames(rf.imp.df(rf.list, name, type = 2))) %>%
    arrange(-MeanDecreaseGini) %>%
    rownames
  n <- min(length(rf.importance), n)
  feature <- rf.importance[1:n]
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
    taxon.val <- seq(min(X[,taxon.name]), max(X[,taxon.name]),len = 100)
    
    for(i in 1:length(taxon.val)){
      newX <- X
      newX[,taxon.name] <- rep(taxon.val[i], nrow(newX))
      y_pred_prob <- predict(fit, newX, type = "prob")
      y_pred_prob_mean <- apply(y_pred_prob, 2, mean)
      val <- c(val, taxon.val[i])
      for(j in 1:label.len){
        prob[[j]] <- c(prob[[j]], y_pred_prob_mean[j])
      }
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
  return(p)
}
