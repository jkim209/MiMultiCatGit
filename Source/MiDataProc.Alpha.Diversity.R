library(tibble)
library(dplyr)
library(forestplot)
library(VGAM)

source("Source/MiDataProc.Data.Input.R")
source("Source/MiDataProc.Data.Upload.R")


data_mani <- function(alpha_div, sam_dat, prim_id, covariate = NULL){
  
  # No covariates
  
  if(is.null(covariate)){
    dat = cbind(sam_dat[,prim_id], alpha_div)
    colnames(dat)[1] <- "prime_var"
  }
  
  else{
    dat = cbind(sam_dat[,prim_id], alpha_div, sam_dat[,covariate])
    colnames(dat)[1] <- "prime_var"
  }
  
  return(dat)
}

decimal_adjust <- function(x){
  x <- as.numeric(x)
  round_num <- round(x, 3)
  decimal <- formatC(round_num, 3, format = "f")
  return(decimal)
}

p.value.0.1_char <- function(x){
  ind.0 <- which(x == "0.000")
  x[ind.0] <- "<.001"
  ind.1 <- which(x == "1.000")
  x[ind.1] <- ">.999"
  return(x)
}

alpha.bin.cat.func <- function(sam.dat, sel.bin.var) {
  bin.cat <- levels(as.factor(unlist(sam.dat[,sel.bin.var])))
  
  return(bin.cat)
}

bin.cat.recode.func.mult <- function(sam.dat, sel.bin.var, ori.cat, rename){
  ind <- list()
  for (i in 1:length(table(sam.dat[,sel.bin.var]))){
    ind <- which(sam.dat[,sel.bin.var] == ori.cat[i])
    sam.dat[ind, sel.bin.var] <- rename[i]
  }
  
  return(sam.dat)
}

no.missing.var <- function(sam.dat, sel.bin.var, covariate = NULL){
  
  # No covariate
  if(is.null(covariate)){
    sam.dat <- sam.dat[!is.na(sam.dat[[sel.bin.var]]), ]
  }
  
  else{
    ind <- sort(unique(which(!is.na(sam.dat[,c(sel.bin.var, covariate)]), arr.ind = TRUE)[,1]))
    sam.dat <- sam.dat[ind, ]
  }
  
  return(sam.dat)
}

no.missing.alpha <- function(alpha.div, sam.dat, sel.bin.var, covariate = NULL){
  
  # No covariate
  if(is.null(covariate)){
    alpha.div <- alpha.div[!is.na(sam.dat[[sel.bin.var]]), ]
  }
  
  else{
    ind <- sort(unique(which(!is.na(sam.dat[,c(sel.bin.var, covariate)]), arr.ind = TRUE)[,1]))
    alpha.div <- alpha.div[ind, ]
  }
  
  return(alpha.div)
}

alpha.mult.pair.cat.ref.func <- function(sel.bin.var, rename.var, sam.dat, alpha.div){
  bin.var <- unlist(sam.dat[,sel.bin.var], use.names = FALSE)
  
  bin.var.vec <- c()
  alpha_div <- c() 
  ind <- list()
  
  for (i in 1:length(rename.var)){
    ind[[i]] <- which(bin.var == rename.var[i])
    
    bin.var.vec <- c(bin.var.vec, bin.var[ind[[i]]])
    alpha_div <- rbind(alpha.div, alpha.div[ind[[i]],])
  }
  
  bin.var.vec <- factor(bin.var.vec)
  return(list(bin.var = bin.var.vec, alpha.div = alpha_div))
}

alpha.f.pair.mult.overall <- function(alpha_mani){
  
  summary_vec <- c()
  
  for (i in 1:(length(alpha_mani)-1)){
    summary_aov <- unlist(summary(aov(alpha_mani[, i+1] ~ prime_var, data = alpha_mani)))
    summary_ind <- c(decimal_adjust(summary_aov[["F value1"]]), decimal_adjust(summary_aov[["Pr(>F)1"]])) 
    summary_vec <- rbind(summary_vec, summary_ind)
  }
  
  rownames(summary_vec) = colnames(alpha_mani)[2:10]
  
  colnames(summary_vec) = c("F", "P.value")
  return(summary_vec)
}


alpha.f.pair.mult.tukey <- function(alpha_mani){
  div_list <- list() 
  
  for (i in 1:(length(alpha_mani)-1)){
    print(i)
    rmanova_result <- aov(alpha_mani[, i+1] ~ prime_var, data = alpha_mani)
    dat <- as.data.frame(TukeyHSD(rmanova_result)$prime_var)
    
    tukey_result <- data.frame(lapply(dat, decimal_adjust), row.names = rownames(dat)) 
    ind_result <- c() 
    
    for (j in 1:nrow(tukey_result)){
      ref_conf_list <- strsplit(rownames(TukeyHSD(rmanova_result)$prime_var), split = "-")[[j]]
      ref <- ref_conf_list[2]
      conf <- ref_conf_list[1]
      ref_mean <- decimal_adjust(mean(alpha_mani[,i+1][alpha_mani$prime_var == ref]))
      conf_mean <- decimal_adjust(mean(alpha_mani[,i+1][alpha_mani$prime_var == conf]))
      tukey_result_ind <- c(ref, conf, ref_mean, conf_mean, tukey_result[j,1:3], p.value.0.1_char(decimal_adjust(tukey_result[j, 4])))
      ind_result <- rbind(ind_result, tukey_result_ind)
    }
    
    rownames(ind_result) <- 1:nrow(ind_result)
    colnames(ind_result) <- c("Ref", "Com", "Mean (Ref) ", "Mean (Com)", "Diff", "Lower", "Upper", "Adj. P.value")
    div_list[[i]] <- as.data.frame(ind_result)
  }
  
  names(div_list) <- names(alpha_mani)[2:10]
  return(div_list)
}

alpha.kruskal.pair.mult.overall <- function(alpha_mani){
  
  summary_vec <- c()
  
  for (i in 1:(length(alpha_mani)-1)){
    summary_kruskal <<- kruskal.test(alpha_mani[, i+1] ~ prime_var, data = alpha_mani)
    summary_ind <- c(decimal_adjust(summary_kruskal[["statistic"]]), decimal_adjust(summary_kruskal[["p.value"]])) 
    summary_vec <- rbind(summary_vec, summary_ind)
  }
  
  rownames(summary_vec) = colnames(alpha_mani)[2:10]
  
  colnames(summary_vec) = c("Chi-sq", "P.value")
  return(summary_vec)
}

alpha.f.pair.mult.dunn <- function(alpha_mani){
  div_list <- list() 
  
  for (i in 1:(length(alpha_mani)-1)){
    name <- colnames(alpha_mani)[i+1]
    
    dunn_result <- dunnTest(alpha_mani[[name]] ~ alpha_mani[["prime_var"]])
    result <- dunn_result$res
    ind_result <- c() 
    
    for (j in 1:nrow(result)){
      ref_conf_list <- result[j,]
      ref_con <- strsplit(ref_conf_list[[1]], " - ")
      ref <- ref_con[[1]][1]
      conf <- ref_con[[1]][2]
      ref_mean <- decimal_adjust(mean(alpha_mani[,i+1][alpha_mani$prime_var == ref]))
      conf_mean <- decimal_adjust(mean(alpha_mani[,i+1][alpha_mani$prime_var == conf]))
      dunn_result_ind <- c(ref, conf, ref_mean, conf_mean, decimal_adjust(ref_conf_list[2]), p.value.0.1_char(decimal_adjust(ref_conf_list [4])))
      ind_result <- rbind(ind_result, dunn_result_ind)
    }
    
    rownames(ind_result) <- 1:nrow(ind_result)
    colnames(ind_result) <- c("Ref", "Com", "Mean (Ref) ", "Mean (Com)", "Z", "Adj. P.value")
    div_list[[i]] <- as.data.frame(ind_result)
  }
  
  names(div_list) <- names(alpha_mani)[2:10]
  return(div_list)
}
# 
# alpha_div <- alpha.div
# alpha_mani <- data.out
# test_result <- p.val.f.test

#0727 here 
alpha.bin.f.mult <- function(alpha_div, alpha_mani, test_result){
  ind <- list()
  time.p.cat <- names(table(alpha_mani[, "prime_var"]))
  for (i in 1:length(time.p.cat)){
    ind[[i]] <- which(alpha_mani[, "prime_var"] == time.p.cat[i])
  }
  
  each_alpha <- list() 
  
  for (j in 1:length(alpha_div)){
    alpha_dat <- list() 
    
    for (i in 1:length(time.p.cat)){
      alpha_dat[[i]] <- alpha_mani[ind[[i]], 1+j]
    }
    
    names(alpha_dat) <- time.p.cat
    each_alpha[[j]] <- alpha_dat
  }
  
  ind.p.sig <- which(as.numeric(test_result) < 0.05)
  
  length <- length(time.p.cat)
  if (length <= 3){
    par(mfrow = c(3, 3), mar = c(4.3, 4.1, 1, 1.9)) 
    
    for (i in 1:length(alpha_div)){
      if (is.element(i, ind.p.sig)){
        xlab.v <- paste("*p:", p.value.0.1_char(decimal_adjust(as.numeric(test_result[i]))), sep = "")
      } else {
        xlab.v <- paste("p:", p.value.0.1_char(decimal_adjust(as.numeric(test_result[i]))), sep = "")
      }
      
      boxplot(each_alpha[[i]], ylab = names(alpha_div)[i], xlab = xlab.v, names = time.p.cat,  notch = TRUE, horizontal = FALSE, col = "#DFC5F2", boxwex=0.7) #col = rainbow(length),
    }
    
  }else if (4 <= length & length <= 6){
    par(mfrow = c(5, 2), mar = c(4.1, 4.1, 1, 1.9)) 
    for (i in 1:length(alpha_div)){
      if (is.element(i, ind.p.sig)){
        xlab.v <- paste("*p:", p.value.0.1_char(decimal_adjust(as.numeric(test_result[i]))), sep = "")
      } else {
        xlab.v <- paste("p:", p.value.0.1_char(decimal_adjust(test_result[i])), sep = "")
      }
      
      boxplot(each_alpha[[i]], ylab = names(alpha_div)[i], xlab = xlab.v, names = time.p.cat, notch = TRUE, col = "#DFC5F2", horizontal = FALSE, boxwex=0.7) #col = rainbow(length),
    }
  }else {
    par(mfrow = c(9, 1), mar = c(4.1, 4.1, 1, 1.9))
    for (i in 1:length(alpha_div)){
      if (is.element(i, ind.p.sig)){
        xlab.v <- paste("*p:", p.value.0.1_char(decimal_adjust(as.numeric(test_result[i]))), sep = "")
      } else {
        xlab.v <- paste("p:", p.value.0.1_char(decimal_adjust(as.numeric(test_result[i]))), sep = "")
      }
      boxplot(each_alpha[[i]], ylab = names(alpha_div)[i], xlab = xlab.v,  names = time.p.cat,  notch = TRUE, col = "#DFC5F2", horizontal = FALSE, boxwex=0.7) #col = rainbow(length), 
    }
  }
}

# Multinomial logit model (nominal) Forest Plot preparation
alpha.nominal.result.prep <- function(logit.model, num.cat){
  
  # Number of categories - 1
  rm.int <- 1:(num.cat - 1)
  
  # Coefficients - remove the rows with intercpets and the columns with z value
  ctable <- coef(summary(logit.model))[-rm.int, c(1,2,4)]
  
  # Confidence interval
  conf.int <- confint(logit.model)[-rm.int,]
  
  # Merge the table
  ctable <- data.frame(cbind(ctable, conf.int))
  
  # Exclude the covariate
  ctable <- head(ctable, 9*length(rm.int))
  
  # Change the rownames to the fitted format
  colnames(ctable) <- c("Estimate", "SE", "P-value", "LLCI", "ULCI")
  
  rownames.end.with <- substr(rownames(ctable), nchar(rownames(ctable)), nchar(rownames(ctable)))
  
  ctable.list <- list()
  
  for(n in unique(rownames.end.with)){
    ind <- which(rownames.end.with == n)
    ctable.list[[n]] <- ctable[ind,]
    rownames(ctable.list[[n]]) <- c("Observed", "Shannon", "Simpson", "InvSimpson", "Fisher", "Chao1", "ACE", "ICE", "PD")
  }
  
  return(ctable.list)
}

# Nominal Forest Plot
alpha.nominal.forest.plot <- function(ctable){
  
  Mean = ctable$Estimate
  SE = ctable$SE
  P_value = ctable$`P-value`
  Lower = ctable$LLCI
  Upper = ctable$ULCI
  
  # Adjust the format of the values
  Estimate <- ifelse(abs(round(Mean, 3)) == 0, "<.001", sprintf("%.3f", round(Mean, 3)))
  SE <- ifelse(abs(round(SE, 3)) == 0, "<.001", sprintf("%.3f", round(SE, 3)))
  P_value <- ifelse(abs(round(P_value, 3)) == 0, "<.001", sprintf("%.3f", round(P_value, 3)))
  
  # Change from data frame to tibble format
  output_base_data <- tibble(mean = Mean,
                             lower = Lower,
                             upper = Upper,
                             alpha_div = rownames(ctable),
                             estimate = as.character(Estimate),
                             se = as.character(SE),
                             p_value = as.character(P_value))
  
  # Add header
  output_header <- tibble(alpha_div = "Alpha Diversity",
                          estimate = "Estimate",
                          se = "SE",
                          p_value = "P-value",
                          summary = TRUE)
  
  # Apply the header to the tibble
  output_df <- bind_rows(output_header,
                         output_base_data)
  
  # Forest Plot
  output_df %>% 
    forestplot(labeltext = c(alpha_div, estimate, se, p_value),
               is.summary = summary,
               xlab = "95% Confidence Interval",
               txt_gp = fpTxtGp(title = gpar(cex=1.35),
                                ticks = gpar(cex=0.75),
                                xlab = gpar(cex = 1.1)),
               hrzl_lines = list("2" = gpar(lty = 2), 
                                 "11" = gpar(lwd = 1, columns = 1:4, col = "#000044")),
               # title = "Multinomial Logistic Regression",
               col = fpColors(box = "#FFBE3C",
                              line = "golden rod"),
               boxsize = 0.15,
               alpha = 0.75)
}

# Proportional odds model Forest Plot preparation
alpha.ordinal.result.prep <- function(logit.model, num.cat){
  # Number of categories - 1
  rm.int <- 1:(num.cat - 1)
  
  # Coefficients - remove the rows with intercpets and the columns with z value
  ctable <- coef(summary(logit.model))[-rm.int, c(1,2,4)]
  
  # Confidence interval
  conf.int <- confint(logit.model)[-rm.int,]
  
  # Merge the table
  ctable <- data.frame(cbind(ctable, conf.int))
  
  # Exclude the covariate
  ctable <- head(ctable, 9)
  
  # Change the rownames to the fitted format
  colnames(ctable) <- c("Estimate", "SE", "P-value", "LLCI", "ULCI")
  
  return(ctable)
}

# Ordinal Forest Plot
alpha.ordinal.forest.plot <- function(ctable){
  
  Mean = ctable$Estimate
  SE = ctable$SE
  P_value = ctable$`P-value`
  Lower = ctable$LLCI
  Upper = ctable$ULCI
  
  # Adjust the format of the values
  Estimate <- ifelse(abs(round(Mean, 3)) == 0, "<.001", sprintf("%.3f", round(Mean, 3)))
  SE <- ifelse(abs(round(SE, 3)) == 0, "<.001", sprintf("%.3f", round(SE, 3)))
  P_value <- ifelse(abs(round(P_value, 3)) == 0, "<.001", sprintf("%.3f", round(P_value, 3)))
  
  # Change from data frame to tibble format
  output_base_data <- tibble(mean = Mean,
                             lower = Lower,
                             upper = Upper,
                             alpha_div = rownames(ctable),
                             estimate = as.character(Estimate),
                             se = as.character(SE),
                             p_value = as.character(P_value))
  
  # Add header
  output_header <- tibble(alpha_div = "Alpha Diversity",
                          estimate = "Estimate",
                          se = "SE",
                          p_value = "P-value",
                          summary = TRUE)
  
  # Apply the header to the tibble
  output_df <- bind_rows(output_header,
                         output_base_data)
  
  # Forest Plot
  output_df %>% 
    forestplot(labeltext = c(alpha_div, estimate, se, p_value),
               is.summary = summary,
               xlab = "95% Confidence Interval",
               txt_gp = fpTxtGp(title = gpar(cex=1.35),
                                ticks = gpar(cex=0.75),
                                xlab = gpar(cex = 1.1)),
               hrzl_lines = list("2" = gpar(lty = 2), 
                                 "11" = gpar(lwd = 1, columns = 1:4, col = "#000044")),
               # title = "Proportional Odds Model",
               col = fpColors(box = "#8f9940",
                              line = "darkgreen"),
               boxsize = 0.15,
               alpha = 0.75)
}

alpha.nominal.global.test <- function(data.out, covariate = NULL){
  
  if(is.null(covariate)){
    alpha.list <- c("Observed", "Shannon", "Simpson", "InvSimpson", "Fisher", "Chao1", "ACE", "ICE", "PD")
    nominal.null <- vglm(prime_var ~ 1, data = data.out, family = multinomial)
    
    nominal.global.summary.df <- data.frame()
    
    for(alpha in alpha.list){
      nominal.fit <- vglm(as.formula(paste0("prime_var", "~", alpha)), data = data.out, family = multinomial)
      # nominal.lrt.summary <- lrtest(nominal.fit, nominal.null)
      nominal.lrt.summary <- lrtest(nominal.null, nominal.fit)
      nominal.global.summary.df <- rbind(nominal.global.summary.df, data.frame(alpha = alpha, subset(nominal.lrt.summary@Body, select = -1)[2,]))
    }
    colnames(nominal.global.summary.df) <- c("Alpha Diversity", "Log Likelihood", "DF", "Chisq", "P-value")
    nominal.global.summary.df$`Log Likelihood` <- round(nominal.global.summary.df$`Log Likelihood`, 3)
    nominal.global.summary.df$Chisq <- round(nominal.global.summary.df$Chisq, 3)
    nominal.global.summary.df$`P-value` <- ifelse(round(nominal.global.summary.df$`P-value`, 3) == 0, "<.001", sprintf("%.3f", round(nominal.global.summary.df$`P-value`, 3)))
  }
  else{
    alpha.list <- c("Observed", "Shannon", "Simpson", "InvSimpson", "Fisher", "Chao1", "ACE", "ICE", "PD")
    nominal.null <- vglm(as.formula(paste("prime_var", "~", paste(c(covariate), collapse = "+"))), data = data.out, family = multinomial)
    
    nominal.global.summary.df <- data.frame()
    
    for(alpha in alpha.list){
      nominal.fit <- vglm(as.formula(paste("prime_var", "~", paste(c(covariate, alpha), collapse = "+"))), data = data.out, family = multinomial)
      # nominal.lrt.summary <- lrtest(nominal.fit, nominal.null)
      nominal.lrt.summary <- lrtest(nominal.null, nominal.fit)
      nominal.global.summary.df <- rbind(nominal.global.summary.df, data.frame(alpha = alpha, subset(nominal.lrt.summary@Body, select = -1)[2,]))
    }
    colnames(nominal.global.summary.df) <- c("Alpha Diversity", "Log Likelihood", "DF", "Chisq", "P-value")
    nominal.global.summary.df$`Log Likelihood` <- round(nominal.global.summary.df$`Log Likelihood`, 3)
    nominal.global.summary.df$Chisq <- round(nominal.global.summary.df$Chisq, 3)
    nominal.global.summary.df$`P-value` <- ifelse(round(nominal.global.summary.df$`P-value`, 3) == 0, "<.001", sprintf("%.3f", round(nominal.global.summary.df$`P-value`, 3)))
  }
  rownames(nominal.global.summary.df) <- NULL
  nominal.global.summary.df$DF <- abs(nominal.global.summary.df$DF)
  return(nominal.global.summary.df)
}
