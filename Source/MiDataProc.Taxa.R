library(VGAM)

taxa.united.func <- function(sel.bin.var, sam.dat, taxa) {
  bin.var <- unlist(sam.dat[,sel.bin.var])
  taxa.out <- taxa
  return(list(bin.var = bin.var, taxa = taxa.out))
}

taxa.data.mani <- function(sam_dat, taxa, prim_id, type, covariate = NULL){
  taxa <- taxa[[type]]
  taxa.list <- list() 
  
  if(is.null(covariate)){
    for (i in 1:length(taxa)){
      taxa.ind <- taxa[[i]]
      new.taxa.ind <- cbind(prim.var = sam_dat[[prim_id]], taxa.ind)
      new.taxa.ind <- new.taxa.ind[!is.na(sam_dat[[prim_id]]), ]
      
      taxa.list[[i]] <- new.taxa.ind
    }
    names(taxa.list) <- names(taxa)
  }
  
  else{
    for(i in 1:length(taxa)){
      taxa.ind <- taxa[[i]]
      new.taxa.ind <- cbind(prim.var = sam_dat[[prim_id]], sam_dat[,covariate], taxa.ind)
      # new.taxa.ind <- new.taxa.ind[!is.na(sam_dat[[prim_id]]), ]
      new.taxa.ind <- na.omit(new.taxa.ind)
      
      taxa.list[[i]] <- new.taxa.ind
    }
    names(taxa.list) <- names(taxa)
  }

  return(taxa.list)
}

taxa.f.pair.mult.overall <- function(taxa_mani, include){
  
  output_list <- list() 
  summary_list <- list() 
  pval_list <- list() 
  
  for (i in 1:(5+include)){
    taxa.ind <- taxa_mani[[i]]
    pval_vec <- c() 
    
    out <- matrix(NA, nrow = length(taxa.ind)-1, ncol = 2)
    for (j in 1:(length(taxa.ind)-1)){
      summary_aov <- unlist(summary(aov(taxa.ind[,j+1] ~ prim.var, data = taxa.ind)))
      p_val <- as.vector(summary_aov[["Pr(>F)1"]])
      pval_vec <- c(pval_vec, p_val)
      summary_ind <- c(as.vector(summary_aov[["F value1"]]), as.vector(summary_aov[["Pr(>F)1"]]))
      out[j,] <- summary_ind
    }
    
    pval_list[[i]] <- decimal_adjust(p.adjust(pval_vec, "BH"))
    rownames(out) <- colnames(taxa.ind[2:length(taxa.ind)])
    colnames(out) <- c("F", "P.value")
    summary_list[[i]] <- out
  }
  
  names(summary_list) <- names(taxa_mani)[1:(5+include)]
  names(pval_list) <- names(taxa_mani)[1:(5+include)]
  output_list$global_test <- summary_list 
  output_list$pval <- pval_list 
  
  return(output_list) 
}

taxa.kruskal.mult.overall <- function(taxa_mani, include){
  
  output_list <- list() 
  summary_list <- list() 
  pval_list <- list() 
  
  for (i in 1:(5+include)){
    taxa.ind <- taxa_mani[[i]]
    pval_vec <- c() 
    
    out <- matrix(NA, nrow = length(taxa.ind)-1, ncol = 2)
    for (j in 1:(length(taxa.ind)-1)){
      
      summary_kruskal <<- kruskal.test(taxa.ind[, j+1] ~ prim.var, data = taxa.ind)
      
      summary_ind <- c(summary_kruskal[["statistic"]], summary_kruskal[["p.value"]])
      p_val <- summary_kruskal[["p.value"]]
      pval_vec <- c(pval_vec, p_val)
      out[j,] <- summary_ind 
    }
    
    pval_list[[i]] <- decimal_adjust(p.adjust(pval_vec, "BH"))
    rownames(out) <- colnames(taxa.ind[2:length(taxa.ind)])
    colnames(out) <- c("F", "P.value")
    summary_list[[i]] <- out
  }
  
  names(summary_list) <- names(taxa_mani)[1:(5+include)]
  names(pval_list) <- names(taxa_mani)[1:(5+include)]
  output_list$global_test <- summary_list 
  output_list$pval <- pval_list 
  
  return(output_list) 
}


q_val_combined_table.anova <- function(result_table, q_val_list){
  new_list_with_q <- list() 
  for (i in 1:length(result_table)){
    new_list_with_q[[i]] <- cbind(result_table[[i]], as.vector(q_val_list[[i]]))
    colnames(new_list_with_q[[i]]) <- c("F", "P.value", "Adj.P.value")
  }
  names(new_list_with_q) <- names(result_table)
  return(new_list_with_q)
}

q_val_combined_table.kruskal <- function(result_table, q_val_list){
  new_list_with_q <- list() 
  for (i in 1:length(result_table)){
    new_list_with_q[[i]] <- cbind(result_table[[i]], as.vector(q_val_list[[i]]))
    colnames(new_list_with_q[[i]]) <- c("Chi-sq", "P.value", "Adj.P.value")
  }
  names(new_list_with_q) <- names(result_table)
  return(new_list_with_q)
}

taxa.f.pair.mult.tukey <- function(taxa_mani, include){
  summary_list <- list() 
  
  for (i in 1:(5+include)){
    taxa.ind <- taxa_mani[[i]]
    each_list <- list() 
    for (j in 1:(length(taxa.ind)-1)){
      rmanova_result <- aov(taxa.ind[,j+1] ~ prim.var , data = taxa.ind)
      tukey_result <- TukeyHSD(rmanova_result)$prim.var
      
      ind_result <- c() 
      for (k in 1:nrow(tukey_result)){
        ref_conf_list <- strsplit(rownames(TukeyHSD(rmanova_result)$prim.var), split = "-")[[k]]
        ref <- ref_conf_list[2]
        conf <- ref_conf_list[1]
        
        ref_mean <- as.vector(mean(taxa.ind[,j+1][taxa.ind$prim.var == ref]))
        conf_mean <- as.vector(mean(taxa.ind[,j+1][taxa.ind$prim.var == conf]))
        
        tukey_result_ind <- c(ref, conf, decimal_adjust(c(ref_mean, conf_mean, tukey_result[k,])))
        ind_result <- rbind(ind_result, tukey_result_ind)
      }
      rownames(ind_result) <- 1:nrow(ind_result)
      colnames(ind_result) <- c("Ref", "Com", "Mean (Ref)", "Mean (Com)", "Diff", "Lower", "Upper", "Adj.P.value")
      each_list[[j]] <- as.data.frame(ind_result)
    }
    names(each_list) <-  colnames(taxa.ind[2:length(taxa.ind)])
    summary_list[[i]] <- each_list 
  }
  names(summary_list) <- names(taxa_mani)[1:(5+include)]
  return(summary_list)
}

taxa.pair.mult.dunn <- function(taxa_mani, include){
  summary_list <- list() 
  
  for (i in 1:(5+include)){
    taxa.ind <- taxa_mani[[i]]
    each_list <- list() 
    for (j in 1:(length(taxa.ind)-1)){
      
      name <- colnames(taxa.ind)[j+1]
      dunn_result <- dunnTest(taxa.ind[[name]] ~ taxa.ind[["prim.var"]])
      result <- dunn_result$res
      
      ind_result <- c() 
      
      for (k in 1:nrow(result)){
        ref_conf_list <- result[k,]
        ref_con <- strsplit(ref_conf_list[[1]], " - ")
        ref <- ref_con[[1]][1]
        conf <- ref_con[[1]][2]
        
        ref_mean <- as.vector(mean(taxa.ind[,j+1][taxa.ind$prim.var == ref]))
        conf_mean <- as.vector(mean(taxa.ind[,j+1][taxa.ind$prim.var == conf]))
        
        dunn_result_ind <- c(ref, conf, decimal_adjust(ref_mean), decimal_adjust(conf_mean), decimal_adjust(ref_conf_list[2]), decimal_adjust(ref_conf_list [4]))
        
        ind_result <- rbind(ind_result, dunn_result_ind)
      }
      rownames(ind_result) <- 1:nrow(ind_result)
      colnames(ind_result) <-c("Ref", "Com", "Mean (Ref) ", "Mean (Com)", "Z", "Adj.P.value")
      each_list[[j]] <- as.data.frame(ind_result)
    }
    
    names(each_list) <-  colnames(taxa.ind[2:length(taxa.ind)])
    summary_list[[i]] <- each_list 
  }
  names(summary_list) <- names(taxa_mani)[1:(5+include)]
  return(summary_list)
}


make_p_val_list <- function(result_table){
  list_total <- list() 
  for (i in 1:length(result_table)){
    list_total_ind <- list()
    for (j in 1:nrow(result_table[[i]][[1]])){
      vec <- NULL
      for (k in 1:length(result_table[[i]])){
        vec <- c(vec, result_table[[i]][[k]]$Adj.P.value[j])
      }
      list_total_ind[[j]] <- vec 
    }
    list_total[[i]] <- list_total_ind
  }
  return(list_total)
}

make_q_val_list <- function(p_val_list){
  q_val_result <- list() 
  for (p in 1:length(p_val_list)){
    q_val_result_ind <- list()
    for (k in 1:length(p_val_list[[p]])){
      q_val_result_ind[[k]] <- p.adjust(p_val_list[[p]][[k]], "BH")
    }
    q_val_result[[p]] <- q_val_result_ind
  }
  return(q_val_result)
}

q_val_dat <- function(q_val_list, result_table, p.val = TRUE){
  q_val_result_dat <- list()
  
  for (i in 1:length(q_val_list)){
    mat <- matrix(unlist(q_val_list[[i]]), ncol = length(q_val_list[[i]]), nrow = length(q_val_list[[i]][[1]]))
    q_val_result_dat[[i]] <- mat 
  }
  
  if(p.val){
    for (k in 1:length(result_table)){
      for (p in 1:length(result_table[[k]])){
        
        result_table[[k]][[p]]$P.value <- p.value.0.1_char(result_table[[k]][[p]]$P.value)
        result_table[[k]][[p]]$Adj.P.value <- p.value.0.1_char(decimal_adjust(q_val_result_dat[[k]][p,]))  
      }
    }
  }else{
    for (k in 1:length(result_table)){
      for (p in 1:length(result_table[[k]])){
        
        result_table[[k]][[p]]$Adj.P.value <- p.value.0.1_char(decimal_adjust(q_val_result_dat[[k]][p,]))  
      }
    }
  }
  
  return(result_table)
}


data_mani_p_taxa <- function(taxa, sam_dat, prim_id, block_id){
  dat = cbind(sam_dat[,prim_id], taxa)
  dat = dat[order(dat[,prim_id]),]
  colnames(dat)[1] <- "prime_var"
  return(dat)
}

taxa.bin.mult<- function(taxa, sam_dat, prim_id, page, q_val_result){
  
  taxon <- taxa[[page]]
  data.mani.taxa <- taxon
  
  time.p.cat <- names(table(data.mani.taxa[,1]))
  
  q_val_result[[page]][!complete.cases(as.numeric(q_val_result[[page]]))] <- 1 
  sig.num <- length(q_val_result[[page]][q_val_result[[page]] < 0.05])
  
  nrow <- ceiling(sig.num/4)
  
  if (nrow > 0){
    par(mfrow = c(nrow, 4))
    
    ind <- list()
    for (i in 1:length(time.p.cat)){
      ind[[i]] <- which(data.mani.taxa[,1] == time.p.cat[i])
    }
    
    each_in_taxon <- list() 
    meaningful = c() 
    taxa_dat_list <- list() 
    id <- 0 
    
    for (j in 1:length(taxon))
    {
      if (j %in% which(q_val_result[[page]] < 0.05)){ 
        meaningful = c(meaningful, j)
        id <- id+1
        taxa_dat <- c() 
        
        for (i in 1:length(time.p.cat)){
          taxa_dat <- cbind(taxa_dat, data.mani.taxa[, 1+j][ind[[i]]])
        }
        
        colnames(taxa_dat) <- time.p.cat
        taxa_dat_list[[id]] <- taxa_dat
        
        for (k in 1:length(meaningful)){
          each_in_taxon[[k]] <- taxa_dat
        }
      }
    }
    
    xlab.v <- c() 
    for (i in 1:length(q_val_result[[page]][q_val_result[[page]]<0.05])){
      xlab.v[i] <- paste0("*q:", p.value.0.1_char(decimal_adjust(q_val_result[[page]][q_val_result[[page]]<0.05][i])), sep = "")
    }
    
    ylab.v <- c() 
    for (i in 1:length(q_val_result[[page]][q_val_result[[page]]<0.05])){
      ylab.v[i] <- rev(strsplit(names(taxon[,2:ncol(taxon)])[q_val_result[[page]] < 0.05][[i]], ";")[[1]])[1]
    }
    
    for (i in 1:length(q_val_result[[page]][q_val_result[[page]]<0.05])){
      boxplot(taxa_dat_list[[i]], main = i, xlab = xlab.v[i], ylab = ylab.v[i], names = time.p.cat, col = "#DFC5F2",  notch = TRUE , boxwex=0.9, horizontal = FALSE)
    }
    
  }else{
    ranks <- c("Phylum", "Class", "Order", "Family", "Genus", "Species")
    plot.new()
    text(x = 0.5, y = 0.5, paste("No significant taxa are found in ", ranks[page], sep = ""), cex = 1.2, col = "black")
  }
}

sig.taxa <- function(p_val){
  
  sig_list <- list() 
  for (i in 1:length(p_val)){
    sig_list[[i]] <- as.numeric(p_val[[i]]) < 0.05 
    sig_list[[i]][!complete.cases(as.numeric(p_val[[i]]))] <- FALSE
  }
  return(sig_list)
}

rank_sig <- function(pairwise_result, sig.taxa.list){
  sig.global.list <- list() 
  for (i in 1:length(pairwise_result)){
    sig.global.list[[i]] <- pairwise_result[[i]][sig.taxa.list[[i]]]
  }
  names(sig.global.list) <- names(pairwise_result)
  return(sig.global.list)
}

list_to_dat <- function(rank_sig_list){
  dat_list <- list() 
  
  for (i in 1:length(rank_sig_list)){
    names(rank_sig_list[[i]]) <- NULL
    dat_list[[i]] <- do.call(rbind.data.frame, rank_sig_list[[i]])
  }
  names(dat_list) <- names(rank_sig_list)
  return(dat_list)
}

taxa.ordinal.result.prep <- function(taxa.mani, name, covariate = NULL){
  level.ctable <- data.frame()
  taxon.df <- taxa.mani[[name]]
  
  if(is.null(covariate)){
    taxon.names <- colnames(taxon.df)[-which(names(taxon.df) %in% c("prim.var"))]
  }
  else{
    taxon.names <- colnames(taxon.df)[-which(names(taxon.df) %in% c("prim.var", covariate))]
  }
  
  for(taxon in taxon.names){
    if(is.null(covariate)){
      df <- subset(taxon.df, select = c("prim.var", taxon))
    }
    else{
      df <- subset(taxon.df, select = c("prim.var", covariate, taxon))
    }
    set.seed(578)
    fit <- try(vglm(prim.var ~ ., data = df, family = cumulative(parallel = TRUE)), silent = TRUE)
    
    ctable <- coef(summary(fit))[, c(1,2,4), drop = FALSE]
    conf.int <- confint(fit)
    ctable <- data.frame(cbind(ctable, conf.int))
    colnames(ctable) <- c("Estimate", "SE", "P-value", "LLCI", "ULCI")
    
    # ctable <- ctable[ncat,]
    ctable <- tail(ctable, 1)
    ctable$Rank <- str_to_title(name)
    ctable$Taxon <- taxon
    
    level.ctable <- rbind(level.ctable, ctable)
  }
  level.ctable$'Q-value' <- p.adjust(level.ctable$'P-value', "BH")
  # result.ctable <- rbind(result.ctable, level.ctable)
  return(level.ctable)
}

taxa.ordinal.result.merge <- function(taxa.mani, level.names, covariate = NULL){
  
  result.ctable <- data.frame()
  
  for(name in level.names){
    level.ctable <- taxa.ordinal.result.prep(taxa.mani, name, covariate)
    result.ctable <- rbind(result.ctable, level.ctable)
  }
  result.ctable <- result.ctable[,c(6, 7, 1, 2, 3, 8, 4, 5)]
  rownames(result.ctable) <- NULL
  
  return(result.ctable)
}

taxa.nominal.result.prep <- function(taxa.mani, name, ncat, covariate = NULL){
  level.ctable <- data.frame()
  taxon.df <- taxa.mani[[name]]
  
  level.ctable <- data.frame()
  taxon.df <- taxa.mani[[name]]
  
  if(is.null(covariate)){
    taxon.names <- colnames(taxon.df)[-which(names(taxon.df) %in% c("prim.var"))]
  }
  else{
    taxon.names <- colnames(taxon.df)[-which(names(taxon.df) %in% c("prim.var", covariate))]
  }
  
  for(taxon in taxon.names){
    if(is.null(covariate)){
      df <- subset(taxon.df, select = c("prim.var", taxon))
    }
    else{
      df <- subset(taxon.df, select = c("prim.var", covariate, taxon))
    }
    set.seed(578)
    fit <- try(vglm(prim.var ~ ., data = df, family = multinomial), silent = TRUE)
    
    ctable <- coef(summary(fit))[, c(1,2,4), drop = FALSE]
    conf.int <- confint(fit)
    ctable <- data.frame(cbind(ctable, conf.int))
    colnames(ctable) <- c("Estimate", "SE", "P-value", "LLCI", "ULCI")
    
    # ctable <- ctable[ncat-1,]
    ctable <- tail(ctable, ncat-1)
    ctable$Rank <- str_to_title(name)
    ctable$Taxon <- taxon
    
    level.ctable <- rbind(level.ctable, ctable)
  }
  
  rownames.end.with <- substr(rownames(level.ctable), nchar(rownames(level.ctable)), nchar(rownames(level.ctable)))
  
  ctable.list <- list()
  
  for(n in unique(rownames.end.with)){
    ind <- which(rownames.end.with == n)
    ctable.list[[n]] <- level.ctable[ind,]
    ctable.list[[n]]$'Q-value' <- p.adjust(ctable.list[[n]]$'P-value', "BH")
    rownames(ctable.list[[n]]) <- taxon.names
  }
  return(ctable.list)
}

taxa.nominal.result.merge <- function(taxa.mani, level.names, ncat, covariate = NULL){
  
  ctable.list.by.taxa <- list()
  result.ctable.list <- list()
  
  for(name in level.names){
    ctable.list.by.taxa[[name]] <- taxa.nominal.result.prep(taxa.mani, name, ncat, covariate)
  }
  for(n in as.character(1:(ncat-1))){
    
    df <- data.frame()
    
    for(name in level.names){
      df <- rbind(df, ctable.list.by.taxa[[name]][[n]])
    }
    
    result.ctable.list[[n]] <- df
    result.ctable.list[[n]] <- result.ctable.list[[n]][,c(6, 7, 1, 2, 3, 8, 4, 5)]
    rownames(result.ctable.list[[n]]) <- NULL
  }
  
  return(result.ctable.list)
}

# How many pages should be required

taxa.forest.plot.pages <- function(result.ctable, mult.test.cor = "TRUE", thres = 0.05) {
  
  # Indices of significant taxon by taxonomic levels
  
  sig.by.rank <- list()
  
  # Total indices of significant taxon
  if(mult.test.cor) total <- sum(result.ctable$`Q-value` < thres)
  else total <- sum(result.ctable$`P-value` < thres)
  
  if(total>80) {
    plot.per.page <- total
    num.pages <- 0
    while(plot.per.page > 80){
      num.pages <- num.pages +1
      plot.per.page <- total%/%num.pages
    }
  } else if(total<=80 & total>=0) {
    num.pages = 1
    plot.per.page <-total
  }
  
  return(num.pages)
}

# Adjust the forest plot pages

taxa.forest.plot.pages1 <- function(result.ctable, thres = 0.05, mult.test.cor = TRUE) {
  sig.by.rank <- list()
  
  # P-G: 5 / P-S: 6
  
  range <- length(unique(result.ctable$Rank))
  
  # Change the report type
  
  # if(report.type == "Est"){
  #   report.txt <- "Est."
  # } else if(report.type == "OR") {
  #   report.txt <- "OR"
  # } else if(report.type == "Effect") {
  #   report.txt <- "Eff."
  # }
  
  # Total indices of significant taxon
  
  if(mult.test.cor) total <- sum(result.ctable$`Q-value` < thres)
  else total <- sum(result.ctable$`P-value` < thres)
  
  # Determine the page number(s)
  
  if(total>80) {
    capacity <- c(40:80)
    plot.per.page <- total
    num.pages <- 0
    num.mod <- 0
    while(plot.per.page > 80){
      num.pages <- num.pages +1
      plot.per.page <- total%/%num.pages
      num.mod <- total - total * (total%/%num.pages) # 수정 mod(total, num.pages)
    }
  } else if(total<=80 & total>=1) {
    num.pages = 1
    plot.per.page <-total
    num.mod <- 0
  }
  
  result <- list()
  if(total >0) {
    text.tab.all <- 0
    actual.plot.per.page <- plot.per.page
    
    info <- result.ctable
    if(mult.test.cor){
      ind <- which(info$'Q-value' < thres)
      info <- info[ind,]
    }
    else{
      ind <- which(info$'P-value' < thres)
      info <- info[ind,]
    }
    info$ID <- 1:nrow(info)
    # info$Estimate <- format(round(info$Estimate,3), nsmall = 3)
    # info$'P-value' <- p.value.0.1_char(info$'P-value')
    # info$'Q-value' <- p.value.0.1_char(info$'Q-value')
    info <- info[,c(9, 1, 2, 3, 4, 5, 6, 7, 8)]
    rownames(info) <- NULL
    
    #format(p.value.0.1(info.data[, c("P.val", "Q.val")]))))
    
    for(p in 1:num.pages) {
      if(p <= num.mod){
        actual.plot.per.page <- plot.per.page+1
        initial <- 1+actual.plot.per.page*(p-1)
      } else {
        initial <- 1+actual.plot.per.page*(p-1)
        actual.plot.per.page <- plot.per.page
      }
      # names(info) <- names(text.tab.all)
      
      result[[p]] <- info[c(initial:(initial+actual.plot.per.page-1)),]
      
      # all.text.tab[[p]] <- rbind(as.matrix(text.tab.all), info[c(initial:(initial+actual.plot.per.page-1)),])
      # colnames(info.ci) <- NULL
      # all.ci.tab[[p]] <- rbind(as.matrix(ci.tab.all), as.matrix(info.ci[c(initial:(initial+actual.plot.per.page-1)),]))
    }
  } else {
    # all.text.tab <- NULL
    # all.ci.tab <- NULL
    result <- NULL
  }
  # return(list(all.text.tab = all.text.tab, all.ci.tab = all.ci.tab))
  return(result)
}

taxa.forest.plot <- function(result.ctable.list, page) { # page.taxa.q.out에 taxa.forest.plot.page1이 들어가야 함
  
  result <- result.ctable.list[[page]]
  
  if(is.null(result)){
    plot.new()
    text(x = 0.5, y = 0.5, "No significant taxa are found.",
         cex = 1.2, col = "black")
  }
  # if(nrow(result) ==  0){
  #   plot.new()
  #   text(x = 0.5, y = 0.5, "No significant taxa are found.", 
  #        cex = 1.2, col = "black")
  # }
  else{
    ID <- as.character(result$ID)
    Rank <- result$Rank
    Taxon <- substr(result$Taxon, 1, 55)
    Mean <- result$Estimate
    P_value <- result$'P-value'
    Q_value <- result$'Q-value'
    Lower = result$LLCI
    Upper = result$ULCI
    
    # Adjust the format of the values
    
    Estimate <- ifelse(abs(round(Mean, 3)) == 0, "<.001", sprintf("%.3f", round(Mean, 3)))
    P_value <- ifelse(abs(round(P_value, 3)) == 0, "<.001", sprintf("%.3f", round(P_value, 3)))
    Q_value <- ifelse(abs(round(Q_value, 3)) == 0, "<.001", sprintf("%.3f", round(Q_value, 3)))
    
    # Change from data frame to tibble format
    # First three for the CI
    # Rest of the columns for the header
    
    output_base_data <- tibble(mean = Mean,
                               lower = Lower,
                               upper = Upper,
                               id = as.character(ID),
                               rank = Rank,
                               taxon_name = Taxon,
                               estimate = as.character(Estimate),
                               p_value = as.character(P_value),
                               q_value = as.character(Q_value))
    
    # Add header
    output_header <- tibble(id = "ID",
                            rank = "Rank",
                            taxon_name = "Taxon",
                            estimate = "Estimate",
                            p_value = "P-value",
                            q_value = "Q-value",
                            summary = TRUE)
    
    # Apply the header to the tibble
    output_df <- bind_rows(output_header,
                           output_base_data)
    
    
    output_df %>% 
      forestplot(labeltext = c(id, rank, taxon_name, estimate, p_value, q_value),
                 is.summary = summary,
                 xlab = "95% Confidence Interval",
                 hrzl_lines = TRUE,
                 new_page = TRUE, 
                 boxsize = 0.18, 
                 grid = 0, 
                 colgap = unit(1, "cm"), 
                 graphwidth = unit(8, "cm"), 
                 lineheight = "lines", 
                 line.margin = unit(0.12, "cm"),
                 col=fpColors(box=rgb(1,0,0,0.5), line="black", summary="red3"), 
                 mar = unit(c(0.5,0,0.5,0), "cm"), #mar = unit(c(blank.space,0,0,0), "cm"),
                 txt_gp=fpTxtGp(label=list(gpar(fontfamily="", cex=0.75), gpar(fontfamily="", cex=0.75)),
                                ticks=gpar(fontfamily="", cex=0.75),
                                xlab=gpar(fontfamily="", cex=0.75)))
  }
}

taxa.nominal.global <- function(taxa.mani, level.names, covariate = NULL){
  taxa.nominal.global.summary.list <- list()
  
  if(is.null(covariate)){
    for(name in level.names){
      data <- taxa.mani[[name]]
      taxon.names <- colnames(subset(data, select = -1))
      
      nominal.global.summary.df <- data.frame()
      
      nominal.null <- vglm(prim.var ~ 1, data = data, family = multinomial)
      
      for(taxon in taxon.names){
        nominal.fit <- vglm(prim.var ~., data = data[,c("prim.var", taxon)], family = multinomial)
        nominal.lrt.summary <- lrtest(nominal.null, nominal.fit)
        nominal.global.summary.df <- rbind(nominal.global.summary.df, data.frame(rank = str_to_title(name), taxon = taxon, subset(nominal.lrt.summary@Body, select = -1)[2,]))
      }
      nominal.global.summary.df$ID <- 1:nrow(nominal.global.summary.df)
      nominal.global.summary.df$Qval <- p.adjust(nominal.global.summary.df$Pr..Chisq., "BH")
      nominal.global.summary.df <- nominal.global.summary.df[,c(7, 1:6, 8)]
      colnames(nominal.global.summary.df) <- c("ID", "Rank", "Taxon", "Log Likelihood", "DF", "Chisq", "P-value", "Q-value")
      rownames(nominal.global.summary.df) <- NULL
      
      taxa.nominal.global.summary.list[[name]] <- nominal.global.summary.df
    }
  }
  else{
    for(name in level.names){
      data <- taxa.mani[[name]]
      taxon.names <- colnames(data[,!(names(data) %in% c("prim.var", covariate))])
      
      nominal.global.summary.df <- data.frame()
      
      nominal.null <- vglm(prim.var ~ ., data = data[,c("prim.var", covariate)], family = multinomial)
      
      for(taxon in taxon.names){
        nominal.fit <- vglm(prim.var ~., data = data[,c("prim.var", covariate, taxon)], family = multinomial)
        nominal.lrt.summary <- lrtest(nominal.null, nominal.fit)
        nominal.global.summary.df <- rbind(nominal.global.summary.df, data.frame(rank = str_to_title(name), taxon = taxon, subset(nominal.lrt.summary@Body, select = -1)[2,]))
      }
      nominal.global.summary.df$ID <- 1:nrow(nominal.global.summary.df)
      nominal.global.summary.df$Qval <- p.adjust(nominal.global.summary.df$Pr..Chisq., "BH")
      nominal.global.summary.df <- nominal.global.summary.df[,c(7, 1:6, 8)]
      colnames(nominal.global.summary.df) <- c("ID", "Rank", "Taxon", "Log Likelihood", "DF", "Chisq", "P-value", "Q-value")
      rownames(nominal.global.summary.df) <- NULL
      
      taxa.nominal.global.summary.list[[name]] <- nominal.global.summary.df
    }
  }
  return(taxa.nominal.global.summary.list)
}

taxa.nominal.global.outcome <- function(taxa.nominal.global.summary.list, level.names, res = 0.05){
  taxa.global.table <- data.frame()
  
  for(name in level.names){
    taxa.global.table <- rbind(taxa.global.table, taxa.nominal.global.summary.list[[name]])
  }
  
  ind <- which(taxa.global.table$`Q-value` < res)
  taxa.global.table <- taxa.global.table[ind,]
  rownames(taxa.global.table) <- NULL
  taxa.global.table$ID <- rownames(taxa.global.table)
  taxa.global.table$`Log Likelihood` <- round(taxa.global.table$`Log Likelihood`, 3)
  taxa.global.table$Chisq <- round(taxa.global.table$Chisq, 3)
  taxa.global.table$DF <- abs(taxa.global.table$DF)
  
  taxa.global.table$`P-value` <- ifelse(round(taxa.global.table$`P-value`, 3) == 0, "<.001", sprintf("%.3f", round(taxa.global.table$`P-value`, 3)))
  taxa.global.table$`Q-value` <- ifelse(round(taxa.global.table$`Q-value`, 3) == 0, "<.001", sprintf("%.3f", round(taxa.global.table$`Q-value`, 3)))
  
  taxa.global.table <- subset(taxa.global.table, select = -ID)
  
  return(taxa.global.table)
}
