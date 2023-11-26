library(phyloseq)
library(proxy)
library(ecodist)
library(GUniFrac)
library(MiRKAT)
library(MiRKATMC)
library(vegan)

source("Source/MiDataProc.Data.Input.R")
source("Source/MiDataProc.Data.Upload.R")

# MiRKATMC ----------------------
Ds.Ks.func <- function(rare.biom, biom.after.qc, is.tree) {
  
  rare.otu.tab <- otu_table(rare.biom)
  no.rare.otu.tab <- otu_table(biom.after.qc)
  
  if (is.tree == "withTree") {
    
    no.rare.tree <- phy_tree(biom.after.qc)
    
    jac <- as.matrix(proxy::dist(t(rare.otu.tab), method = "Jaccard"))
    bc <- as.matrix(bcdist(t(rare.otu.tab)))
    unifs <- GUniFrac(t(no.rare.otu.tab ), no.rare.tree, alpha = c(0.5, 1))$unifracs
    u.unif <- unifs[, , "d_UW"]
    g.unif <- unifs[, , "d_0.5"]
    w.unif <- unifs[, , "d_1"]
    
    jac.k <- D2K(jac)
    bc.k <- D2K(bc)
    u.unif.k <- D2K(u.unif)
    g.unif.k <- D2K(g.unif)
    w.unif.k <- D2K(w.unif)
    
    rownames(jac.k) <- colnames(jac.k) <- colnames(rare.otu.tab)
    rownames(bc.k) <- colnames(bc.k) <- colnames(rare.otu.tab)
    rownames(u.unif.k) <- colnames(u.unif.k) <- colnames(rare.otu.tab)
    rownames(g.unif.k) <- colnames(g.unif.k) <- colnames(rare.otu.tab)
    rownames(w.unif.k) <- colnames(w.unif.k) <- colnames(rare.otu.tab)
    
    return(list(Ds = list(Jaccard = jac, Bray.Curtis = bc, U.UniFrac = u.unif, G.UniFrac = g.unif, W.UniFrac = w.unif),
                Ks = list(Jaccard = jac.k, Bray.Curtis = bc.k, U.UniFrac = u.unif.k, G.UniFrac = g.unif.k, W.UniFrac = w.unif.k)))
    
  } 
  else if (is.tree == "withoutTree") {
    
    jac <- as.matrix(proxy::dist(t(rare.otu.tab), method = "Jaccard"))
    bc <- as.matrix(bcdist(t(rare.otu.tab)))
    
    jac.k <- D2K(jac)
    bc.k <- D2K(bc)
    
    rownames(jac.k) <- colnames(jac.k) <- colnames(rare.otu.tab)
    rownames(bc.k) <- colnames(bc.k) <- colnames(rare.otu.tab)
    
    return(list(Ds = list(Jaccard = jac, Bray.Curtis = bc), Ks = list(Jaccard = jac.k, Bray.Curtis = bc.k)))
  }
}

p.value.0.1 <- function(x, round.x = 3) {
  x <- format(round(x, 3), digits = 3)
  ind.0 <- which(x == "0.000" | x == 0)
  x[ind.0] <- "<.001"
  ind.1 <- which(x == "1.000" | x == 1)
  x[ind.1] <- ">.999"
  return(x)
}

beta.PCoA.plot <- function(Ds, pvs, y){
  par(mfrow = c(3, 2))
  for (i in 1:length(Ds)) {
    if (pvs[i] < 0.05) { # Error in if: argument is of length zero pvs 확인 필요.
      sub.tit <- paste("*p:", p.value.0.1(pvs[i]), sep="")
    }
    if (pvs[i] >= 0.05) {
      sub.tit <- paste("p:", p.value.0.1(pvs[i]), sep="")
    }
    mod <- betadisper(as.dist(Ds[[i]]), y)
    plot(mod, ellipse = TRUE, hull = FALSE, main = names(Ds)[i], xlab="PC 1", ylab="PC 2",
         sub = sub.tit, col = 1:nlevels(y), cex=1.5)
  }
  
  plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
  legend("center", title = NULL, legend = levels(y), fil = c(1:nlevels(y), cex=2.5, box.lty=0), 
         bty = "n", cex=1.5)
}

reduced_beta <- function(beta.dat, sam.dat, prim_id){
  
  mat_vec <- match(rownames(beta.dat[[1]]), rownames(beta.dat[[1]]))
  dat <- list() 
  if(length(mat_vec[mat_vec == FALSE]) != 0){
    
    for (i in 1:length(beta.dat)){
      ind <- match(rownames(beta.dat[[i]], rownames(sam.dat)))
      dat[[i]] <- beta.dat[[i]][ind, ind]
    }
  }else{
    for (i in 1:length(beta.dat)){
      dat[[i]] <- beta.dat[[i]]
    }
  }
  names(dat) <- names(beta.dat)
  return(dat)
}

beta.permanova.multi <- function(num_perm = 3000, div_num, beta_div, sam_dat, prim_id, method_adj = "BH"){
  
  beta_list <- list() 
  
  ind <- list() 
  time.p.cat <- names(table(sam_dat[,prim_id]))
  for (i in 1:length(time.p.cat)){
    ind[[i]] <- which(sam_dat[,prim_id] == time.p.cat[i])
  }
  
  comb <- combn(length(time.p.cat), 2)
  out <- matrix(NA, ncol(comb), 5)
  out_2 <- matrix(NA, ncol(comb), 5)
  beta_div_ind <- beta_div[[div_num]]
  
  for (j in 1:ncol(comb)){
    combination <- comb[,j]
    combined_ind <- c(ind[[combination[[1]]]], ind[[combination[[2]]]])
    
    beta.dat.pair <- beta_div_ind[combined_ind, combined_ind]
    sam.dat.pair <- sam_dat[combined_ind,]
    
    perm <- how(nperm = num_perm)
    
    fit <- adonis2(beta.dat.pair ~ sam.dat.pair, data = data.frame(sam.dat.pair = sam.dat.pair[[prim_id]]), permutations = perm)
    
    out[j,] <- c(time.p.cat[combination[1]], time.p.cat[combination[2]], fit[1, "F"], fit[1,"R2"], fit[1,"Pr(>F)"])
    out_2[j,] <- c(time.p.cat[combination[1]], time.p.cat[combination[2]], decimal_adjust(fit[1, "F"]), decimal_adjust(fit[1,"R2"]), fit[1,"Pr(>F)"])
  }
  
  q_val <- p.adjust(out[,5], method = method_adj)
  out <- as.data.frame(cbind(out[,1:3], out[,5], q_val))
  
  q_val_2 <- p.value.0.1_char(decimal_adjust(p.adjust(out_2[,5], method = method_adj, n = length(out_2[,5]))))  
  out_2 <- as.data.frame(cbind(out_2[,1:3],  q_val_2))
  
  colnames(out) <- c("Ref", "Com", "F", "P.value", "Adj. P.value")
  colnames(out_2) <- c("Ref", "Com", "F", "Adj. P.value")
  
  beta_list$download <- out
  beta_list$table <- out_2 
  
  return(beta_list)
}

# beta_div <- dat_1
# sam_dat <- dat_2 
# prim_id <- dat_3 
# num_perm = 3000 
# method_adj = "BH"
# download = FALSE

permanova.mult.united <- function(num_perm = 3000, beta_div, sam_dat, prim_id, method_adj = "BH", download = FALSE){
  total_list <- list()
  download_list <- list() 
  table_list <- list()
  
  for (i in 1:length(beta_div)){
    div_list <- beta.permanova.multi(num_perm = 10000, i, beta_div, sam_dat, prim_id, method_adj)
    download_list[[i]] <- div_list$download 
    table_list[[i]] <- div_list$table
  }
  
  names(download_list) <- names(beta_div) 
  names(table_list) <- names(beta_div)
  
  total_list$download <- download_list
  total_list$table <- table_list
  return(total_list)
}

beta.permanova.mult.glob <- function(num_perm = 3000, beta_div, sam_dat, prim_id, download = FALSE){
  
  global_result <- list() 
  p_val_result <- c() 
  
  out <- matrix(NA, length(beta_div), 3)
  
  for (i in 1:length(names(beta_div))){
    perm <- how(nperm = num_perm)
    beta_div_ind <- beta_div[[i]]
    
    fit <- adonis2(beta_div_ind ~ sam_dat_pair, data = data.frame(sam_dat_pair = sam_dat[[prim_id]]), permutations = perm)
    
    if (download){
      out[i,] <- c(fit[1, "F"], fit[1,"R2"], fit[1,"Pr(>F)"])
      each_p_result <- decimal_adjust(fit[1, "Pr(>F)"])
      p_val_result <- c(p_val_result, each_p_result)
    }else{
      out[i,] <- c(decimal_adjust(fit[1, "F"]), decimal_adjust(fit[1,"R2"]), decimal_adjust(fit[1,"Pr(>F)"]))
      each_p_result <- decimal_adjust(fit[1, "Pr(>F)"])
      p_val_result <- c(p_val_result, each_p_result)
    }
  }
  
  colnames(out) <- c("F", "R2", "P.value")
  rownames(out) <- names(beta_div)
  global_result$table <- out
  global_result$p.val <- p_val_result 
  
  return(global_result)
}
