gen_hybrid = function(H1, H2, alpha, jumps){
  nHAP <- min(ncol(H1), ncol(H2))
  if (nrow(H1) != nrow(H2)){
    stop("Ancestral populations should contain the same number of markers.")
  }
  nSNP <- nrow(H1)
  n.jumps <- sum(jumps)
  n.chunks <- n.jumps + 1 
  beg <- vector(mode = "numeric", length = n.chunks)
  end <- vector(mode = "numeric", length = n.chunks)
  beg[1] <- 1
  end[n.chunks] <- nSNP
  jumps.loc <- which(jumps == 1)
  beg[-1] <- jumps.loc
  end[1:n.jumps] <- pmax(jumps.loc - 1, 1)
  idx.father <- sample(1:nHAP, size = n.chunks, replace = TRUE)
  idx.mother <- sample(1:nHAP, size = n.chunks, replace = TRUE)
  haplotype.1 <- vector(mode = "numeric", length = nSNP)
  haplotype.2 <- vector(mode = "numeric", length = nSNP)
  for (i in 1:n.chunks){
    p <- 1 - alpha
    nbino <- rbinom(1, 2, prob = p)
    if (nbino == 2){
      haplotype.1[beg[i]:end[i]] <- H1[beg[i]:end[i], idx.father[i]]
      haplotype.2[beg[i]:end[i]] <- H1[beg[i]:end[i], idx.mother[i]]
    } else if (nbino == 1){
      haplotype.1[beg[i]:end[i]] <- H1[beg[i]:end[i], idx.father[i]]
      haplotype.2[beg[i]:end[i]] <- H2[beg[i]:end[i], idx.mother[i]]
    } else if (nbino == 0){
      haplotype.1[beg[i]:end[i]] <- H2[beg[i]:end[i], idx.father[i]]
      haplotype.2[beg[i]:end[i]] <- H2[beg[i]:end[i], idx.mother[i]]
    }
  }
  return(list(h1 = haplotype.1, h2 = haplotype.2))
}

gen_matrix = function(H1, H2, alpha = 0.5, gen_map, n.hyb = ncol(H1) / 2, lambda = 1.0){
  H <- matrix(0, nrow = nrow(H1), ncol = (2 * n.hyb))
  for (i in 1:n.hyb){
    jumps <- jumps_from_map(gen_map, lambda = lambda)  
    h <- gen_hybrid(H1, H2, alpha, jumps)
    H[, (2 * i - 1)] <- h$h1
    H[, (2 * i)] <- h$h2
  }
  return(H)
}

create_input_pcadapt = function(H1, H2, H){
  G1 <- haplo_to_geno(H1)  
  G2 <- haplo_to_geno(H2)
  G <- haplo_to_geno(H)
  nIND <- ncol(G1) + ncol(G2) + ncol(G)
  Gmat <- matrix(0, nrow = nrow(G), ncol = nIND)
  Gmat[, 1:ncol(G1)] <- G1
  Gmat[, (ncol(G1) + 1):(ncol(G1) + ncol(G2))] <- G2
  Gmat[, (ncol(G1) + ncol(G2) + 1):(ncol(G1) + ncol(G2) + ncol(G))] <- G
  return(Gmat)
}

create_input_hapmix = function(H1, H2, H, gen_map, phys_map, chr, 
                               ancstrl.1 = "ANC1", ancstrl.2 = "ANC2", 
                               admxd = "AA", rates.digits = 8){
  cat("Creating the genotype file for ancestral population 1...")
  write.table(t(H1), paste0(ancstrl.1, ".lfmm"), col.names = FALSE, row.names = FALSE)
  cat("DONE\n")
  LEA::lfmm2geno(paste0(ancstrl.1, ".lfmm"), paste0(ancstrl.1, ".geno"), force = TRUE)
  cat("Creating the genotype file for ancestral population 2...")
  write.table(t(H2), paste0(ancstrl.2, ".lfmm"), col.names = FALSE, row.names = FALSE)
  cat("DONE\n")
  LEA::lfmm2geno(paste0(ancstrl.2, ".lfmm"), paste0(ancstrl.2, ".geno"), force = TRUE)
  cat("Creating the genotype file for the admixed population...")
  G <- haplo_to_geno(H)
  write.table(t(G), paste0(admxd, ".lfmm"), col.names = FALSE, row.names = FALSE)
  cat("DONE\n")
  LEA::lfmm2geno(paste0(admxd, ".lfmm"), paste0(admxd, ".geno"), force = TRUE)
  
  dt <- array(dim = c(nrow(H1), ncol = 4))
  dt[, 1] <- paste0("rs", 1:nrow(H1))
  dt[, 2] <- chr
  dt[, 3] <- format(round(gen_map, rates.digits), nsmall = rates.digits)
  dt[, 4] <- phys_map
  cat("Creating the snpfiles...")
  write.table(dt, paste0(ancstrl.1, "snpfile"), col.names = FALSE, row.names = FALSE, quote = FALSE)
  write.table(dt, paste0(ancstrl.2, "snpfile"), col.names = FALSE, row.names = FALSE, quote = FALSE)
  write.table(dt, paste0(admxd, "snpfile"), col.names = FALSE, row.names = FALSE, quote = FALSE)
  dt.ind <- as.data.frame(matrix(nrow = ncol(G), ncol = 3))
  dt.ind[[1]] <- paste0(rep("IND_", ncol(G)), 1:ncol(G))
  dt.ind[[2]] <- "F"
  dt.ind[[3]] <- "Hybrid"
  write.table(dt.ind, paste0(admxd, ".ind"), col.names = FALSE, row.names = FALSE, quote = FALSE)
  cat("DONE\n")
  cat("Creating the rates file...")
  l1 <- paste0(":sites:", nrow(G))
  l2 <- stringr::str_c(phys_map, collapse = " ")
  l3 <- stringr::str_c(format(gen_map, scientific = FALSE, nsmall = rates.digits), 
                       collapse = " ")
  con <- file("rates")
  writeLines(paste0(l1, "\n", l2, "\n", l3), con = con)
  close(con = con)
  cat("DONE\n")
}