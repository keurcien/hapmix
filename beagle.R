# Extract haplotypes from Beagle outputs

library(vcfR)

setwd("~/Documents/thesis/git/Introgression/populus/Original_data_set_ch6_12_15/")

# Read Beagle output with vcfR
vcf <- read.vcfR("test_beagle.vcf") 
pop <- read.table("populus_original.txt", header = TRUE)$POP 

recombinationRate <- 0.05 
anc1 <- 1
anc2 <- 3

# Duplicate labels
hap.pop <- vector(length = 2 * length(pop), mode = "numeric")
hap.pop[seq(1, length(hap.pop), by = 2)] <- pop
hap.pop[seq(2, length(hap.pop), by = 2)] <- pop

# Use vcfR to extract haplotypes, positions and ref/alt alleles
pos <- getPOS(vcf)
hap <- extract.haps(vcf)
ref <- getREF(vcf)
alt <- getALT(vcf)

# Convert positions to genetic distance
gen <- recombinationRate * pos / 1000000
rates <- rbind(pos, gen)

# Separate ancestral haplotypes
nSNP <- nrow(hap)
nIND <- ncol(hap) / 2
hap.1 <- hap[, (hap.pop == anc1)]
hap.2 <- hap[, (hap.pop == anc2)]
hap.int.1 <- array(0, dim = dim(hap.1))
hap.int.2 <- array(0, dim = dim(hap.2))

for (j in 1:ncol(hap.1)){
  hap.int.1[, j] <- 1 * (alt == hap.1[, j])
}

for (j in 1:ncol(hap.2)){
  hap.int.2[, j] <- 1 * (alt == hap.2[, j])
}

# Convert ancestral haplotypes to genotypes
geno.int.1 <- array(0, dim = c(nrow(hap.int.1), ncol(hap.int.1) / 2))
geno.int.2 <- array(0, dim = c(nrow(hap.int.2), ncol(hap.int.2) / 2))

for (j in 1:ncol(geno.int.1)){
  geno.int.1[, j] <- hap.int.1[, (2 * j -1)] + hap.int.1[, (2 * j)]
}

for (j in 1:ncol(geno.int.2)){
  geno.int.2[, j] <- hap.int.2[, (2 * j -1)] + hap.int.2[, (2 * j)]
}

write.table(rates, "populus.map", col.names = FALSE, row.names = FALSE)
write.table(hap.int.1, "H_1", col.names = FALSE, row.names = FALSE)
write.table(hap.int.2, "H_2", col.names = FALSE, row.names = FALSE)
write.table(geno.int.1, "G_1", col.names = FALSE, row.names = FALSE)
write.table(geno.int.2, "G_2", col.names = FALSE, row.names = FALSE)
