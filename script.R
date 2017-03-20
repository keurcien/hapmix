library(vcfR)
library(data.table)
library(Rcpp)

setwd("~/Documents/thesis/git/Introgression/populus/Original_data_set_ch6_12_15/")
sourceCpp("~/Documents/thesis/git/hapmix/simulate.cpp")
source("~/Documents/thesis/git/hapmix/simulate.R")

gen_map <- t(as.matrix(fread("populus.map")))
H1 <- as.matrix(fread("H_1"))
H2 <- as.matrix(fread("H_2"))

nSNP <- 30000

H1 <- H1[1:nSNP, ]
H2 <- H2[1:nSNP, ]
H3 <- gen_matrix(H1, H2, alpha = 0.3, n.hyb = 50, gen_map = gen_map[1:nSNP, 2], lambda = 10)

create_input_hapmix(H1 = H1, H2 = H2, H = H3, 
                    gen_map = gen_map[1:nSNP, 2], 
                    phys_map = gen_map[1:nSNP, 1], chr = 6, 
                    ancstrl.1 = "ANC1", ancstrl.2 = "ANC2", 
                    admxd = "AA", rates.digits = 8)

G <- create_input_pcadapt(H1, H2, H3)
pop <- c(rep("ANC1", 25), rep("ANC2", 25), rep("Hybrid", 50))
x <- pcadapt(G, K = 2)
plot(x, option = "scores", pop = pop)