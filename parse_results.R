setwd("~/Documents/thesis/git/Introgression/HapmixRelease/populus/RUN")
require(RcppRoll)

chr <- 6
filename <- "HYBRID.LOCALANC."
nIND <- 25

anc.i = function(loc.i){
  res <- which(loc.i == max(loc.i)) - 1  
  return(2 - res)
}

ancestries = function(loc){
  apply(loc, MARGIN = 1, FUN = anc.i)
}

loc.anc <- read.table(paste0(filename, 0, ".", chr))
total.anc <- ancestries(loc.anc) / (2 * nIND)

for (k in 1:(nIND - 1)){
  loc.anc <- read.table(paste0(filename, k, ".", chr))
  total.anc <- total.anc + ancestries(loc.anc) / (2 * nIND)
}

stat.hapmix <- roll_mean(total.anc, n = 100)
plot(stat.hapmix)

G <- create_input_pcadapt(H1, H2, H3)
x <- pcadapt(G, K = 2, min.maf = 0.001)

f <- scan.intro(G, K = 1, pop = pop, ancstrl.1 = "Trichocarpa", ancstrl.2 = "Balsamifera",
                admxd = "Hybrid", min.maf = 0.001, window.size = 1000)

class(f) <- NULL

