demo_geno <- read.delim(file = "haplotype_dataset_2.txt",
                        header = TRUE, stringsAsFactors = TRUE)
geo_dist_km <- demo_geno[, "Geo.dist.km."]
mean_haplo_het <- demo_geno[, "MeanHaploHet"]
region <- demo_geno[, "Region"]
region <- region[!is.na(geo_dist_km)]
geo_dist_km <- geo_dist_km[!is.na(geo_dist_km)]
mean_haplo_het <- mean_haplo_het[!is.na(mean_haplo_het)]
selectTF <- function (x,y) {
  o <- {}
  for (i in 1:length(x)) {
    if (x[i] == y) {
      o[i] <- TRUE
    }
    else {
      o[i] <- FALSE
    }
  }
  return(o)
}
selectIndex <- function (x,y) {
  o <- {}
  for (i in 1:length(x)) {
    if (x[i] == y) {
      o <- append(o,i)
    }
  }
  return(o)
}
##
geo_dist_km[selectTF(region,"Central/South Asia")]
geo_dist_km[region %in% "Central/South Asia"]
##
geo_dist_km[selectIndex(region,"Africa")]
geo_dist_km[which(region=="Africa")]

### Calculate Pearson's correlation coefficient ===========================
#rS1vsT.wt <- rS1vsT.wt / sum(rS1vsT.wt)
rhoPearson.weighted <- function (x, y, w) {
  stopifnot(length(x) == dim(y)[2] )
  w <- w / sum(w)
  # Center x and y, using the weighted means
  x <- x - sum(x * w)
  ty <- t(y - colSums(t(y) * w))
  # Compute the variance
  vx <- sum(w * x * x)
  vy <- colSums(w * ty * ty)
  # Compute the covariance
  vxy <- colSums(ty * x * w)
  # Compute the correlation
  return(vxy / sqrt(vx * vy))
}
cor.test(head(dataFrame$temp, -1), head(dataFrame$rSpecies_1, -1), method = "pearson")
cor.test(dataFrame.fit$temp, dataFrame.fit$rSpecies_1, method = "pearson")
rhoPearson.weighted(head(pop_data$temp, -1), head(pop_data$rSpecies_1, -1), rS1vsT.wt)
