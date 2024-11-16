#!/usr/bin/env RScript
# chmod +x Script.R
rm(list = ls(all = TRUE)) # Remove previous objects from globe environment
tryCatch(assign("last.warning", NULL, envir = baseenv()), error = function(e) {}) # Remove previous warnings
tryCatch(dev.off(dev.list()["RStudioGD"]), error = function(e) {}) # Remove previous plots in RStudio
loadPackage <- function(x, ...) {
  x <- c(x, ...)
  un <- x[!(x %in% installed.packages()[, "Package"])]
  if (length(un)) {
    invisible(readline(prompt = paste("Press [Enter] to install package \"", paste(un, collapse = ", "), "\"", sep = "")))
    install.packages(un, dependencies = TRUE)
  }
  sapply(x, require, character.only = TRUE)
}
getmCq <- function(dataFrame, threshold = 0.7, c = 2, outlier = TRUE, outlierTest = grubbs.test, outlier.omit = TRUE) {
  # dataFrame$Cq <- lapply(dataFrame$Cq, function(x) {gsub("Undetermined", "0", x)})
  target <- dataFrame$Target[!duplicated(dataFrame$Target)]
  sample <- dataFrame$Sample[!duplicated(dataFrame$Sample)]
  mCq <- matrix(nrow = length(sample), ncol = length(target))
  rownames(mCq) <- sample
  colnames(mCq) <- target
  for (i in 1:length(target)) {
    for (j in 1:length(sample)) {
      Cq <- dataFrame$Cq[dataFrame$Target == target[i] & dataFrame$Sample == sample[j]]
      if (any(grepl("Undetermined", Cq))) {
        Cq <- Cq[-which(Cq == "Undetermined")]
      }
      Cq <- as.numeric(Cq)
      if ((!(max(Cq) - min(Cq) == 0)) && (length(Cq) > 0)) {
        Cq <- as.numeric(Cq)
        if (outlier) {
          t <- outlierTest(Cq)
          p <- t$p.value
          if (p < 0.05) {
            if (strsplit(t$alternative, "\\s")[[1]][1] == "highest") {
              o <- max(Cq)
            } else {
              o <- min(Cq)
            }
            w <- which(Cq == o)
            if (outlier.omit) {
              Cq <- Cq[-w]
              o <- "Omitted"
            }
            warning(paste(target[i], sample[j], w, "is an outlier,", "p =", round(t$p.value, digits = 3), o))
          }
        }
        a <- c()
        p <- c()
        for (k in 1:length(Cq)) {
          o <- 0
          for (l in 1:(length(Cq) - 1)) {
            diff <- abs(Cq[k] - Cq[-k][l])
            if (diff > threshold) {
              o <- o + 1
            }
            if ((!any(p == Cq[-k][l])) & (diff > threshold)) {
              p[length(p) + 1] <- Cq[k]
              warning(paste(target[i], sample[j], k, "v.s.", which(Cq == Cq[-k][l]), "Cq diff =", round(diff, digits = 5)))
            }
          }
          if (o == c) {
            warning(paste(target[i], sample[j], k, "Omitted"))
          } else {
            a[length(a) + 1] <- Cq[k]
          }
        }
        mCq[j, i] <- mean(a)
      }
    }
  }
  return(mCq)
}
getdCq <- function(mCq, e = 1) { # e: endogenous reference
  dCq <- mCq - mCq[, 1]
  return(dCq[-nrow(dCq), -1])
}
getddCq <- function(dCq, c = 1, delete.c = FALSE) {
  for (i in 1:(nrow(dCq))) {
    if (i != c) {
      dCq[i, ] <- dCq[i, ] - dCq[c, ]
    }
  }
  if (delete.c) {
    dCq <- dCq[-c]
  } else {
    dCq[c, ] <- dCq[c, ] - dCq[c, ]
  }
  return(dCq)
}
getdf <- function(matrix) {
  df <- as.data.frame(matrix)
  df <- cbind(sample = rownames(df), df)
  df <- melt(df, id = "sample")
  sample <- df$sample[!duplicated(df$sample)]
  df$group <- df$sample
  for (i in 1:length(sample)) {
    df$sample[df$sample == sample[i]] <- i
  }
  return(df)
}
plotGroupedBarChart <- function(df, label = "default", omit = c(), y.action = NULL) {
  colour <- c("#00a3ff", "#5dd934", "#f9bb00", "#ff2501", "#d41876", "#929292", "#b0643c", "#2994ac", "#3f7cad")
  l <- length(omit)
  if (!l == 0) {
    colour <- colour[-omit]
    for (i in 1:length(omit)) {
      df <- df[!df$sample == as.character(i), ]
    }
  }
  group <- df$group[!duplicated(df$group)]
  colour <- colour[1:length(group)]
  if (length(label) == length(group)) {
    group <- label
  }
  if (!is.null(y.action)) {
    func <- get(y.action)
    ggplot(df, aes(x = variable, y = func(value), fill = sample)) +
      geom_bar(stat = "identity", position = "dodge") +
      scale_fill_manual(name = "Treatment", label = group, values = colour) +
      xlab("Gene") +
      ylab(paste("Relative fold change (", y.action, ")", sep = "")) +
      theme_classic()
  } else {
    ggplot(df, aes(x = variable, y = value, fill = sample)) +
      geom_bar(stat = "identity", position = "dodge") +
      scale_fill_manual(name = "Treatment", label = group, values = colour) +
      xlab("Gene") +
      ylab("Relative fold change") +
      theme_claxssic()
  }
}
loadPackage("readxl", "ggplot2", "reshape2", "outliers")
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
path.1 <- "R1 repeat RA BMS DEAB-sfrp1a nr2f5 cyp26b1-2023-07-25-15205_20230725_142259_20230725 151152.xls"
dataFrame.1 <- read_excel(path.1, sheet = "Results", range = cell_limits(c(23, 1), c(NA, NA)))
mCq.1 <- getmCq(dataFrame.1)
dCq.1 <- getdCq(mCq.1)
ddCq.1 <- getddCq(dCq.1)
ddCq.DMSO.1 <- getddCq(dCq.1, c = 2)
crtl.1 <- 2^-ddCq.1
DMSO.1 <- 2^-ddCq.DMSO.1
plotGroupedBarChart(getdf(DMSO.1), c("Contol", "RA", "BMS961", "DEAB"), y.action = "log2", omit = c(1, 2))
plotGroupedBarChart(getdf(crtl.1), c("Contol", "DMSO", "RA", "BMS961", "DEAB"), y.action = "log2", omit = 1)
path.2 <- "Run 2 RA BMS965 DEAB sfrp1a nr2f5 cyp26b1_2023-07-27-17405_20230727_164251_20230727 173143.xls"
dataFrame.2 <- read_excel(path.2, sheet = "Results", range = cell_limits(c(23, 1), c(NA, NA)))
mCq.2 <- getmCq(dataFrame.2)
dCq.2 <- getdCq(mCq.2)
ddCq.2 <- getddCq(dCq.2)
ddCq.DMSO.2 <- getddCq(dCq.2, c = 2)
crtl.2 <- 2^-ddCq.2
DMSO.2 <- 2^-ddCq.DMSO.2
plotGroupedBarChart(getdf(DMSO.2), c("Contol", "RA", "BMS961", "DEAB"), y.action = "log2", omit = c(1, 2))
plotGroupedBarChart(getdf(crtl.2), c("Contol", "DMSO", "RA", "BMS961", "DEAB"), y.action = "log2", omit = 1)
path.3 <- "Repeat R3 RA BMS961 DEAB sfrp1a nr2f5 cyp26b1_2023-08-03-161020_20230803_151243_20230803 160122.xls"
dataFrame.3 <- read_excel(path.3, sheet = "Results", range = cell_limits(c(23, 1), c(NA, NA)))
mCq.3 <- getmCq(dataFrame.3, outlier = FALSE)
dCq.3 <- getdCq(mCq.3)
ddCq.3 <- getddCq(dCq.3)
ddCq.DMSO.3 <- getddCq(dCq.3, c = 2)
crtl.3 <- 2^-ddCq.3
DMSO.3 <- 2^-ddCq.DMSO.3
plotGroupedBarChart(getdf(DMSO.3), c("Contol", "RA", "BMS961", "DEAB"), y.action = "log2", omit = c(1, 2))
plotGroupedBarChart(getdf(crtl.3), c("Contol", "DMSO", "RA", "BMS961", "DEAB"), y.action = "log2", omit = 1)
path.4 <- "R4 RA BMS961 DEAB sfrp1a nr2f5 cyp26b1_2023-08-03-172323_20230803_162544_20230803 171420.xls"
dataFrame.4 <- read_excel(path.4, sheet = "Results", range = cell_limits(c(23, 1), c(NA, NA)))
mCq.4 <- getmCq(dataFrame.4, outlier = FALSE)
dCq.4 <- getdCq(mCq.4)
ddCq.4 <- getddCq(dCq.4)
ddCq.DMSO.4 <- getddCq(dCq.4, c = 2)
crtl.4 <- 2^-ddCq.4
DMSO.4 <- 2^-ddCq.DMSO.4
plotGroupedBarChart(getdf(DMSO.4), c("Contol", "RA", "BMS961", "DEAB"), y.action = "log2", omit = c(1, 2))
plotGroupedBarChart(getdf(crtl.4), c("Contol", "DMSO", "RA", "BMS961", "DEAB"), y.action = "log2", omit = 1)
crtl.1
DMSO.1
crtl.2
DMSO.2
crtl.3
DMSO.3
crtl.4
DMSO.4
