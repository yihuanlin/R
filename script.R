#!/usr/bin/env RScript
# §1. Clear workspace and define data frame and public functions ----------
rm(list = ls(all = TRUE)) # Remove previous objects from globe environment
tryCatch(assign("last.warning", NULL, envir = baseenv()), error = function(e){}) # Remove previous warnings
tryCatch(dev.off(dev.list()["RStudioGD"]), error = function(e){}) # Remove previous plots in RStudio
loadPackage <- function(x, ...) {
  x <- c(x, ...)
  un <- x[!(x %in% installed.packages()[, "Package"])]
  if (length(un)) {
    invisible(readline(prompt = paste("Press [Enter] to install package \"", paste(un, collapse = ", "), "\"", sep = "")))
    install.packages(un, dependencies = TRUE)
  }
  sapply(x, require, character.only = TRUE)
}
loadPackage("ggplot2", "car", "sandwich", "lmtest", "caret")
# Require "ggplot2" for ggplot(), "car" for ncvTest(), vif(), outlierTest() and avPlots(), "sandwich" for vcovHC(), "lmtest" for bgtest() and coeftest(), "caret" for train()
getADNormarlityTest <- function (x, warning = TRUE) {
  # Modified from ad.test of nortest package "Tests for Normality", GPLv3
  # Author: Juergen Gross [aut], Uwe Ligges [aut, cre]
  x.name <- deparse(substitute(x))
  x <- sort(x[complete.cases(x)])
  n <- length(x)
  if (n < 8) 
    stop("Sample size must be greater than 7")
  logp1 <- pnorm((x-mean(x)) / sd(x), log.p = TRUE)
  logp2 <- pnorm(-(x-mean(x)) / sd(x), log.p = TRUE)
  h <-  (2 * seq(1:n) - 1) * (logp1 + rev(logp2)) 
  A <- -n - mean(h)
  AA <- (1 + 0.75 / n + 2.25 / n^2) * A
  if (AA < 0.2) {
    p <- 1 - exp(-13.436 + 101.14 * AA - 223.73 * AA^2)
  }
  else if (AA < 0.34) {
    p <- 1 - exp(-8.318 + 42.796 * AA - 59.938 * AA^2)
  }
  else if (AA < 0.6) {
    p <- exp(0.9177 - 4.279 * AA - 1.38 * AA^2)
  }
  else if (AA < 10) {
    p <- exp(1.2937 - 5.709 * AA + 0.0186 * AA^2)
  } 
  else p <- 3.7e-24
  if (p < 0.05&&warning) warning(x.name, " does not follow a normal distribution, p = ", sprintf("%.3f", p))
  result <- list(statistic = c(A = A), p.value = p, method = "Anderson-Darling normality test", data.name = x.name)
  class(result) <- "htest"
  return(result)
}
linearModel.getDeltaRsquare <- function (model) {
  # Modified from getDeltaRsquare of rockchalk package, GPLv3
  # Author: Paul E. Johnson [aut, cre], Gabor Grothendieck [ctb]
  modeldrop1 <- drop1(model)
  RSS <- modeldrop1[1, "RSS"]
  deltaSS <- modeldrop1[, "Sum of Sq"]
  SST = sum((model$model[, 1] - mean(model$model[, 1]))^2)
  deltaRsquare <- deltaSS/SST
  names(deltaRsquare) <- row.names(modeldrop1)
  omit <- is.na(deltaRsquare)
  deltaRsquare <- deltaRsquare[-omit]
  cat("The deltaR-square values: the change in the R-square\n      observed when a single term is removed.", 
      fill = TRUE)
  cat("Same as the square of the 'semi-partial correlation coefficient'", 
      fill = TRUE)
  print(as.data.frame(deltaRsquare))
  invisible(deltaRsquare)
}
dataFrame.getADNormarlityTestTable <- function (df, col, ..., warning = TRUE) {
  x <- c(col, ...)
  e <- matrix(ncol = 2, nrow = length(x))
  for (i in 1:length(x)) {
    r <- getADNormarlityTest(df[,colnames(df) == x[i]], warning)
    e[i,1] <- r$statistic
    e[i,2] <- r$p.value
  }
  rownames(e) <- x
  colnames(e) <- c("A", "p")
  return(e)
}
getConfidenceInterval <- function(x, SE, df) {
  ME <- qt(0.975, df) * SE # Get 95% confidence interval by default
  return(paste(
    sprintf("%.4f", x), " ± ", # Round to 3 decimal places
    sprintf("%.4f", SE), "; from ", # Return ± standard error,
    sprintf("%.4f", x - ME), " to ",
    sprintf("%.4f", x + ME), sep="" # and 95% confidence interval
  ))
}
mean.printConfidenceInterval <- function(x) {
  x <- x[complete.cases(x)]
  N <- length(x)
  mean <- mean(x)
  SE <- sd(x) / sqrt(N)
  df <- N - 1
  cat(paste("Mean ± SE:", getConfidenceInterval(mean, SE, df)))
}
mean.getConfidenceInterval <- function(x) {
  capture.output(mean.printConfidenceInterval(x))
}
querySelector <- function (char = "", inl = NULL, exl = NULL) {
  ls <- ls(pattern = char, envir = .GlobalEnv)
  if(is.null(inl)) {
    if(!is.null(exl)) {
      for (i in 1:length(exl)) {
        ls <- ls[-grep(exl[i], ls)]
      }
    }
  } else if(is.null(exl)) {
    for (i in 1:length(inl)) {
      ls <- ls[grep(inl[i], ls)]
    }
  } else {
    for (i in 1:length(inl)) {
      ls <- ls[grep(inl[i], ls)]
    }
    for (i in 1:length(exl)) {
      ls <- ls[-grep(exl[i], ls)]
    }
  }
  return(ls)
}
linearModel.getADNormarlityTestTable.querySelector <- function (char = "", ..., warning = TRUE) {
  ls <- querySelector(char, ...)
  e <- matrix(ncol = 2, nrow = length(ls))
  for (i in 1:length(ls)) {
    r <- getADNormarlityTest(eval(parse(text = ls[i]), envir = .GlobalEnv)$residuals, warning)
    e[i,1] <- r$statistic
    e[i,2] <- r$p.value
  }
  rownames(e) <- ls
  colnames(e) <- c("A", "p")
  return(e)
}
linearModel.getBestFitModelNameByRSquared.querySelector <- function (char = "", adj.r.squared = FALSE, ...) {
  ls <- querySelector(char, ...)
  if (adj.r.squared) e <- summary.lm(eval(parse(text = ls[1]), envir = .GlobalEnv))$adj.r.squared
  else e <- summary.lm(eval(parse(text = ls[1]), envir = .GlobalEnv))$r.squared
  c <- 1
  for (i in 2:length(ls)) {
    if (adj.r.squared) x <- summary.lm(eval(parse(text = ls[i]), envir = .GlobalEnv))$adj.r.squared
    else x <- summary.lm(eval(parse(text = ls[i]), envir = .GlobalEnv))$r.squared
    if (x > e) {
      e <- x
      c <- i
    }
  }
  return(ls[c])
}
linearModel.getLargestedModelName.querySelector <- function (char = "", ...) {
  ls <- querySelector(char, ...)
  su <- summary.lm(eval(parse(text = ls[1]), envir = .GlobalEnv))
  df <- su$df[2]
  c <- 1
  for (i in 2:length(ls)) {
    su <- summary.lm(eval(parse(text = ls[i]), envir = .GlobalEnv))
    x <- su$df[2]
    r <- su$adj.r.squared
    if (x < df) {
      df <- x
      c <- i
    } else if ((x = df)&&(su$adj.r.squared > r)) {
      c <- i
    }
  }
  return(ls[c])
}
linearModel.getBPHomoscedasticityTestList.querySelector <- function (char = "", warning = TRUE, ...) {
  ls <- querySelector(char, ...)
  e <- list()
  for (i in 1:length(ls)) {
    e[i] <- list(ncvTest(eval(parse(text = ls[i]), envir = .GlobalEnv)))
    p <- e[[i]]$p
    if (p < 0.05&&warning) warning(ls[i], " is heteroscedastic, p = ", sprintf("%.3f", p))
  }
  names(e) <- ls
  return(e)
}
linearModel.getANOVAList.querySelector <- function (char = "", c = "bestFitByRSquared", ...) {
  # Perform ANOVA each with the best fit
  ls <- querySelector(char, ...)
  switch(
    c,
    "shortest" = c <- ls[1],
    "largest" = c <- linearModel.getLargestedModelName.querySelector(char, ...),
    "bestFitByRSquared" = c <- linearModel.getBestFitModelNameByRSquared.querySelector(char, ...)
  )
  ls <- ls[is.na(match(ls, c))]
  e <- list()
  for (i in 1:length(ls)) {
    e[i] <- list(eval(parse(text = paste("anova(", c, ",", ls[i], ")")), envir = .GlobalEnv))
  }
  names(e) <- ls
  return(e)
}
linearModel.getANOVA.querySelector <- function (char = "", c = "bestFitByRSquared", ...) {
  # Perform ANOVA with all models together
  ls <- querySelector(char, ...)
  switch(
    c,
    "shortest" = c <- ls[1],
    "largest" = c <- linearModel.getLargestedModelName.querySelector(char, ...),
    "bestFitByRSquared" = c <- linearModel.getBestFitModelNameByRSquared.querySelector(char, ...)
  )
  ls <- ls[is.na(match(ls, c))]
  return(eval(parse(text = paste("anova(", c, ",", paste(ls, collapse = ","), ")")), envir = .GlobalEnv))
}
linearModel.getVIFList.querySelector <- function (char = "", ...) {
  ls <- querySelector(char, ...)
  e <- list()
  for (i in 1:length(ls)) {
    try(e[i] <- list(vif(eval(parse(text = ls[i]), envir = .GlobalEnv))))
  }
  names(e) <- ls
  return(e)
}
linearModel.getOutlierTestList.querySelector <- function (char = "", ..., warning = TRUE) {
  ls <- querySelector(char, ...)
  e <- list()
  for (i in 1:length(ls)) {
    try({t <- outlierTest(eval(parse(text = ls[i]), envir = .GlobalEnv))
    e[i] <- list(t)
    if (t$bonf.p < 0.05&&warning) warning(ls[i], " has outlier(s)")})
  }
  names(e) <- ls
  return(e)
}
linearModel.getSignificantAutocorrelationBGTestList.querySelector <- function (char = "", ..., order = 9, warning = TRUE) {
  ls <- querySelector(char, ...)
  e <- list()
  ln <- c()
  for (i in 1:length(ls)) {
    for (n in 1:order) tryCatch({
      lm <- eval(parse(text = ls[i]), envir = .GlobalEnv)
      bg <- bgtest(lm, n)
      p <- bg$p.value
      if (p < 0.05) {
        e[length(e) + 1] <- list(bg)
        ln <- c(ln, ls[i])
        if (warning) warning(ls[i], " is autocorrelated, p = ", sprintf("%.3f", p))
      }
    }, error = function(error) {})}
  names(e) <- ln
  return(e)
}
linearModel.getModelSelectionTable.querySelector <- function (char = "", warning = FALSE, rank = "AIC", ...) {
  suppressWarnings({ls <- querySelector(char, ...)
  a <- eval(parse(text = paste("AIC(", paste(ls, collapse = ","), ")")), envir = .GlobalEnv)
  b <- eval(parse(text = paste("BIC(", paste(ls, collapse = ","), ")")), envir = .GlobalEnv)
  a[3] <- b[2]
  if (is.null(rank)) return(a)
  if (rank == "AIC") a <- a[order(a[2]), ]
  if (rank == "BIC") a <- a[order(a[3]), ]
  })
  return(a)
}
linearModel.getLOOCVList.querySelector <- function (char = "", ...) {
  ls <- querySelector(char, ...)
  e <- list()
  for (i in 1:length(ls)) {
    c <- eval(parse(text = paste(ls[i], "$call")), envir = .GlobalEnv)
    f <- c$formula
    e[i] <- list(eval(parse(text = paste("train(", f[2], "~", paste(f[3:length(f)], collapse = ""), ", data =", c$data, ", method = 'lm', trControl = trainControl(method = 'LOOCV'), na.action = na.omit)$results")), envir = .GlobalEnv))
  }
  names(e) <- ls
  return(e)
}
linearModel.getCombinedSelectionTable <- function (table, ls, rm = FALSE, rm.intercept = TRUE) {
  n <- names(ls)
  c <- ncol(table)
  tb <- deparse(substitute(table))
  l <- deparse(substitute(ls))
  for (i in 1:length(n)) {
    su <- summary(eval(parse(text = n[i]), envir = .GlobalEnv))
    t <- cbind(ls[[i]],su$r.squared, su$adj.r.squared)
    table[rownames(table) == n[i],(c + 1):(c + length(t))] <- t
    table[rownames(table) == n[i],1] <- su$df[2]
  }
  if (rm.intercept) table <- table[ , ! names(table) %in% "intercept"]
  if (rm) eval(parse(text = paste("rm(", tb, ",", l, ", pos = '.GlobalEnv')")), envir = .GlobalEnv)
  return(table)
}
getResults.querySelector <- function (char = "", ...) {
  ls <- querySelector(char, ...)
  e <- list()
  for (i in 1:length(ls)) {
    e[i] <- list(eval(parse(text = ls[i]), envir = .GlobalEnv))
  }
  names(e) <- ls
  return(e)
}
dataFrame <- data.frame(
  year = c(2000:2020),
  temp = c(0, 0.1, -0.3, -0.2, -0.3, -0.4, 0.1, -0.1, 0.6, -0.5, 0.4, -0.2, 0, 0.1, 0, 0.2, -0.6, 0.8, 0.4, 0.5, NA),
  species_1 = as.integer(c(57, 61, 68, 90, 110, 147, 225, 208, 205, 172, 250, 199, 267, 290, 344, 410, 389, 545, 253, 268, 239)),
  species_2 = as.integer(c(38, 51, 67, 62, 84, 84, 112, 123, 104, 122, 118, 117, 138, 148, 137, 179, 191, 194, 266, 351, 314))
)
#dataFrame <- read.csv("pop_data.csv") # Uncomment this line to load data frame from "./pop_data.csv"

# §2. Data visualisation --------------------------------------------------
S1S2vsYr.plot <- ggplot(dataFrame) +
  scale_x_continuous(limits = c(1999.6, 2020.4), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 930), expand = c(0, 0), sec.axis = sec_axis(~ (. - 400) / 400, name = "Temperature anomaly (°C)")) +
  geom_point(aes(x = year, y = temp * 400 + 400), colour = "#f9bb00", size = 1.5, na.rm = TRUE) +
  geom_line(aes(x = year, y = temp * 400 + 400), colour = "#f9bb00", size = 0.5, na.rm = TRUE) +
  geom_smooth(aes(x = year, y = species_1), colour = "#00a3ff", method = "glm", formula = y~x, method.args = list(family = Gamma)) +
  geom_smooth(aes(x = year, y = species_2), colour = "#5dd934", method = "glm", formula = y~x, method.args = list(family = Gamma)) +
  geom_point(aes(x = year, y = species_1), colour = "#00a3ff", size = 1.5) +
  geom_point(aes(x = year, y = species_2), colour = "#5dd934", size = 1.5) +
  geom_segment(aes(x = 2013.2, y = 850, xend = 2015.2, yend = 850),
               size = 1, colour = "#00a3ff") +
  geom_point(aes(x = 2014.2, y = 850), size = 1.5, colour = "#00a3ff") +
  geom_text(aes(x = 2017, y = 850, label = "Species 1"),
            stat = "unique", size = 4, colour = "#00a3ff") +
  geom_segment(aes(x = 2013.2, y = 800, xend = 2015.2, yend = 800),
               size = 1, colour = "#5dd934") +
  geom_point(aes(x = 2014.2, y = 800), size = 1.5, colour = "#5dd934") +
  geom_text(aes(x = 2017, y = 800, label = "Species 2"),
            stat = "unique", size = 4, colour = "#5dd934") +
  geom_segment(aes(x = 2013.2, y = 900, xend = 2015.2, yend = 900),
               size = 1, colour = "#f9bb00") +
  geom_point(aes(x = 2014.2, y = 900), size = 1.5, colour = "#f9bb00") +
  geom_text(aes(x = 2017.4, y = 900, label = "Temperature"),
            stat = "unique", size = 4, colour = "#f9bb00") +
  geom_text(aes(x = 2015.7, y = 545, label = "2017"),
            stat = "unique", color = "#00a3ff") +
  ylab("Population Size") +
  xlab("Time (Year)") +
  theme_minimal()

# §3. Calculate intrinsic rate of increase using exponential model --------
exponentialModel.getDiscreteGrowthRate <- function (N) {
  r <- c()
  for (i in 1:length(N)-1){
    r[i] <- log(N[i+1]/N[i])
  }
  r[length(N)] <- NA
  return(r)
}
dataFrame.fit <- dataFrame
dataFrame.fit$rSpecies_1 <- exponentialModel.getDiscreteGrowthRate(dataFrame$species_1)
dataFrame.fit$rSpecies_2 <- exponentialModel.getDiscreteGrowthRate(dataFrame$species_2)
dataFrame.fit <- head(dataFrame.fit, -1)
## §3.1. Check for normality with Anderson-Darling test -------------------
normarlityAD.table <- dataFrame.getADNormarlityTestTable(dataFrame.fit, "temp", "rSpecies_1", "rSpecies_2") 
# In exponential growth model, population size is unlikely to be normally distributed

# §4. Linear temperature dependence of growth rate ------------------------
## §4.1. Regression visualisation -----------------------------------------
rS1rS2vsT.plot <- ggplot(dataFrame.fit) + # Remove the last row containing NA
  geom_hline(yintercept = 0, colour = "gray") +
  geom_vline(xintercept = 0, colour = "gray") +
  scale_x_continuous(limits = c(-0.63, 0.83), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-0.8, 0.6), expand = c(0, 0)) +
  geom_smooth(aes(x = temp, y = rSpecies_1), colour = "#00a3ff", method = "lm", formula = y~x) +
  geom_smooth(aes(x = temp, y = rSpecies_2), colour = "#5dd934", method = "lm", formula = y~x) +
  geom_point(aes(x = temp, y = rSpecies_1), colour = "#00a3ff", size = 1.5) +
  geom_point(aes(x = temp, y = rSpecies_2), colour = "#5dd934", size = 1.5) +
  geom_segment(aes(x = 0.45, y = 0.55, xend = 0.57, yend = 0.55),
               size = 1, colour = "#00a3ff") +
  geom_point(aes(x = 0.51, y = 0.55), colour = "#00a3ff") +
  geom_text(aes(x = 0.68, y = 0.55, label = "Species 1"),
            stat = "unique", color = "#00a3ff") +
  geom_text(aes(x = 0.72, y = -0.767, label = "2017"),
            stat = "unique", color = "#00a3ff") +
  geom_segment(aes(x = 0.45, y = 0.47, xend = 0.57, yend = 0.47),
               size = 1, colour = "#5dd934") +
  geom_point(aes(x = 0.51, y = 0.47), colour = "#5dd934") +
  geom_text(aes(x = 0.68, y = 0.47, label = "Species 2"),
            stat = "unique", color = "#5dd934") +
  geom_boxplot(aes(x = temp, y = -0.5), notch = TRUE, size = 0.7, width = 0.2, colour = "#f9bb00", fill = "#f9bb0070") +
  ylab(bquote("Intrinsic rate of increase " (year^-1))) +
  xlab("Temperature anomaly (°C)") +
  theme_minimal()
## §4.2. Calculate regression parameters ----------------------------------
rS1vsT.lm <- lm(rSpecies_1 ~ temp, data = dataFrame.fit)
rS1vsT.su <- summary.lm(rS1vsT.lm)
rS2vsT.lm <- lm(rSpecies_2 ~ temp, data = dataFrame.fit)
rS2vsT.su <- summary.lm(rS2vsT.lm)
## §4.3. Climate change ---------------------------------------------------
temp.t <- t.test(dataFrame$temp)

# §5. Assuming constant intrinsic growth rates ----------------------------
rS1.mean.ci <- mean.getConfidenceInterval(dataFrame.fit$rSpecies_1)
rS2.mean.ci <- mean.getConfidenceInterval(dataFrame.fit$rSpecies_2)
dataFrame$tS1 <- dataFrame$year - 1997
dataFrame$tS2 <- dataFrame$year - 1999
S1.lm <- lm(log(species_1) ~ tS1, data = dataFrame)
S1.su <- summary(S1.lm)
S2.lm <- lm(log(species_2) ~ tS2, data = dataFrame)
S2.su <- summary(S2.lm)
linearModel.compareCoefficientByWelchTTest <- function (lm, y) {
  s <- summary(lm)
  xn <- s$df[1] + s$df[2] - 2
  xv <- s$coefficients[2,2]
  yn <- length(y) - 1
  yv <- var(y)^2/(yn + 1)
  t <- (s$coefficients[2,1] - mean(y))/sqrt(xv + yv)
  df <- (xv + yv)^2/(xv^2/xn + yv^2/yn)
  p <- 2*pt(abs(t), df, lower.tail = FALSE)
  return(paste("t = ", sprintf("%.3f", t), ", df = ", sprintf("%.3f", df), ", p = ", sprintf("%.3f", p), sep = ""))
}
rS1.t <- linearModel.compareCoefficientByWelchTTest(S1.lm, dataFrame.fit$rSpecies_1)
rS2.t <- linearModel.compareCoefficientByWelchTTest(S2.lm, dataFrame.fit$rSpecies_2)

# §6. Competitive Lotka-Volterra model ------------------------------------
# Expect linear dependence of rate of one species on the population size 
# of the other under competitive Lotka-Volterra model
rS1vsS1S2T.lm <- lm(rSpecies_1 ~ species_1 + species_2 + temp, data = dataFrame.fit)
rS1vsS1S2T.su <- summary.lm(rS1vsS1S2T.lm)
rS2vsS1S2T.lm <- lm(rSpecies_2 ~ species_1 + species_2 + temp, data = dataFrame.fit)
rS2vsS1S2T.su <- summary.lm(rS2vsS1S2T.lm)

# §7. Introduce time lags -------------------------------------------------
# Result: up to 9 lags 9 (df = 1) can be significant, thus create models
# with number of lags from 0 to the maximum (9), test for significance
dataFrame.introduceLags <- function (df, x, ..., max = TRUE, pattern = ".L") {
  x <- c(x, ...)
  if (max)  m <- floor((nrow(df) - 3) / 2)
  else  m <- floor((nrow(df) - length(x) - 1) / (length(x) + 1))
  # Each regressors cost 1 df, each addition of lag per regressors lose 1 
  # 1 df, each addition of lags in all regressors lose 1 df because 
  # discontinuous data were removed, (Intercept) cost 1 df. df should > 0
  r <- nrow(df)
  for (v in x) {
    for (i in 1:m) {
      n <- ncol(df) + 1
      df[,n] <- df[,colnames(df) == v][c((r + 1):(r + i), 1:(r - i))]
      colnames(df)[n] <- paste(v, i, sep = pattern)
    }
  }
  return(df)
}
linearModel.getSignificantLagsList <- function (df, y, x1, x2, x3, df.min = 1, m = "max", gls = FALSE, vcovHC = FALSE, ncv.action = NULL) {
  if (gls) type <- "gls" else type <- "lm"
  x <- c(x1, x2, x3)
  z <- list()
  df.name <- deparse(substitute(df))
  e <- c()
  k <- 0
  for (i in 1:length(x)) {
    d <- ls (df, pattern = x[i])
    f <- c(1)
    for (g in 1:length(d)) {
      f[g+1] <- d[g]
    }
    z[[i]] <- f
  }
  dof <- nrow(df) - df.min + 1
  if (m == "max")  m <- floor((nrow(df) - 3) / 2) + 2
  else if (m == "min") m <- floor((nrow(df) - length(x) - 1) / (length(x) + 1)) + 2
  h <- 0
  for (a in 1:m) {
    for (b in 1:m) {
      for (c in 1:m) {
        if (a + b + c + max(a, b, c) - 4 < dof) {
          # Each regressor cost 1 df, each of a, b, c = number of lags + 2
          # = number of regressor, max(a, b, c) -1 return number of
          # observations removed + 1, (Intercept) cost 1 df. df should > 0
          k <- k + 1
          t <- paste(type, "(", y, " ~ ",
                     paste(paste(z[[1]][1:a], collapse = " + "),
                           paste(z[[2]][1:b], collapse = " + "),
                           paste(z[[3]][1:c], collapse = " + "), sep = " + "),
                     ", data = ", df.name, ", na.action = na.omit)", sep = "")
          lm = eval(parse(text = t), envir = .GlobalEnv)
          ho <- TRUE
          if (!is.null(ncv.action)) {
            n <- ncvTest(lm)$p
            if (n < 0.05) {
              if (ncv.action[1] == "gls") {
                t <- paste("gls(", y, " ~ ",
                           paste(paste(z[[1]][1:a], collapse = " + "),
                                 paste(z[[2]][1:b], collapse = " + "),
                                 paste(z[[3]][1:c], collapse = " + "), sep = " + "),
                           ", data = ", df.name, ", na.action = na.omit)", sep = "")
                lm = eval(parse(text = t), envir = .GlobalEnv)
                d <- summary(lm)$tTable
              } else if (ncv.action[1] == "vcovHC") {
                d <- coeftest(lm, vcov = vcovHC(lm, type = ncv.action[2]))
              } else message(t, " is heteroscedastic, p = ", n)
              h <- h + 1
              ho <- FALSE
            }
          }
          if (ho) {
            if (gls) d <- summary(lm)$tTable
            else if (vcovHC == FALSE) d <- summary.lm(lm)$coefficients
            else d <- coeftest(lm, vcov = vcovHC(lm, type = vcovHC))
          }
          r <- d[,4]
          if (a > 1) pa <- r[rownames(d) == z[[1]][a]] else pa <- 1
          if (b > 1) pb <- r[rownames(d) == z[[2]][b]] else pb <- 1
          if (c > 1) pc <- r[rownames(d) == z[[3]][c]] else pc <- 1
          if ((pa < 0.05)&&(pb < 0.05)&&(pc < 0.05)) e[length(e) + 1] <- t
        }
      }
    }
  }
  if (!is.null(ncv.action)) e[length(e) + 1] <- paste("A total of", k, "models are generated, among which", h, 'are(is) heteroscedastic')
  else e[length(e) + 1] <- paste("A total of", k, "models are generated")
  return(e)
}
dataFrame.fit <- dataFrame.introduceLags(dataFrame.fit, "temp", "species_1", "species_2")
linearModel.getSignificantLagsList(dataFrame.fit, "rSpecies_1" , "temp", "species_1", "species_2", df.min = 2, ncv.action = c("vcovHC", "HC4m"))
# A total of 359 models are generated with df > 1, "HC4m" estimator for robust 
# covariance matrix is applied to heteroscedastic models detected by 
# ncvTest() as suggested by Cribari-Neto and Da Silva, 2011
# [1]: the same as rS1vsS1S2T.lm
summary.lm(lm(rSpecies_1 ~ temp + species_1 + species_2, data = dataFrame.fit))
# [2]: t-test shows that regressor species_2:species_2.L5 is insignificant, thus removed
summary.lm(lm(rSpecies_1 ~ temp + species_1 + species_2 + species_2.L1 + species_2.L2 + species_2.L3 + species_2.L4 + species_2.L5 + species_2.L6, data = dataFrame.fit))
rS1Lag1.lm <- lm(rSpecies_1 ~ temp + species_1 + species_2.L6, data = dataFrame.fit)
rS1Lag1.su <- summary.lm(rS1Lag1.lm)
# [3]: t-test shows that regressor species_2 without lag is insignificant, thus removed
summary.lm(lm(rSpecies_1 ~ temp + temp.L1 + temp.L2 + species_1 + species_1.L1 + species_1.L2 + species_1.L3 + species_1.L4 + species_2 + species_2.L1 + species_2.L2, data = dataFrame.fit))
rS1Lag2.lm <- lm(rSpecies_1 ~ temp + temp.L1 + temp.L2 + species_1 + species_1.L1 + species_1.L2 + species_1.L3 + species_1.L4 + species_2.L1 + species_2.L2, data = dataFrame.fit)
rS1Lag2.su <- summary.lm(rS1Lag2.lm)
# For rSpecies_2:
linearModel.getSignificantLagsList(dataFrame.fit, "rSpecies_2" , "temp", "species_1", "species_2", df.min = 2, ncv.action = c("vcovHC", "HC4m"))
# No model selected

# §8. Power analysis and quality control ----------------------------------
normarlityAD.table <- rbind(normarlityAD.table, linearModel.getADNormarlityTestTable.querySelector(".lm")) # Check for residual normality
Lag.vif.ls <- linearModel.getVIFList.querySelector(".lm", inl = "Lag")
outliers.ls <- linearModel.getOutlierTestList.querySelector(".lm", warning = FALSE)
autocorrelationBG.ls <- linearModel.getSignificantAutocorrelationBGTestList.querySelector(".lm", order = 9, warning = FALSE) # Up to 9 lags otherwise df = 0
## §8.1. Akaike and Bayesian information criterion (AIC, BIC) -------------
fitnessIC.table <- linearModel.getModelSelectionTable.querySelector("lm", rank = "AIC")
# Stepwise procedure: removing any regressor fail to improve AIC of the following models:
step(rS1Lag1.lm)
step(rS1Lag2.lm)
step(rS1vsS1S2T.lm)
step(rS2vsS1S2T.lm) # (removing Species 1)
rS2vsS2T.lm <- lm(rSpecies_2 ~ species_2 + temp, data = dataFrame.fit)
rS2vsS2T.su <- summary(rS2vsS2T.lm)
## §8.2. Get the best fit by comparing Adjusted R-squared values ----------
# Cohen 1988: R2 < 0.02 - Very weak; 0.02 <= R2 < 0.13 - Weak; 0.13 <= R2 < 0.26 - Moderate; R2 >= 0.26 - Substantial
rS1.bestFitByAdjR.name <- linearModel.getBestFitModelNameByRSquared.querySelector("rS1vs", inl = ".lm", adj.r.squared = TRUE)
rS1.bestFitByAdjR.lm <- eval(parse(text = rS1.bestFitByAdjR.name))
avPlots(rS1.bestFitByAdjR.lm, col.lines = "#00a3ff", layout = c(1,3), main = "Species 1 Rate of Increase Added Variable Plots") # Added variable plots
rS2.bestFitByAdjR.name <- linearModel.getBestFitModelNameByRSquared.querySelector("rS2vs", inl = ".lm", adj.r.squared = TRUE)
rS2.bestFitByAdjR.lm <- eval(parse(text = rS2.bestFitByAdjR.name))
avPlots(rS2.bestFitByAdjR.lm, col.lines = "#5dd934", layout = c(1,2), main = "Species 2 Rate of Increase Added Variable Plots")
rS1.bestFitByAdjR.su <- summary.lm(rS1.bestFitByAdjR.lm)
rS2.bestFitByAdjR.su <- summary.lm(rS2.bestFitByAdjR.lm)
## §8.3. Analysis of variance (ANOVA) -------------------------------------
# For every model that does not involve lags, compare SSR (sum of squared
# residuals) with the largest (with all regressors) model to test if
# each addition of regressors is significant in reducing SSR
#rS1.anova.ls <- linearModel.getANOVAList.querySelector("rS1vs", inl = ".lm") # Do ANOVA Test for every linear model of rS1 without lags
#rS2.anova.ls <- linearModel.getANOVAList.querySelector("rS2vs", inl = ".lm")
## §8.4. Heteroskedasticity-consistent standard errors --------------------
# Check homoscedasticity with Breusch–Pagan test
homoscedasticityBP.ls <- linearModel.getBPHomoscedasticityTestList.querySelector(".lm", warning = FALSE)
# Since many models of rS1 are heteroscedastic (homoscedasticityBP.ls),
# heteroscedasticity robust t-tests were also performed
rS1vsT.t <- coeftest(rS1vsT.lm, vcov = vcovHC(rS1vsT.lm, "HC4m"))
## §8.5. Calculate Lotka-Volterra model parameters ------------------------
# rS1 = r1(K1-S1-αS2)/K1 + γT + Ɛ(error) <=> rS1 = (-r1/K1)*S1 + (-r1α/K1)*S2 + r1 + γT + Ɛ
linearModel.getLotkaVolterraModelParameters <- function (lm, targetSpecies, interactingSpecies, errorPropagation = TRUE) {
  su <- summary.lm(lm)$coefficients
  r <- su[,1][1] # When the error term Ɛ is ignored, (Intercept) = r
  l1 <- grep(targetSpecies, rownames(su))
  s1 <- sum(su[,1][l1])
  l2 <- grep(interactingSpecies, rownames(su))
  s2 <- sum(su[,1][l2])
  K <- -r / s1
  α <- s2 / s1
  if (errorPropagation) {
    # Linear combinations: f = Σai*xi <=> σf^2 = Σ(ai^2 + σi^2)
    if (length(l1) > 1) {
      s1.se <- 0
      for (i in l1) {
        s1.se <- s1.se + (su[,1][i]*su[,2][i])^2
      }
      s1.se <- sqrt(s1.se)
    } else s1.se <- su[,2][l1]
    if (length(l2) > 1) {
      s2.se <- 0
      for (t in l2) {
        s2.se <- s2.se + (su[,1][t]*su[,2][t])^2
      }
      s2.se <- sqrt(s2.se)
    } else s2.se <- su[,2][l2]
    # Non-linear combinations: SD is estimated from a simplified Taylor
    # expansion. Since their sample sizes are the same, this method also
    # works for SE, given that SE = SD / sqrt(n)
    r.se <- su[,2][1]
    K.se <- abs(K * sqrt((r.se / r)^2 + (s1.se / s1)^2))
    α.se <- abs(α * sqrt((s2.se / s2)^2 + (s1.se / s1)^2))
    return(paste("r = ", sprintf("%.3f", r), " ± ", sprintf("%.3f", r.se), ", K = ", sprintf("%.3f", K), " ± ", sprintf("%.3f", K.se), ", α = ", sprintf("%.3f", α), " ± ", sprintf("%.3f", α.se), sep = ""))
  } else return(paste("r = ", r, ", K = ", K, ", α = ", α, sep = ""))
}
rS1vsS1S2T.lv <- linearModel.getLotkaVolterraModelParameters(rS1vsS1S2T.lm, "species_1", "species_2", errorPropagation = TRUE)
rS1Lag1.lv <- linearModel.getLotkaVolterraModelParameters(rS1Lag1.lm, "species_1", "species_2", errorPropagation = TRUE)
rS1Lag2.lv <- linearModel.getLotkaVolterraModelParameters(rS1Lag2.lm, "species_1", "species_2", errorPropagation = TRUE)
rS2vsS1S2T.lv <- linearModel.getLotkaVolterraModelParameters(rS2vsS1S2T.lm, "species_2", "species_1", errorPropagation = TRUE)
# §8.6. Leave One Out Cross-Validation (LOOCV) ----------------------------
# (This might take some time)
fitnessLOOCV.ls <- linearModel.getLOOCVList.querySelector(".lm")
fitness.table <- linearModel.getCombinedSelectionTable(fitnessIC.table, fitnessLOOCV.ls, rm = TRUE)
# §8.7.  Squared Semi-partial Correlation ---------------------------------
rS1vsS1S2T.pcor <- linearModel.getDeltaRsquare(rS1vsS1S2T.lm)
rS2vsS1S2T.pcor <- linearModel.getDeltaRsquare(rS2vsS1S2T.lm)
rS1Lag1.pcor <- linearModel.getDeltaRsquare(rS1Lag1.lm)

# §9. Print results ------------------------------------------------------
S1S2vsYr.plot
normarlityAD.table
getResults.querySelector(".ls")
getResults.querySelector("temp.")
rS1.t
rS1.mean.ci
S1.su
rS2.t
rS2.mean.ci
S2.su
getResults.querySelector("vsT", exl = ".lm")
getResults.querySelector("vsS1S2T", exl = ".lm")
rS2vsS2T.su
getResults.querySelector(".lv")
getResults.querySelector("rS1Lag1", exl = ".lm")
fitness.table