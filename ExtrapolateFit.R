# Idea is to extrapolate the fit to an very high age.
# Best way would be to go into the fitting function but that is risky.
# Therefore we take the finished analysis and add a few years.
library(tidyverse)

extrapolateAnalysisFrom109 <- function (analysis) {
  par <- analysis$resultList$gradResults$par
  ages <- seq(110, 125)

  if (analysis$resultList$gradFunc$name == "GompertzGPD" || analysis$resultList$gradFunc$name == "MakehamGPD") {
    rates <- analysis$resultList$gradFunc$hazardFuncGPD(par, ages)
    mod <- log(rates)
    extrapolatedData <- data.frame(Age = ages, mod = mod, model = analysis$resultList$gradFunc$name)
    extrapolatedData <- bind_rows(analysis$plottingData, extrapolatedData)
  } else if (grepl("Hermite", analysis$resultList$gradFunc$name, fixed=TRUE)) {
    # type, X0, X1
    modelSpecifications <- analysis$resultList$modelSpecifications
    rates <- hermite(modelSpecifications$X0, modelSpecifications$X1, par, ages)
    mod <- log(rates)
    # Overwrite max value
    maxMod <- -10000
    for (i in seq_along(mod)) {
      if (mod[i] > maxMod) {
        maxMod <- mod[i]
      } else {
        mod[i] <- maxMod
      }
    }
    extrapolatedData <- data.frame(Age = ages, mod = mod, model = analysis$resultList$gradFunc$name)
    extrapolatedData <- bind_rows(analysis$plottingData, extrapolatedData)
  } else {
    rates <- analysis$resultList$gradFunc$func(par, ages)
    mod <- log(rates)
    extrapolatedData <- data.frame(Age = ages, mod = mod, model = analysis$resultList$gradFunc$name)
    extrapolatedData <- bind_rows(analysis$plottingData, extrapolatedData)
  }

  analysis$extrapolatedData <- extrapolatedData
  return(analysis)
}

hermite <- function (x0, x1, para, ages) {
  para <- c(if_else(!is.na(para["hermiteAlpha"]),para["hermiteAlpha"],0),
            if_else(!is.na(para["hermiteOmega"]),para["hermiteOmega"],0),
            if_else(!is.na(para["hermiteM0"]),para["hermiteM0"],0),
            if_else(!is.na(para["hermiteM1"]),para["hermiteM1"],0),
            if_else(!is.na(para["hermiteV"]),para["hermiteV"],0))

  t <- (ages - x0) / (x1 - x0)

  h00x <- function(t) {
    return((1.0 + 2 * t) * (1.0 - t) * (1.0 - t))
  }
  h10x <- function(t) {
    return(t * (1.0 - t) * (1.0 - t))
  }
  h01x <- function(t) {
    return(t * t * (3.0 - 2 * t))
  }
  h11x <- function(t) {
    return(t * t * (t - 1.0))
  }
  hQuarticx <- function(t) {
    return(16 * t ^ 2 * (1 - t) ^ 2)
  }

  hermiteMu <- function(para, t) {
    return(exp(
      para[1] * h00x(t) + para[2] * h01x(t) +
        para[3] * h10x(t) + para[4] * h11x(t) +
        para[5] * hQuarticx(t)
    ))
  }

  rates <- hermiteMu(para, t)
}