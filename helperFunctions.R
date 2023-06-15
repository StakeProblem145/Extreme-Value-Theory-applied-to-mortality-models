# Some helper functions
findFullCohortYears <- function(dataFrame) {
  start <- 0
  end <- 0
  startl <- 0
  endl <- 0
  for (c in sort(unique(dataFrame$Cohort))) {
    current <-
      length(filter(dataFrame, Cohort == c & Gender == "Female")$Age)
    if (current > startl) {
      start <- c
      startl <- current
    }
    if (current >= endl) {
      end <- c
      endl <- current
    }
  }
  return(c(start, end))
}



bestGompertzGPDThresholdAge <- function(dataFrame) {
  returnDf <-
    setNames(data.frame(matrix(ncol = 2, nrow = 0)), c("tA", "AIC"))
  for (tA in c(85:103)) {
    AIC <-
      graduateGompertzGPD(dataFrame, thresholdAge = tA)$modelResults$AIC
    returnDf[nrow(returnDf) + 1, ] <- c(tA, AIC)
  }
  return(returnDf)
}

bestMakehamGPDThresholdAge <- function(dataFrame) {
  returnDf <-
    setNames(data.frame(matrix(ncol = 2, nrow = 0)), c("tA", "AIC"))
  for (tA in c(85:103)) {
    AIC <-
      graduateMakehamGPD(dataFrame, thresholdAge = tA)$modelResults$AIC
    returnDf[nrow(returnDf) + 1, ] <- c(tA, AIC)
  }
  return(returnDf)
}

bestHermiteII <- function(dataFrame, range = c(105:125)) {
  returnDf <-
    setNames(data.frame(matrix(ncol = 2, nrow = 0)), c("X1", "AIC"))
  for (X1 in range) {
    temp <- graduateHermite(dataFrame, "II", X1 = X1)
    AIC <- temp$modelResults$AIC
    returnDf[nrow(returnDf) + 1, ] <- c(X1, AIC)
  }
  return(returnDf)
}

bestHermiteV <- function(dataFrame, range = c(105:125)) {
  returnDf <-
    setNames(data.frame(matrix(ncol = 2, nrow = 0)), c("X1", "AIC"))
  for (X1 in range) {
    temp <- graduateHermite(dataFrame, "V", X1 = X1)
    AIC <- temp$modelResults$AIC
    returnDf[nrow(returnDf) + 1,] <- c(X1, AIC)
  }
  return(returnDf)
}

addFittedColumnTo109 <- function(analysisData) {
  temp <- analysisData$plottingData %>%
    mutate(fitted = if_else(Age <= 104, "Fit", "Extrapolation"))
  A <- filter(temp, fitted == "Fit")
  dup <- slice_tail(A, n = 1)
  B <- filter(temp, fitted == "Extrapolation")
  B <- bind_rows(dup, B)
  B$fitted <- "Extrapolation"
  analysisData$plottingData <- bind_rows(A, B)
  
  temp <- analysisData$residualData %>%
    mutate(fitted = if_else(Age <= 104, "Fit", "Extrapolation"))
  A <- filter(temp, fitted == "Fit")
  dup <- slice_tail(A, n = 1)
  B <- filter(temp, fitted == "Extrapolation")
  B <- bind_rows(dup, B)
  B$fitted <- "Extrapolation"
  analysisData$residualData <- bind_rows(A, B)

  if("extrapolatedData" %in% names(analysisData)){
    temp <- analysisData$extrapolatedData %>%
      mutate(fitted = if_else(Age <= 104, "Fit", "Extrapolation"))
    A <- filter(temp, fitted == "Fit")
    dup <- slice_tail(A, n = 1)
    B <- filter(temp, fitted == "Extrapolation")
    B <- bind_rows(dup, B)
    B$fitted <- "Extrapolation"
    analysisData$extrapolatedData <- bind_rows(A, B)
  }
  
  return(analysisData)
}