require(tidyverse)

##### CONSTANTS #####
# DO NOT CHANGE THESE CONSTANTS. USE THEM AS A TEMPLATE AND THEN AS AN ARGUMENT
INIT_PARAMETER_CLA <-
  list(
    initAlpha = -10,
    initBeta = 0.1,
    initEpsilon = -10,
    initRho = 0
  )
INIT_PARAMETER_HMT <-
  list(
    hermiteAlpha = -10,
    hermiteOmega = -1,
    hermiteM0 = 0,
    hermiteM1 = 0,
    hermiteV = 0
  )
INIT_PARAMETER_GPD <-
  list(
    initAlpha = -10,
    initBeta = 0.1,
    initEpsilon = -10,
    initRho = 0,
    initXi = 0.1
  )

OPTIM_CONTROL_PARAMETERS <-
  list(
    optimMethod = "BFGS",
    gradMethod = "Richardson",
    # Maps to factr instead of reltol because of L-BFGS-B
    relTolerance = 1e-20,
    maxIteration = 1000,
    minimumMu = 1e-20
  )

OPTIM_CONSTRAINTS_CLA <- list(
  alphaLower = -Inf,
  betaLower = -Inf,
  epsilonLower = -Inf,
  rhoLower = -Inf,
  alphaUpper = Inf,
  betaUpper = Inf,
  epsilonUpper = Inf,
  rhoUpper = Inf
)

OPTIM_CONSTRAINTS_HMT <- list(
  hermiteAlphaLower = -Inf,
  hermiteOmegaLower = -Inf,
  hermiteM0Lower = -Inf,
  hermiteM1Lower = -Inf,
  hermiteV = -Inf,
  hermiteAlphaUpper = Inf,
  hermiteOmegaUpper = Inf,
  hermiteM0Upper = Inf,
  hermiteM1Upper = Inf,
  hermiteV = Inf
)


### Obviously dont change these
INTEGRATED_HAZARD_FUNCTIONS <- list(
  gompertzFunc = function(para, ages) {
    ((exp(para["beta"]) - 1) / para["beta"]) * exp(para["alpha"] + para["beta"] * ages)
  },
  makehamFunc = function(para, ages) {
    exp(para["epsilon"]) + ((exp(para["beta"]) - 1) / para["beta"]) * exp(para["alpha"] + para["beta"] * ages)
  },
  perksFunc = function(para, ages) {
    (1 / para["beta"]) * log((1 + exp(para["alpha"] + para["beta"] * (ages + 1)))
                             / (1 + exp(para["alpha"] + para["beta"] * ages)))
  },
  beardFunc = function(para, ages) {
    (exp(-para["rho"]) / para["beta"]) * log((1 + exp(para["alpha"] + para["rho"] + para["beta"] * (ages + 1)))
                                             / (1 + exp(para["alpha"] + para["rho"] + para["beta"] * ages)))
  },
  makehamPerksFunc = function(para, ages) {
    exp(para["epsilon"]) + ((1 - exp(para["epsilon"])) / para["beta"]) * log((1 + exp(para["alpha"] + para["beta"] * (ages + 1)))
                                                                             / (1 + exp(para["alpha"] + para["beta"] * ages)))
  },
  makehamBeardFunc = function(para, ages) {
    exp(para["epsilon"]) + ((exp(-para["rho"]) - exp(para["epsilon"])) / para["beta"]) *
      log((1 + exp(para["alpha"] + para["rho"] + para["beta"] * (ages + 1))) /
            (1 + exp(para["alpha"] + para["rho"] + para["beta"] * ages)))
  }
)

HAZARD_FUNCTIONS <- list(
  gompertzFunc = function(para, ages) {
    exp(para["alpha"] + para["beta"] * ages)
  },
  makehamFunc = function(para, ages) {
    exp(para["epsilon"]) + exp(para["alpha"] + para["beta"] * ages)
  },
  perksFunc = function(para, ages) {
    exp(para["alpha"] + para["beta"] * ages) / (1 + exp(para["alpha"] + para["beta"] * ages))
  },
  beardFunc = function(para, ages) {
    exp(para["alpha"] + para["beta"] * ages) /
      (1 + exp(para["alpha"] + para["rho"] + para["beta"] * ages))
  },
  makehamPerksFunc = function(para, ages) {
    (exp(para["epsilon"]) + exp(para["alpha"] + para["beta"] * ages)) / (1 + exp(para["alpha"] + para["beta"] * ages))
  },
  makehamBeardFunc = function(para, ages) {
    (exp(para["epsilon"]) + exp(para["alpha"] + para["beta"] * ages)) /
      (1 + exp(para["alpha"] + para["rho"] + para["beta"] * ages))
  }
)




##### GRADUATE RESULT ANALYSIS #####

#' Statistical Analysis of a Graduation Result
#'
#' @param resultList
createGradAnalysis <- function(resultList) {
  # Required packages
  require(dplyr)
  
  # Make sure data is in the right order and extract Vectors
  gradData <- resultList$gradData %>%
    arrange(.data$Age)
  gradResults <- resultList$gradResults
  hessian <- resultList$hessian
  rates <- resultList$rates
  
  ages <- gradData$Age
  decs <- gradData$Deaths
  exps <- gradData$Exposure
  
  
  logLike <- gradResults$value
  fittedParameters <- gradResults$par
  nPars <- length(fittedParameters)
  
  # Analysis
  paraSD <- tryCatch({
    sqrt(diag(solve(-hessian, tol = 1e-20)))
  },
  error = function(con) {
    print("No inverse of Hessian calculated")
    print(con)
    return(NULL)
  })
  
  ZScore <- if (!is.null(paraSD)) {
    fittedParameters / paraSD
  } else {
    NULL
  }
  
  pVal <- if (!is.null(ZScore)) {
    2 * pnorm(abs(ZScore), 0, 1, lower.tail = FALSE)
  } else {
    NULL
  }
  
  # Calculate deviance and residuals
  A <- decs
  E <- rates * exps
  dev <- sum(2 * (ifelse(A == 0, 0, A * log(A / E)) - (A - E)))
  devRes <-
    sign(A - E) * sqrt(2 * (ifelse(A == 0, 0, A * log(A / E)) - (A - E)))
  stRes <- (A - E) / sqrt(E)
  
  # Determine dispersion coefficient
  dis <- dev / (length(ages) - nPars)
  
  # Chi-squared calculations
  chi <- sum(((A - E) ^ 2) / ifelse(E == 0, 1, E))
  chiP <- pchisq(chi, df = length(ages) - nPars, lower.tail = FALSE)
  
  # Information criteria
  AIC <- -2 * logLike + 2 * nPars
  BIC <- -2 * logLike + log(length(ages)) * nPars
  
  # Signs test
  if (!any(is.na(devRes))) {
    signsP <- sum(devRes > 0)
    signsN <- sum(devRes < 0)
    signsA <- if (signsP <= signsN) {
      "less"
    } else {
      "greater"
    }
    signsTest <- binom.test(signsP, signsP + signsN,
                            alternative = signsA)$p.value
  } else {
    signsP <- NULL
    signsN <- NULL
    signsTest <- NULL
  }
  
  
  # Runs test
  runs <- length(rle(sign(devRes))$lengths)
  runsP <-
    randtests::runs.test(devRes, "left.sided", 0, "exact", FALSE)$p.value
  
  # Calculate 95% confidence intervals
  lower <- log(qpois(0.025, A) / exps)
  upper <- log(qpois(0.975, A) / exps)
  
  modelResults <-
    list(
      Parameters = split(unname(fittedParameters), names(fittedParameters)),
      SD = if (!is.null(paraSD)) {
        split(unname(paraSD), names(fittedParameters))
      },
      Hessian = hessian,
      ZStat = if (!is.null(paraSD)) {
        split(unname(ZScore), names(fittedParameters))
      },
      PVals = if (!is.null(paraSD)) {
        split(unname(pVal), names(fittedParameters))
      },
      LogLikelihood = logLike,
      AIC = AIC,
      BIC = BIC,
      dis = dis,
      chi = chi,
      chiP = chiP,
      signsP = signsP,
      signsN = signsN,
      signsTest = signsTest,
      runs = runs,
      runsP = runsP
    )
  
  # Plot data crude versus fitted central mortality Rates
  obs <- log(decs / exps)
  mod <- log(rates)
  
  # if(length(mod) == 1) mod <- rep(mod, length(obs))
  
  plottingData <- data.frame(Age = ages, upper, lower, obs, mod, model = resultList$gradFunc$name)
  residualData <- data.frame(Age = ages, devRes, stRes, model = resultList$gradFunc$name)
  
  FinalOutput <-
    list(
      gradData,
      modelResults,
      plottingData,
      residualData,
      gradFunc = resultList$gradFunc
    )
  names(FinalOutput) <-
    c("gradData",
      "modelResults",
      "plottingData",
      "residualData",
      "gradFunc")
  return(FinalOutput)
}


applyFittedCLAModelOnData <- function(resultList, newData, analysis = TRUE) {
  # Required packages
  require(dplyr)
  
  hessian <- resultList$hessian
  
  # Make sure data is in the right order and extract Vectors
  gradData <- newData %>%
    arrange(.data$Age)
  
  ages <- gradData$Age
  decs <- gradData$Deaths
  exps <- gradData$Exposure
  
  rates <- calculateRatesCLAModel(resultList, ages)
  # Plot data crude versus fitted central mortality Rates
  obs <- log(decs / exps)
  mod <- log(rates)
  
  if(!analysis) {
    plottingData <- data.frame(Age = ages, obs, mod)
    return(
      list(
        gradData = gradData,
        plottingData = plottingData
      )
    )
  }
  
  logLike <- resultList$gradResults$value
  fittedParameters <- resultList$gradResults$par
  nPars <- length(fittedParameters)
  
  # Analysis
  paraSD <- tryCatch({
    sqrt(diag(solve(-resultList$hessian, tol = 1e-20)))
  },
  error = function(con) {
    print("No inverse of Hessian calculated")
    print(con)
    return(NULL)
  })
  
  ZScore <- if (!is.null(paraSD)) {
    fittedParameters / paraSD
  } else {
    NULL
  }
  
  pVal <- if (!is.null(ZScore)) {
    2 * pnorm(abs(ZScore), 0, 1, lower.tail = FALSE)
  } else {
    NULL
  }
  
  # Calculate deviance and residuals
  A <- decs
  E <- rates * exps
  dev <- sum(2 * (ifelse(A == 0, 0, A * log(A / E)) - (A - E)))
  devRes <-
    sign(A - E) * sqrt(2 * (ifelse(A == 0, 0, A * log(A / E)) - (A - E)))
  stRes <- (A - E) / sqrt(E)
  
  # Determine dispersion coefficient
  dis <- dev / (length(ages) - nPars)
  
  # Chi-squared calculations
  chi <- sum(((A - E) ^ 2) / ifelse(E == 0, 1, E))
  chiP <- pchisq(chi, df = length(ages) - nPars, lower.tail = FALSE)
  
  # Information criteria
  AIC <- -2 * logLike + 2 * nPars
  BIC <- -2 * logLike + log(length(ages)) * nPars
  
  # Signs test
  if (!any(is.na(devRes))) {
    signsP <- sum(devRes > 0)
    signsN <- sum(devRes < 0)
    signsA <- if (signsP <= signsN) {
      "less"
    } else {
      "greater"
    }
    signsTest <- binom.test(signsP, signsP + signsN,
                            alternative = signsA)$p.value
  } else {
    signsP <- NULL
    signsN <- NULL
    signsTest <- NULL
  }
  
  
  # Runs test
  runs <- length(rle(sign(devRes))$lengths)
  runsP <-
    randtests::runs.test(devRes, "left.sided", 0, "exact", FALSE)$p.value
  
  # Calculate 95% confidence intervals
  lower <- log(qpois(0.025, A) / exps)
  upper <- log(qpois(0.975, A) / exps)
  
  modelResults <-
    list(
      Parameters = split(unname(fittedParameters), names(fittedParameters)),
      SD = if (!is.null(paraSD)) {
        split(unname(paraSD), names(fittedParameters))
      },
      Hessian = resultList$hessian,
      ZStat = if (!is.null(paraSD)) {
        split(unname(ZScore), names(fittedParameters))
      },
      PVals = if (!is.null(paraSD)) {
        split(unname(pVal), names(fittedParameters))
      },
      LogLikelihood = logLike,
      AIC = AIC,
      BIC = BIC,
      dis = dis,
      chi = chi,
      chiP = chiP,
      signsP = signsP,
      signsN = signsN,
      signsTest = signsTest,
      runs = runs,
      runsP = runsP
    )
  
  # if(length(mod) == 1) mod <- rep(mod, length(obs))
  
  plottingData <- data.frame(Age = ages, upper, lower, obs, mod, model = resultList$gradFunc$name)
  residualData <- data.frame(Age = ages, devRes, stRes, model = resultList$gradFunc$name)
  
  FinalOutput <-
    list(
      gradData,
      modelResults,
      plottingData,
      residualData,
      gradFunc = resultList$gradFunc
    )
  names(FinalOutput) <-
    c("gradData",
      "modelResults",
      "plottingData",
      "residualData",
      "gradFunc")
  return(FinalOutput)
}

calculateRatesCLAModel <- function(resultList, ages) {
  return(resultList$gradFunc$func(unlist(resultList$gradResults$par), ages))
}


applyFittedEVTModelOnData <- function(resultList, newData, analysis = TRUE) {
  # Required packages
  require(dplyr)
  
  hessian <- resultList$hessian
  
  # Make sure data is in the right order and extract Vectors
  gradData <- newData %>%
    arrange(.data$Age)
  
  ages <- gradData$Age
  decs <- gradData$Deaths
  exps <- gradData$Exposure
  
  rates <- calculateRatesEVTModel(resultList, ages)
  # Plot data crude versus fitted central mortality Rates
  obs <- log(decs / exps)
  mod <- log(rates)
  
  if(!analysis) {
    plottingData <- data.frame(Age = ages, obs, mod)
    return(
      list(
        gradData = gradData,
        plottingData = plottingData
      )
    )
  }
  
  logLike <- resultList$gradResults$value
  fittedParameters <- resultList$gradResults$par
  nPars <- length(fittedParameters)
  
  # Analysis
  paraSD <- tryCatch({
    sqrt(diag(solve(-resultList$hessian, tol = 1e-20)))
  },
  error = function(con) {
    print("No inverse of Hessian calculated")
    print(con)
    return(NULL)
  })
  
  ZScore <- if (!is.null(paraSD)) {
    fittedParameters / paraSD
  } else {
    NULL
  }
  
  pVal <- if (!is.null(ZScore)) {
    2 * pnorm(abs(ZScore), 0, 1, lower.tail = FALSE)
  } else {
    NULL
  }
  
  # Calculate deviance and residuals
  A <- decs
  E <- rates * exps
  dev <- sum(2 * (ifelse(A == 0, 0, A * log(A / E)) - (A - E)))
  devRes <-
    sign(A - E) * sqrt(2 * (ifelse(A == 0, 0, A * log(A / E)) - (A - E)))
  stRes <- (A - E) / sqrt(E)
  
  # Determine dispersion coefficient
  dis <- dev / (length(ages) - nPars)
  
  # Chi-squared calculations
  chi <- sum(((A - E) ^ 2) / ifelse(E == 0, 1, E))
  chiP <- pchisq(chi, df = length(ages) - nPars, lower.tail = FALSE)
  
  # Information criteria
  AIC <- -2 * logLike + 2 * nPars
  BIC <- -2 * logLike + log(length(ages)) * nPars
  
  # Signs test
  if (!any(is.na(devRes))) {
    signsP <- sum(devRes > 0)
    signsN <- sum(devRes < 0)
    signsA <- if (signsP <= signsN) {
      "less"
    } else {
      "greater"
    }
    signsTest <- binom.test(signsP, signsP + signsN,
                            alternative = signsA)$p.value
  } else {
    signsP <- NULL
    signsN <- NULL
    signsTest <- NULL
  }
  
  
  # Runs test
  runs <- length(rle(sign(devRes))$lengths)
  runsP <-
    randtests::runs.test(devRes, "left.sided", 0, "exact", FALSE)$p.value
  
  # Calculate 95% confidence intervals
  lower <- log(qpois(0.025, A) / exps)
  upper <- log(qpois(0.975, A) / exps)
  
  modelResults <-
    list(
      Parameters = split(unname(fittedParameters), names(fittedParameters)),
      SD = if (!is.null(paraSD)) {
        split(unname(paraSD), names(fittedParameters))
      },
      Hessian = resultList$hessian,
      ZStat = if (!is.null(paraSD)) {
        split(unname(ZScore), names(fittedParameters))
      },
      PVals = if (!is.null(paraSD)) {
        split(unname(pVal), names(fittedParameters))
      },
      LogLikelihood = logLike,
      AIC = AIC,
      BIC = BIC,
      dis = dis,
      chi = chi,
      chiP = chiP,
      signsP = signsP,
      signsN = signsN,
      signsTest = signsTest,
      runs = runs,
      runsP = runsP
    )
  
  # if(length(mod) == 1) mod <- rep(mod, length(obs))
  
  plottingData <- data.frame(Age = ages, upper, lower, obs, mod, model = resultList$gradFunc$name)
  residualData <- data.frame(Age = ages, devRes, stRes, model = resultList$gradFunc$name)
  
  FinalOutput <-
    list(
      gradData,
      modelResults,
      plottingData,
      residualData,
      gradFunc = resultList$gradFunc
    )
  names(FinalOutput) <-
    c("gradData",
      "modelResults",
      "plottingData",
      "residualData",
      "gradFunc")
  return(FinalOutput)
}

calculateRatesEVTModel <- function(resultList, ages) {
  agesB <- sort(subset(ages, ages < resultList$thresholdAge))
  agesA <- sort(subset(ages, ages >= resultList$thresholdAge))
  
  rates <- c(
    resultList$gradFunc$hazardFunc(unlist(resultList$gradResults$par), ages = agesB),
    resultList$gradFunc$hazardFuncGPD(unlist(resultList$gradResults$par), ages = agesA)
  )
  return(rates)
}



##### CLASSIC GRADUATE FUNCTIONS #####
graduatePoissonCLA <- function(data,
                               model,
                               initParameter = INIT_PARAMETER_CLA,
                               constraintPara = OPTIM_CONSTRAINTS_CLA,
                               optimControl = OPTIM_CONTROL_PARAMETERS,
                               integratedHazard = TRUE,
                               analysis = TRUE) {
  if (model == "Gompertz") {
    graduateGompertz(data,
                     initParameter,
                     constraintPara,
                     optimControl,
                     integratedHazard,
                     analysis)
  }
  if (model == "Makeham") {
    graduateMakeham(data,
                    initParameter,
                    constraintPara,
                    optimControl,
                    integratedHazard,
                    analysis)
  }
  if (model == "Perks") {
    graduatePerks(data,
                  initParameter,
                  constraintPara,
                  optimControl,
                  integratedHazard,
                  analysis)
  }
  if (model == "Beard") {
    graduateBeard(data,
                  initParameter,
                  constraintPara,
                  optimControl,
                  integratedHazard,
                  analysis)
  }
  if (model == "MakehamPerks") {
    graduateMakehamPerks(data,
                         initParameter,
                         constraintPara,
                         optimControl,
                         integratedHazard,
                         analysis)
  }
  if (model == "MakehamBeard") {
    graduateMakehamBeard(data,
                         initParameter,
                         constraintPara,
                         optimControl,
                         integratedHazard,
                         analysis)
  }
}


graduateGompertz <- function(gradData,
                             initParameter = INIT_PARAMETER_CLA,
                             constraintPara = OPTIM_CONSTRAINTS_CLA,
                             optimControl = OPTIM_CONTROL_PARAMETERS,
                             integratedHazard = TRUE,
                             analysis = TRUE) {
  # Required packages
  require(numDeriv)
  require(dplyr)
  
  # Make sure data is in the right order and extract Vectors
  gradData <- gradData %>%
    arrange(.data$Age)
  
  ages <- gradData$Age
  decs <- gradData$Deaths
  exps <- gradData$Exposure
  
  # Initial Parameters
  paraInit <- c(alpha = initParameter$initAlpha,
                beta = initParameter$initBeta)
  
  lowerConstraint <- c(constraintPara$alphaLower,
                       constraintPara$betaLower)
  
  upperConstraint <- c(constraintPara$alphaUpper,
                       constraintPara$betaUpper)
  
  # Choose hazard function
  if (integratedHazard) {
    hazardFunc <- INTEGRATED_HAZARD_FUNCTIONS$gompertzFunc
  } else if (!integratedHazard) {
    hazardFunc <- HAZARD_FUNCTIONS$gompertzFunc
  }
  
  
  # Optim
  logLikeFunc <- function(para) {
    mu <-
      pmax(optimControl$minimumMu, hazardFunc(para = para, ages = ages))
    sum(ifelse(exps == 0, 0,
               -exps * mu + decs * log(exps * mu) - lfactorial(decs)))
  }
  
  grr <- function(para) {
    numDeriv::grad(logLikeFunc, para, method = optimControl$gradMethod)
  }
  
  gradResults <- optim(
    par = paraInit,
    gr = grr,
    fn = logLikeFunc,
    method = optimControl$optimMethod,
    lower = lowerConstraint,
    upper = upperConstraint,
    control = list(
      fnscale = -1,
      reltol = optimControl$relTolerance,
      maxit = optimControl$maxIteration
    )
  )
  
  calcHessian <-
    numDeriv::hessian(logLikeFunc, gradResults$par, method = optimControl$gradMethod)
  
  # TODO This should not be here
  rates <- hazardFunc(para = gradResults$par, ages = ages)
  
  result <-
    list(
      gradData = gradData,
      gradResults = gradResults,
      hessian = calcHessian,
      gradFunc = list(name = "Gompertz", func = hazardFunc),
      rates = rates
    )
  
  if (!analysis) {
    return(result)
  }
  
  return(createGradAnalysis(result))
}

graduateMakeham <- function(gradData,
                            initParameter = INIT_PARAMETER_CLA,
                            constraintPara = OPTIM_CONSTRAINTS_CLA,
                            optimControl = OPTIM_CONTROL_PARAMETERS,
                            integratedHazard = TRUE,
                            analysis = TRUE) {
  # Call Gompertz for initial Parameters
  result <- graduateGompertz(
    gradData,
    initParameter,
    constraintPara,
    optimControl,
    integratedHazard,
    analysis = FALSE
  )
  
  # Required packages
  require(numDeriv)
  require(dplyr)
  
  # Make sure data is in the right order and extract Vectors
  gradData <- gradData %>%
    arrange(.data$Age)
  
  ages <- gradData$Age
  decs <- gradData$Deaths
  exps <- gradData$Exposure
  
  
  # Initial Parameters
  paraInit <- c(result$gradResult$par["alpha"],
                result$gradResult$par["beta"],
                epsilon = initParameter$initEpsilon)
  
  lowerConstraint <- c(
    constraintPara$alphaLower,
    constraintPara$betaLower,
    constraintPara$epsilonLower
  )
  upperConstraint <- c(
    constraintPara$alphaUpper,
    constraintPara$betaUpper,
    constraintPara$epsilonUpper
  )
  
  # Choose hazard function
  if (integratedHazard) {
    hazardFunc <- INTEGRATED_HAZARD_FUNCTIONS$makehamFunc
  } else if (!integratedHazard) {
    hazardFunc <- HAZARD_FUNCTIONS$makehamFunc
  }
  
  # Optim
  logLikeFunc <- function(para) {
    mu <-
      pmax(optimControl$minimumMu, hazardFunc(para = para, ages = ages))
    sum(ifelse(exps == 0, 0,
               -exps * mu + decs * log(exps * mu) - lfactorial(decs)))
  }
  
  grr <- function(para) {
    numDeriv::grad(logLikeFunc, para, method = optimControl$gradMethod)
  }
  
  gradResults <- optim(
    par = paraInit,
    gr = grr,
    fn = logLikeFunc,
    method = optimControl$optimMethod,
    lower = lowerConstraint,
    upper = upperConstraint,
    control = list(
      fnscale = -1,
      reltol = optimControl$relTolerance,
      maxit = optimControl$maxIteration
    )
  )
  
  calcHessian <-
    numDeriv::hessian(logLikeFunc, gradResults$par, method = optimControl$gradMethod)
  
  # TODO This should not be here
  rates <- hazardFunc(para = gradResults$par, ages = ages)
  
  result <-
    list(
      gradData = gradData,
      gradResults = gradResults,
      hessian = calcHessian,
      gradFunc = list(name = "Makeham", func = hazardFunc),
      rates = rates
    )
  
  if (!analysis) {
    return(result)
  }
  
  return(createGradAnalysis(result))
}


graduatePerks <- function(gradData,
                          initParameter = INIT_PARAMETER_CLA,
                          constraintPara = OPTIM_CONSTRAINTS_CLA,
                          optimControl = OPTIM_CONTROL_PARAMETERS,
                          integratedHazard = TRUE,
                          analysis = TRUE) {
  # Required packages
  require(numDeriv)
  require(dplyr)
  
  # Make sure data is in the right order and extract Vectors
  gradData <- gradData %>%
    arrange(.data$Age)
  
  ages <- gradData$Age
  decs <- gradData$Deaths
  exps <- gradData$Exposure
  
  # Initial Parameters
  paraInit <- c(alpha = initParameter$initAlpha,
                beta = initParameter$initBeta)
  
  lowerConstraint <- c(constraintPara$alphaLower,
                       constraintPara$betaLower)
  upperConstraint <- c(constraintPara$alphaUpper,
                       constraintPara$betaUpper)
  
  # Choose hazard function
  if (integratedHazard) {
    hazardFunc <- INTEGRATED_HAZARD_FUNCTIONS$perksFunc
  } else if (!integratedHazard) {
    hazardFunc <- HAZARD_FUNCTIONS$perksFunc
  }
  
  # Optim
  logLikeFunc <- function(para) {
    mu <-
      pmax(optimControl$minimumMu, hazardFunc(para = para, ages = ages))
    sum(ifelse(exps == 0, 0,
               -exps * mu + decs * log(exps * mu) - lfactorial(decs)))
  }
  
  grr <- function(para) {
    numDeriv::grad(logLikeFunc, para, method = optimControl$gradMethod)
  }
  
  gradResults <- optim(
    par = paraInit,
    gr = grr,
    fn = logLikeFunc,
    method = optimControl$optimMethod,
    lower = lowerConstraint,
    upper = upperConstraint,
    control = list(
      fnscale = -1,
      reltol = optimControl$relTolerance,
      maxit = optimControl$maxIteration
    )
  )
  
  calcHessian <-
    numDeriv::hessian(logLikeFunc, gradResults$par, method = optimControl$gradMethod)
  
  # TODO This should not be here
  rates <- hazardFunc(para = gradResults$par, ages = ages)
  
  result <-
    list(
      gradData = gradData,
      gradResults = gradResults,
      hessian = calcHessian,
      gradFunc = list(name = "Perks", func = hazardFunc),
      rates = rates
    )
  
  if (!analysis) {
    return(result)
  }
  
  return(createGradAnalysis(result))
}

graduateBeard <- function(gradData,
                          initParameter = INIT_PARAMETER_CLA,
                          constraintPara = OPTIM_CONSTRAINTS_CLA,
                          optimControl = OPTIM_CONTROL_PARAMETERS,
                          integratedHazard = TRUE,
                          analysis = TRUE) {
  # Call Perks for initial Parameters
  result <- graduatePerks(
    gradData,
    initParameter,
    constraintPara,
    optimControl,
    integratedHazard,
    analysis = FALSE
  )
  
  # Required packages
  require(numDeriv)
  require(dplyr)
  
  # Make sure data is in the right order and extract Vectors
  gradData <- gradData %>%
    arrange(.data$Age)
  
  ages <- gradData$Age
  decs <- gradData$Deaths
  exps <- gradData$Exposure
  
  # Initial Parameters
  paraInit <- c(result$gradResult$par["alpha"],
                result$gradResult$par["beta"],
                rho = initParameter$initRho)
  
  lowerConstraint <- c(constraintPara$alphaLower,
                       constraintPara$betaLower,
                       constraintPara$rhoLower)
  upperConstraint <- c(constraintPara$alphaUpper,
                       constraintPara$betaUpper,
                       constraintPara$rhoUpper)
  
  # Choose hazard function
  if (integratedHazard) {
    hazardFunc <- INTEGRATED_HAZARD_FUNCTIONS$beardFunc
  } else if (!integratedHazard) {
    hazardFunc <- HAZARD_FUNCTIONS$beardFunc
  }
  
  # Optim
  logLikeFunc <- function(para) {
    mu <-
      pmax(optimControl$minimumMu, hazardFunc(para = para, ages = ages))
    sum(ifelse(exps == 0, 0,
               -exps * mu + decs * log(exps * mu) - lfactorial(decs)))
  }
  
  grr <- function(para) {
    numDeriv::grad(logLikeFunc, para, method = optimControl$gradMethod)
  }
  
  gradResults <- optim(
    par = paraInit,
    gr = grr,
    fn = logLikeFunc,
    method = optimControl$optimMethod,
    lower = lowerConstraint,
    upper = upperConstraint,
    control = list(
      fnscale = -1,
      reltol = optimControl$relTolerance,
      maxit = optimControl$maxIteration
    )
  )
  
  calcHessian <-
    numDeriv::hessian(logLikeFunc, gradResults$par, method = optimControl$gradMethod)
  
  # TODO This should not be here
  rates <- hazardFunc(para = gradResults$par, ages = ages)
  
  result <-
    list(
      gradData = gradData,
      gradResults = gradResults,
      hessian = calcHessian,
      gradFunc = list(name = "Beard", func = hazardFunc),
      rates = rates
    )
  
  if (!analysis) {
    return(result)
  }
  
  return(createGradAnalysis(result))
}

graduateMakehamPerks <- function(gradData,
                                 initParameter = INIT_PARAMETER_CLA,
                                 constraintPara = OPTIM_CONSTRAINTS_CLA,
                                 optimControl = OPTIM_CONTROL_PARAMETERS,
                                 integratedHazard = TRUE,
                                 analysis = TRUE) {
  # Call Perks for initial Parameters
  result <- graduatePerks(
    gradData,
    initParameter,
    constraintPara,
    optimControl,
    integratedHazard,
    analysis = FALSE
  )
  
  # Required packages
  require(numDeriv)
  require(dplyr)
  
  # Make sure data is in the right order and extract Vectors
  gradData <- gradData %>%
    arrange(.data$Age)
  
  ages <- gradData$Age
  decs <- gradData$Deaths
  exps <- gradData$Exposure
  
  
  # Initial Parameters
  paraInit <- c(result$gradResult$par["alpha"],
                result$gradResult$par["beta"],
                epsilon = initParameter$initEpsilon)
  
  lowerConstraint <- c(
    constraintPara$alphaLower,
    constraintPara$betaLower,
    constraintPara$epsilonLower
  )
  upperConstraint <- c(
    constraintPara$alphaUpper,
    constraintPara$betaUpper,
    constraintPara$epsilonUpper
  )
  
  # Choose hazard function
  if (integratedHazard) {
    hazardFunc <- INTEGRATED_HAZARD_FUNCTIONS$makehamPerksFunc
  } else if (!integratedHazard) {
    hazardFunc <- HAZARD_FUNCTIONS$makehamPerksFunc
  }
  
  
  # Optim
  logLikeFunc <- function(para) {
    mu <-
      pmax(optimControl$minimumMu, hazardFunc(para = para, ages = ages))
    sum(ifelse(exps == 0, 0,
               -exps * mu + decs * log(exps * mu) - lfactorial(decs)))
  }
  
  grr <- function(para) {
    numDeriv::grad(logLikeFunc, para, method = optimControl$gradMethod)
  }
  
  gradResults <- optim(
    par = paraInit,
    gr = grr,
    fn = logLikeFunc,
    method = optimControl$optimMethod,
    lower = lowerConstraint,
    upper = upperConstraint,
    control = list(
      fnscale = -1,
      reltol = optimControl$relTolerance,
      maxit = optimControl$maxIteration
    )
  )
  
  calcHessian <-
    numDeriv::hessian(logLikeFunc, gradResults$par, method = optimControl$gradMethod)
  
  # TODO This should not be here
  rates <- hazardFunc(para = gradResults$par, ages = ages)
  
  result <-
    list(
      gradData = gradData,
      gradResults = gradResults,
      hessian = calcHessian,
      gradFunc = list(name = "MakehamPerks", func = hazardFunc),
      rates = rates
    )
  
  if (!analysis) {
    return(result)
  }
  
  return(createGradAnalysis(result))
}

graduateMakehamBeard <- function(gradData,
                                 initParameter = INIT_PARAMETER_CLA,
                                 constraintPara = OPTIM_CONSTRAINTS_CLA,
                                 optimControl = OPTIM_CONTROL_PARAMETERS,
                                 integratedHazard = TRUE,
                                 analysis = TRUE) {
  # Call MakehamPerks for initial Parameters
  result <- graduateMakehamPerks(
    gradData,
    initParameter,
    constraintPara,
    optimControl,
    integratedHazard,
    analysis = FALSE
  )
  
  # Required packages
  require(numDeriv)
  require(dplyr)
  
  # Make sure data is in the right order and extract Vectors
  gradData <- gradData %>%
    arrange(.data$Age)
  
  ages <- gradData$Age
  decs <- gradData$Deaths
  exps <- gradData$Exposure
  
  # Initial Parameters
  paraInit <- c(
    result$gradResult$par["alpha"],
    result$gradResult$par["beta"],
    result$gradResult$par["epsilon"],
    rho = initParameter$initRho
  )
  
  lowerConstraint <- c(
    constraintPara$alphaLower,
    constraintPara$betaLower,
    constraintPara$epsilonLower,
    constraintPara$rhoLower
  )
  upperConstraint <- c(
    constraintPara$alphaUpper,
    constraintPara$betaUpper,
    constraintPara$epsilonUpper,
    constraintPara$rhoUpper
  )
  
  # Choose hazard function
  if (integratedHazard) {
    hazardFunc <- INTEGRATED_HAZARD_FUNCTIONS$makehamBeardFunc
  } else if (!integratedHazard) {
    hazardFunc <- INTEGRATED_HAZARD_FUNCTIONS$makehamBeardFunc
  }
  
  # Optim
  logLikeFunc <- function(para) {
    mu <-
      pmax(optimControl$minimumMu, hazardFunc(para = para, ages = ages))
    sum(ifelse(exps == 0, 0,
               -exps * mu + decs * log(exps * mu) - lfactorial(decs)))
  }
  
  grr <- function(para) {
    numDeriv::grad(logLikeFunc, para, method = optimControl$gradMethod)
  }
  
  gradResults <- optim(
    par = paraInit,
    gr = grr,
    fn = logLikeFunc,
    method = optimControl$optimMethod,
    lower = lowerConstraint,
    upper = upperConstraint,
    control = list(
      fnscale = -1,
      reltol = optimControl$relTolerance,
      maxit = optimControl$maxIteration
    )
  )
  
  calcHessian <-
    numDeriv::hessian(logLikeFunc, gradResults$par, method = optimControl$gradMethod)
  
  # TODO This should not be here
  rates <- hazardFunc(para = gradResults$par, ages = ages)
  
  result <-
    list(
      gradData = gradData,
      gradResults = gradResults,
      hessian = calcHessian,
      gradFunc = list(name = "MakehamBeard", func = hazardFunc),
      rates = rates
    )
  
  if (!analysis) {
    return(result)
  }
  
  return(createGradAnalysis(result))
}





##### HERMITE GRADUATE FUNCTIONS #####
graduateHermite <- function(gradData,
                            type,
                            X0 = 40,
                            X1 = 110,
                            initParameter = INIT_PARAMETER_HMT,
                            constraintPara = OPTIM_CONSTRAINTS_HMT,
                            optimControl = OPTIM_CONTROL_PARAMETERS,
                            integratedHazard = TRUE,
                            analysis = TRUE) {
  # Required packages
  require(numDeriv)
  require(dplyr)
  
  # Make sure data is in the right order and extract Vectors
  gradData <- gradData %>%
    arrange(.data$Age)
  
  ages <- gradData$Age
  decs <- gradData$Deaths
  exps <- gradData$Exposure
  
  # Initial Parameters
  paraInit <-
    c(hermiteAlpha = INIT_PARAMETER_HMT$hermiteAlpha,
      hermiteOmega = INIT_PARAMETER_HMT$hermiteOmega)
  if (type == "II" || type == "IV" || type == "V") {
    paraInit <- c(paraInit, hermiteM0 = INIT_PARAMETER_HMT$hermiteM0)
  }
  if (type == "III" || type == "IV") {
    paraInit <- c(paraInit, hermiteM1 = INIT_PARAMETER_HMT$hermiteM1)
  }
  if (type == "V") {
    paraInit <- c(paraInit, hermiteV = INIT_PARAMETER_HMT$hermiteV)
  }
  
  lowerConstraint <- c(
    constraintPara$hermiteAlphaLower,
    constraintPara$hermiteOmegaLower,
    constraintPara$hermiteM0Lower,
    constraintPara$hermiteM1Lower
  )
  upperConstraint <- c(
    constraintPara$hermiteAlphaUpper,
    constraintPara$hermiteOmegaUpper,
    constraintPara$hermiteM0Upper,
    constraintPara$hermiteM1Upper
  )
  
  
  # Hermite Functions
  scaleAge <- function(x, x0, x1) {
    (x - x0) / (x1 - x0)
  }
  h00x <- function(x, x0, x1) {
    t <- scaleAge(x, x0, x1)
    return((1.0 + 2 * t) * (1.0 - t) * (1.0 - t))
  }
  h10x <- function(x, x0, x1) {
    t <- scaleAge(x, x0, x1)
    return(t * (1.0 - t) * (1.0 - t))
  }
  h01x <- function(x, x0, x1) {
    t <- scaleAge(x, x0, x1)
    return(t * t * (3.0 - 2 * t))
  }
  h11x <- function(x, x0, x1) {
    t <- scaleAge(x, x0, x1)
    return(t * t * (t - 1.0))
  }
  hQuarticx <- function(x, x0, x1) {
    t <- scaleAge(x, x0, x1)
    return(16 * t ^ 2 * (1 - t) ^ 2)
  }
  
  pH <- function(para) {
    if (type == "I") {
      return(c(para, 0, 0, 0))
    }
    if (type == "II") {
      return(c(para, 0, 0))
    }
    if (type == "III") {
      return(c(para[1:2], 0, para[3], 0))
    }
    if (type == "IV") {
      return(c(para, 0))
    }
    if (type == "V") {
      return(c(para[1:3], 0, para[4]))
    }
    return(para)
  }
  
  
  hermiteMu <- function(para, x) {
    return(exp(
      para[1] * h00x(x, X0, X1) + para[2] * h01x(x, X0, X1) +
        para[3] * h10x(x, X0, X1) + para[4] * h11x(x, X0, X1) +
        para[5] * hQuarticx(x, X0, X1)
    ))
  }
  
  logLikeFunc <- function(para) {
    mu <-
      pmax(optimControl$minimumMu, hermiteMu(para = pH(para), x = ages))
    sum(ifelse(exps == 0, 0,
               -exps * mu + decs * log(exps * mu) - lfactorial(decs)))
  }
  
  grr <- function(para) {
    numDeriv::grad(logLikeFunc, para, method = optimControl$gradMethod)
  }
  
  gradResults <- optim(
    par = paraInit,
    fn = logLikeFunc,
    gr = grr,
    method = optimControl$optimMethod,
    lower = lowerConstraint,
    upper = upperConstraint,
    control = list(
      fnscale = -1,
      reltol = optimControl$relTolerance,
      maxit = optimControl$maxIteration
    )
  )
  
  calcHessian <-
    numDeriv::hessian(logLikeFunc, gradResults$par, method = optimControl$gradMethod)
  
  # TODO This should not be here
  rates <- hermiteMu(p = pH(gradResults$par), x = ages)
  
  result <-
    list(
      gradData = gradData,
      modelSpecifications = list(type = type, X0 = X0, X1 = X1),
      gradResults = gradResults,
      hessian = calcHessian,
      gradFunc = list(name = paste("Hermite ", type), func = hermiteMu),
      rates = rates
    )
  
  if (!analysis) {
    return(result)
  }
  
  return(createGradAnalysis(result))
}


##### EVT GRADUATE FUNCTIONS #####
graduateGompertzGPD <- function(gradData,
                                thresholdAge,
                                initParameter = INIT_PARAMETER_GPD,
                                constraintPara = OPTIM_CONSTRAINTS_CLA,
                                optimControl = OPTIM_CONTROL_PARAMETERS,
                                integratedHazard = TRUE,
                                analysis = TRUE) {
  # Required packages
  require(numDeriv)
  require(dplyr)
  
  # Make sure data is in the right order and extract Vectors
  gradData <- gradData %>%
    arrange(.data$Age)
  
  dataBelowThreshold <-
    filter(gradData, .data$Age < thresholdAge) %>%
    arrange(Age)
  dataAboveThreshold <-
    filter(gradData, .data$Age >= thresholdAge) %>%
    arrange(Age)
  
  agesB <- dataBelowThreshold$Age
  decsB <- dataBelowThreshold$Deaths
  expsB <- dataBelowThreshold$Exposure
  
  agesA <- dataAboveThreshold$Age
  decsA <- dataAboveThreshold$Deaths
  expsA <- dataAboveThreshold$Exposure
  
  # Initial Parameters
  paraInit <- c(
    alpha = initParameter$initAlpha,
    beta = initParameter$initBeta,
    xi = initParameter$initXi
  )
  
  lowerConstraint <- c(constraintPara$alphaLower,
                       constraintPara$betaLower)
  
  upperConstraint <- c(constraintPara$alphaUpper,
                       constraintPara$betaUpper)
  
  
  # Choose hazard function
  if (integratedHazard) {
    hazardFunc <- INTEGRATED_HAZARD_FUNCTIONS$gompertzFunc
    hazardFuncGPD <- function(para, ages) {
      (log(abs(
        para["xi"] * (ages + 1 - thresholdAge) + exp(-para["alpha"] - para["beta"] * thresholdAge)
      )) / para["xi"]) -
        (log(abs(
          para["xi"] * (ages - thresholdAge) + exp(-para["alpha"] - para["beta"] * thresholdAge)
        )) / para["xi"])
    }
  } else if (!integratedHazard) {
    hazardFunc <- HAZARD_FUNCTIONS$gompertzFunc
    hazardFuncGPD <- function(para, ages) {
      1 / (exp(-para["alpha"] - para["beta"] * thresholdAge) + para["xi"] * (ages - thresholdAge))
    }
  }
  
  
  # Optim
  logLikeFunc <- function(para) {
    muGomp <-
      pmax(optimControl$minimumMu,
           hazardFunc(para = para, ages = agesB))
    muGPD <-
      pmax(optimControl$minimumMu,
           hazardFuncGPD(para = para, ages = agesA))
    sum(ifelse(
      expsB == 0,
      0,
      -expsB * muGomp + decsB * log(expsB * muGomp) - lfactorial(decsB)
    )) +
      sum(ifelse(
        expsA == 0,
        0,
        -expsA * muGPD + decsA * log(expsA * muGPD) - lfactorial(decsA)
      ))
  }
  
  grr <- function(para) {
    numDeriv::grad(logLikeFunc, para, method = optimControl$gradMethod)
  }
  
  gradResults <- optim(
    par = paraInit,
    gr = grr,
    fn = logLikeFunc,
    method = optimControl$optimMethod,
    lower = lowerConstraint,
    upper = upperConstraint,
    control = list(
      fnscale = -1,
      reltol = optimControl$relTolerance,
      maxit = optimControl$maxIteration
    )
  )
  
  calcHessian <-
    numDeriv::hessian(logLikeFunc, gradResults$par, method = optimControl$gradMethod)
  
  rates <- c(
    hazardFunc(para = gradResults$par, ages = agesB),
    hazardFuncGPD(para = gradResults$par, ages = agesA)
  )
  
  result <-
    list(
      gradData = gradData,
      gradResults = gradResults,
      hessian = calcHessian,
      gradFunc = list(name = "GompertzGPD", hazardFunc = hazardFunc, hazardFuncGPD = hazardFuncGPD),
      thresholdAge = thresholdAge,
      rates = rates
    )
  
  if (!analysis) {
    return(result)
  }
  
  return(createGradAnalysis(result))
}




graduateMakehamGPD <- function(gradData,
                               thresholdAge,
                               initParameter = INIT_PARAMETER_GPD,
                               constraintPara = OPTIM_CONSTRAINTS_CLA,
                               optimControl = OPTIM_CONTROL_PARAMETERS,
                               integratedHazard = TRUE,
                               analysis = TRUE) {
  # Required packages
  require(numDeriv)
  require(dplyr)
  
  # Make sure data is in the right order and extract Vectors
  gradData <- gradData %>%
    arrange(.data$Age)
  
  dataBelowThreshold <-
    filter(gradData, .data$Age < thresholdAge) %>%
    arrange(Age)
  dataAboveThreshold <-
    filter(gradData, .data$Age >= thresholdAge) %>%
    arrange(Age)
  
  agesB <- dataBelowThreshold$Age
  decsB <- dataBelowThreshold$Deaths
  expsB <- dataBelowThreshold$Exposure
  
  agesA <- dataAboveThreshold$Age
  decsA <- dataAboveThreshold$Deaths
  expsA <- dataAboveThreshold$Exposure
  
  # Initial Parameters
  paraInit <- c(
    alpha = initParameter$initAlpha,
    beta = initParameter$initBeta,
    epsilon = initParameter$initEpsilon,
    xi = initParameter$initXi
  )
  
  lowerConstraint <- c(constraintPara$alphaLower,
                       constraintPara$betaLower)
  
  upperConstraint <- c(constraintPara$alphaUpper,
                       constraintPara$betaUpper)
  
  
  # Choose hazard function
  if (integratedHazard) {
    hazardFunc <- INTEGRATED_HAZARD_FUNCTIONS$makehamFunc
    hazardFuncGPD <- function(para, ages) {
      (log(abs(
        para["xi"] * (ages + 1 - thresholdAge) + 1/(exp(para["epsilon"]) + exp(para["alpha"] + para["beta"] * thresholdAge))
      )) / para["xi"]) -
        (log(abs(
          para["xi"] * (ages - thresholdAge) + 1/(exp(para["epsilon"]) + exp(para["alpha"] + para["beta"] * thresholdAge))
        )) / para["xi"])
    }
  } else if (!integratedHazard) {
    hazardFunc <- HAZARD_FUNCTIONS$makehamFunc
    hazardFuncGPD <- function(para, ages) {
      1 / (1/(exp(para["epsilon"]) + exp(para["alpha"] + para["beta"] * thresholdAge)) + para["xi"] * (ages - thresholdAge))
    }
  }
  
  
  # Optim
  logLikeFunc <- function(para) {
    muGomp <-
      pmax(optimControl$minimumMu,
           hazardFunc(para = para, ages = agesB))
    muGPD <-
      pmax(optimControl$minimumMu,
           hazardFuncGPD(para = para, ages = agesA))
    sum(ifelse(
      expsB == 0,
      0,
      -expsB * muGomp + decsB * log(expsB * muGomp) - lfactorial(decsB)
    )) +
      sum(ifelse(
        expsA == 0,
        0,
        -expsA * muGPD + decsA * log(expsA * muGPD) - lfactorial(decsA)
      ))
  }
  
  grr <- function(para) {
    numDeriv::grad(logLikeFunc, para, method = optimControl$gradMethod)
  }
  
  gradResults <- optim(
    par = paraInit,
    gr = grr,
    fn = logLikeFunc,
    method = optimControl$optimMethod,
    lower = lowerConstraint,
    upper = upperConstraint,
    control = list(
      fnscale = -1,
      reltol = optimControl$relTolerance,
      maxit = optimControl$maxIteration
    )
  )
  
  calcHessian <-
    numDeriv::hessian(logLikeFunc, gradResults$par, method = optimControl$gradMethod)
  
  rates <- c(
    hazardFunc(para = gradResults$par, ages = agesB),
    hazardFuncGPD(para = gradResults$par, ages = agesA)
  )
  
  result <-
    list(
      gradData = gradData,
      gradResults = gradResults,
      hessian = calcHessian,
      gradFunc = list(name = "MakehamGPD", hazardFunc = hazardFunc, hazardFuncGPD = hazardFuncGPD),
      thresholdAge = thresholdAge,
      rates = rates
    )
  
  if (!analysis) {
    return(result)
  }
  
  return(createGradAnalysis(result))
}
