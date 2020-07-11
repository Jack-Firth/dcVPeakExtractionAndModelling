setwd("C:/Users/owner/Documents/Projects/dcVExtraction/")

allXAxisLabel = "Voltage / V"
allYAxisLabel = expression(paste("Current (", mu, "A)"))

separatePeaks = function(passedData) {
  # Converts loaded data into two separate columns for forward and reverse sweeps
  print("Separating forward and reverse peaks")
  returnDataSet = c()

  if (abs(passedData[3, 2]) < 0.001) {
    passedData[, 2] = passedData[, 2] * 1000000
  }

  for (i in (1:(nrow(passedData) / 2))) {
    for (j in (nrow(passedData) / 2):nrow(passedData)) {
      if (i != j) {
        if (passedData[i, 1] == passedData[j, 1]) {
          newRow = c(passedData[i, 1], passedData[i, 2], passedData[j, 2])
          returnDataSet = rbind(newRow, returnDataSet)
        }
      }
    }
  }
  print("Separated forward and reverse scans")

  plot(x = returnDataSet[, 1], y = returnDataSet[, 2], col = "blue", type = "l",
      ylim = c(min(returnDataSet[, 3]), max(returnDataSet[, 2])),
      main = "Separated Plots",
      xlab = allXAxisLabel,
      ylab = allYAxisLabel,
      frame.plot = FALSE)
  lines(x = returnDataSet[, 1], y = returnDataSet[, 3], col = "red")
  legend("topleft", legend = c("Forward Scan", "Reverse Scan"),
        col = c("blue", "red"),
        lty = c(1, 1),
        bty = "n")

  toContinue = FALSE

  while (toContinue == FALSE) {
    if (readline(prompt = "Continue? (y/n) ") == "y") {
      toContinue = TRUE
    }
  }
  return(returnDataSet)
}

selectCropRegion = function(passedData) {
  plot(x = passedData[, 1], y = passedData[, 2], type = "l", xlab = allXAxisLabel, ylab = allYAxisLabel, frame.plot = FALSE)
  toContinue = FALSE

  while (toContinue == FALSE) {
    returnData = c()
    tempData = c()
    upperLimit = max(passedData[, 1])
    lowerLimit = min(passedData[, 1])

    upperLimit = as.numeric(readline("Set upper limit: "))

    for (i in (1:length(passedData[, 1]))) {
      if (passedData[i, 1] <= upperLimit) {
        tempData = rbind(passedData[i,], tempData)
      }
    }
    plot(x = tempData[, 1], y = tempData[, 2], type = "l", xlab = allXAxisLabel, ylab = allYAxisLabel, frame.plot = FALSE)

    lowerLimit = as.numeric(readline("Set lower limit: "))

    for (i in (1:length(tempData[, 1]))) {
      if (tempData[i, 1] >= lowerLimit) {
        returnData = rbind(tempData[i,], returnData)
      }
    }
    plot(x = returnData[, 1], y = returnData[, 2], type = "l", xlab = allXAxisLabel, ylab = allYAxisLabel, frame.plot = FALSE)

    if (readline(prompt = "Change Limits? (y/n) ") == "n") {
      toContinue = TRUE
    }

  }
  return(returnData)
}

fitPolynomial = function(passedData) {
  baselineSubtraction = passedData

  toContinue = FALSE

  while (toContinue == FALSE) {
    fittingPolynomialData = c()
    upperLimit = as.numeric(readline("Set Upper Signal Limit: "))
    lowerLimit = as.numeric(readline("Set Lower Signal Limit: "))

    for (i in (1:length(passedData[, 1]))) {
      if (passedData[i, 1] >= upperLimit) {
        fittingPolynomialData = rbind(passedData[i,], fittingPolynomialData)
      }
      else if (passedData[i, 1] <= lowerLimit) {
        fittingPolynomialData = rbind(passedData[i,], fittingPolynomialData)
      }
    }

    plot(passedData, type = "l", xlab = allXAxisLabel, ylab = allYAxisLabel, frame.plot = FALSE)
    lines(fittingPolynomialData, col = "red")

    if (readline("Change Limits? (y/n) ") == "n") {
      toContinue = TRUE
    }
  }
  n = 2
  toContinue = FALSE
  while (toContinue == FALSE) {
    testVolts = fittingPolynomialData[, 1]
    polyOrder = as.numeric(readline("Enter Polynomial Order: "))
    polyFunction = lm(fittingPolynomialData[, 2] ~ poly(testVolts, polyOrder, raw = TRUE))
    testVolts = passedData[, 1]
    lines(x = testVolts, y = predict(polyFunction, data.frame(x = testVolts)), lty = n)
    n = n + 1

    if (readline("Try another polynomial? (y/n) ") == "n") {
      baselineSubtraction[, 2] = (baselineSubtraction[, 2] - predict(polyFunction, data.frame(x = testVolts)))
      if (mean(fittingPolynomialData[, 2]) > 0) {
        plot(passedData, type = "l", main = "Baseline Subtraction", ylim = c(min(baselineSubtraction), max(passedData)), xlab = allXAxisLabel, ylab = allYAxisLabel, frame.plot = FALSE)
        legend(x = "right", y = 0.4, legend = c("RAW Data", "Polynomial Fitting Region", "Polynomial Baseline", "Baseline Subtraction"), lty = c(1, 1, 1, 1), col = c("black", "red", "blue", "purple"), bty = "n")
      }
      else if (mean(fittingPolynomialData[, 2]) < 0) {
        plot(passedData, type = "l", main = "Baseline Subtraction", ylim = c(min(passedData), 0.05), xlab = allXAxisLabel, ylab = allYAxisLabel, frame.plot = FALSE)
        legend("bottomright", legend = c("RAW Data", "Polynomial Fitting Region", "Polynomial Baseline", "Baseline Subtraction"), lty = c(1, 1, 1, 1), col = c("black", "red", "blue", "purple"), bty = "n")
      }
      lines(fittingPolynomialData, col = "red")
      lines(x = testVolts, y = predict(polyFunction, data.frame(x = testVolts)), lty = 1, col = "blue")
      lines(smooth.spline(baselineSubtraction), col = "purple")

      if (readline(prompt = "Continue? (y/n) ") == "y") {
        toContinue = TRUE
        message("Polynomial Baseline Function:")
        print(summary(polyFunction))
      }
    }
  }
  return(baselineSubtraction)
}

peakWidthAtHalfHeight = function(passedData, peaks) {
  # Column 2 and 4 of peaks contain the Anodic and Cathodic peaks respectively
  # Divide each of these values by 2, then search passedData$y for == value, noting voltage
  # Find difference in voltages and report (peak width at half height)
  points = c()

  for (i in (1:length(peaks[, 2]))) {
    for (j in (1:length(passedData$y))) {
      if (round(peaks[i, 2] / 2, digits = 2) == round(passedData$y[j], digits = 2)) {
        points = c(points, passedData$x[j])
      }
    }
    print(paste("Approximate peak width at half height:", (max(points) - min(points))))
    print(paste("Voltages found:", points))
    points = c()
  }
}

findPeaks = function(passedData) {
  current = passedData$y
  voltage = passedData$x

  if (mean(current) > 0) {
    peakDetection = 0.05
    peakVolts = c(voltage[which(diff(sign(diff(current))) == -2) + 1])
    peakAmps = c(current[which(diff(sign(diff(current))) == -2) + 1])
  }

  else if (mean(current) < 0) {
    peakDetection = -0.05
    peakVolts = c(voltage[which(diff(sign(diff(current))) == +2) + 1])
    peakAmps = c(current[which(diff(sign(diff(current))) == +2) + 1])
  }

  toContinue = FALSE

  while (toContinue == FALSE) {
    listOfPeaks = c()
    for (i in (1:length(peakAmps))) {
      if (abs(peakAmps[i]) > abs(peakDetection)) {
        newRow = c(peakVolts[i], peakAmps[i])
        listOfPeaks = rbind(newRow, listOfPeaks)
      }
    }
    print(paste("Number of peaks size >", peakDetection, ":", nrow(listOfPeaks)))
    colnames(listOfPeaks) = c("Voltage", "Current")
    print(listOfPeaks)
    if (readline("Continue? (y/n) ") == "y") {
      toContinue = TRUE
    } else {
      peakDetection = as.numeric(readline("Enter new peak detection limit: "))
    }
  }
  return(listOfPeaks)
}

redoxPairs = function(dataCV) {
  # Separates Anodic and Cathodic peaks using function
  dataCV = separatePeaks(dataCV)

  anodicData = cbind(dataCV[, 1], dataCV[, 2]) # Change passed to anodic
  colnames(anodicData) = c("Voltage", "Current")
  #cathodicData = cbind(dataCV[,1], dataCV[,3])
  #colnames(cathodicData) = c("Voltage", "Current")

  anodicData = selectCropRegion(anodicData)
  anodicBaseline = fitPolynomial(anodicData)
  anodicBaseline = smooth.spline(anodicBaseline)
  anodicPeaks = findPeaks(anodicBaseline)

  #cathodicData = selectCropRegion(cathodicData)
  #cathodicBaseline = fitPolynomial(cathodicData)
  #cathodicBaseline = smooth.spline(cathodicBaseline)
  #cathodicPeaks = findPeaks(cathodicBaseline)

  redoxPeaks = cbind(anodicPeaks[, 1], anodicPeaks[, 2]) #, cathodicPeaks[,1], cathodicPeaks[,2])
  colnames(redoxPeaks) = c("Anodic Voltage", "Anodic Current") #, "Cathodic Voltage", "Cathodic Current")
  print("Redox Peaks: ")
  print(redoxPeaks)

  plot(x = dataCV[, 1], y = dataCV[, 2], type = "l", col = "blue", ylim = c(min(dataCV[, 3]), max(dataCV[, 2])), main = "Peak Extraction", xlab = allXAxisLabel, ylab = allYAxisLabel, frame.plot = FALSE)
  lines(anodicBaseline, col = "blue", lwd = 2)
  #lines(x = dataCV[,1], y = dataCV[,3], col = "red")
  #lines(cathodicBaseline, col = "red", lwd = 2)
  #legend("topleft", legend = c("Forward Sweep", "Reverse Sweep", "Oxidative Baseline Extraction", "Reductive Baseline Extraction"), lty = c(1,1,1,1), lwd = c(1,1,2,2), col = c("blue", "red", "blue", "red"), bty = "n")

  toContinue = FALSE

  # while(toContinue == FALSE)
  #{
  #line = readline(prompt = "Continue? (y/n) ")

  #if(line == "y"){
  toContinue = TRUE
  #}
  #}
  fitsAnodic = modelPeaks(anodicBaseline, redoxPeaks)
  #fitsCathodic = modelPeaks(cathodicBaseline, redoxPeaks)

  plot(x = dataCV[, 1], y = dataCV[, 2], type = "l", col = "black", ylim = c(min(dataCV[, 3]), max(dataCV[, 2])), main = "Peak Extraction with Modelling", xlab = allXAxisLabel, ylab = allYAxisLabel, frame.plot = FALSE)
  lines(anodicBaseline, col = "blue", lwd = 2)
  lines(x = dataCV[, 1], y = dataCV[, 3], col = "black")
  lines(cathodicBaseline, col = "red", lwd = 2)
  generateCurrent(fitsAnodic[[1]], fitsAnodic[[2]], fitsAnodic[[3]], fitsAnodic[[4]], FALSE, TRUE, "grey", "green", FALSE)
  generateCurrent(fitsCathodic[[1]], fitsCathodic[[2]], fitsCathodic[[3]], fitsCathodic[[4]], FALSE, TRUE, "grey", "green", FALSE)
  legend("topleft", legend = c("Input dcV", "Oxidative Baseline Extraction", "Reductive Baseline Extraction", "Modelled Peaks", "Combined Modelled Peaks"), lty = c(1, 1, 1, 1, 1), lwd = c(1, 1, 1, 1, 1), col = c("black", "blue", "red", "grey", "green"), bty = "n")

  plot(anodicBaseline, type = "l", col = "blue", ylim = c(1.2 * min(cathodicBaseline$y), 1.2 * max(anodicBaseline$y)), main = "Peak Extraction with Modelling", xlab = allXAxisLabel, ylab = allYAxisLabel, frame.plot = FALSE)
  lines(cathodicBaseline, col = "red", lwd = 2)
  generateCurrent(fitsAnodic[[1]], fitsAnodic[[2]], fitsAnodic[[3]], fitsAnodic[[4]], FALSE, TRUE, "grey", "green", FALSE)
  generateCurrent(fitsCathodic[[1]], fitsCathodic[[2]], fitsCathodic[[3]], fitsCathodic[[4]], FALSE, TRUE, "grey", "green", FALSE)
  legend("topright", legend = c("Oxidative Baseline Extraction", "Reductive Baseline Extraction", "Modelled Peaks", "Combined Modelled Peaks"), lty = c(1, 1, 1, 1), lwd = c(1, 1, 1, 1), col = c("blue", "red", "grey", "green"), bty = "n")
}

generateCurrent = function(passedData, peaks, parameters, electrodeCoverage, toPlot, toLine, colPeaks, colCombined, optRun) {
  gasConstant = parameters[1]
  faradConstant = parameters[2]
  scanRate = parameters[3]
  elecArea = parameters[4]
  temp = parameters[5]
  nElec = parameters[6]
  nAppElec = parameters[7]
  avgRedPot = parameters[8]
  elecCoverage = c(electrodeCoverage)
  R2 = c()
  summativeCurrent = c()
  avgError = 0

  if (toPlot == TRUE) {
    plot(smooth.spline(passedData), type = "l", main = "Modelled Peaks", xlab = allXAxisLabel, ylab = allYAxisLabel, ylim = c(min(passedData$y), (1.2 * max(passedData$y))), frame.plot = FALSE)
  }

  for (i in (1:length(peaks[, 1])))
    # Models current and calculates error between RAW and Model
  {
    covXArea = ((elecCoverage[i] * (10 ^ -12)) * elecArea)
    if (optRun == FALSE) {
      oxPeakPot = peaks[i, 1]
      #redPeakPot = peaks[i,3]

      if (mean(passedData$y) > 0) {
        collectionOfConstants = (nElec * nAppElec * (faradConstant ^ 2) * scanRate * covXArea) / (gasConstant * temp)
        avgRedPot = oxPeakPot
      }
      else {
        collectionOfConstants = -(nElec * nAppElec * (faradConstant ^ 2) * scanRate * covXArea) / (gasConstant * temp)
        avgRedPot = redPeakPot
      }
    }
    else {
      avgRedPot = peaks[i, 1]
      nAppElec = peaks[i, 2]
      if (mean(passedData$y) > 0) {
        collectionOfConstants = (nElec * nAppElec * (faradConstant ^ 2) * scanRate * covXArea) / (gasConstant * temp)
      }
      else {
        collectionOfConstants = -(nElec * nAppElec * (faradConstant ^ 2) * scanRate * covXArea) / (gasConstant * temp)
      }
    }
    #print(paste("Redox Potential for Peak",i,":", avgRedPot))
    modelledCurrent = c()
    for (j in (1:length(passedData$x)))
      # Predicts current
    {
      topHalfOfFraction = exp((nAppElec * faradConstant * (passedData$x[j] - avgRedPot)) / (gasConstant * temp))
      bottomHalfOfFraction = (1 + topHalfOfFraction) ^ 2

      predictedCurrent = collectionOfConstants * (topHalfOfFraction / (bottomHalfOfFraction) * 1000000)
      modelledCurrent = rbind(modelledCurrent, c(passedData$x[j], predictedCurrent))
    }
    if (toLine == TRUE) {
      lines(x = modelledCurrent[, 1], y = modelledCurrent[, 2], col = colPeaks, lty = 1)
    }
    if (i == 1) {
      summativeCurrent = modelledCurrent
    }
    else {
      summativeCurrent[, 2] = (summativeCurrent[, 2] + modelledCurrent[, 2])
      if (toLine == TRUE) {
        lines(summativeCurrent, col = colCombined, lwd = 2)
      }
    }
  }

  if (toPlot == TRUE) {
    legend("right", legend = c("Baseline Subtraction", "Modelled Peaks", "Combined Model Peaks"), lty = c(1, 1, 1), col = c("black", colPeaks, colCombined), bty = "n")
  }
  for (i in 1:length(summativeCurrent[, 2])) {
    R2[i] = (passedData$y[i] - summativeCurrent[i, 2]) ^ 2
  }
  return(R2)
}
# Generates current, plots modelled peaks and returns avg fit error

optimiseCoverage = function(electrodeCoverage, passedData, parameters, peaks) {
  gasConstant = parameters[1]
  faradConstant = parameters[2]
  scanRate = parameters[3]
  elecArea = parameters[4]
  temp = parameters[5]
  nElec = parameters[6]
  nAppElec = parameters[7]
  nOptElec = c()
  avgRedPot = parameters[8]
  startR2 = parameters[9]
  storedSumR2 = sum(startR2)
  elecCoverage = electrodeCoverage
  varsToOpt = c()
  itDataFrame = c()

  for (i in 1:length(peaks[, 1])) {
    oxPeakPot = peaks[i, 1]
    #redPeakPot = peaks[i,3]
    nOptElec[i] = nAppElec

    if (mean(passedData$y) > 0) {
      avgRedPot[i] = oxPeakPot
    }
    else {
      avgRedPot[i] = redPeakPot
    }
  }

  for (i in 1:length(avgRedPot)) {
    varsToOpt = rbind(varsToOpt, c(avgRedPot[i], nOptElec[i], elecCoverage[i]))
  }

  modVector = sample.int(3, length(varsToOpt[1,]), replace = TRUE)

  for (i in 1:50000) {
    optExport = c()
    passRedPot = c()
    passElec = c()
    passElecCoverage = c()
    optimisingVars = c()

    for (j in 1:length(varsToOpt[, 1])) {
      modVector = modVector - 2
      optimisingVars = varsToOpt + (varsToOpt * modVector) / 100


      passRedPot = c(passRedPot, optimisingVars[j, 1])
      passElec = c(passElec, optimisingVars[j, 2])
      passElecCoverage = c(passElecCoverage, optimisingVars[j, 3])
    }

    optExport = cbind(passRedPot, passElec)
    exportParameters = c(gasConstant, faradConstant, scanRate, elecArea, temp, nElec, passElec, passRedPot)

    if (i %% 100 == 0) {
      R2 = generateCurrent(passedData, optExport, exportParameters, passElecCoverage, TRUE, TRUE, "red", "blue", TRUE)
    }
    else {
      R2 = generateCurrent(passedData, optExport, exportParameters, passElecCoverage, FALSE, FALSE, "red", "blue", TRUE)
    }

    itDataFrame = rbind(itDataFrame, c(i, sum(R2)))

    if (sum(R2) < storedSumR2) {
      message("######################################################################")
      storedSumR2 = sum(R2)
      print(varsToOpt)
      varsToOpt = optimisingVars
      print(varsToOpt)
      modVector = modVector * 3
    }
    else {
      optimisingVars = varsToOpt
      modVector = sample.int(3, length(varsToOpt[1,]), replace = TRUE)
    }
  }
  plot(x = itDataFrame[, 1], y = itDataFrame[, 2], xlab = "Iteration", ylab = "Sum of R2")
  return(varsToOpt)
}

modelPeaks = function(passedData, peaks) {
  gasConstant = 8.3144621
  faradConstant = 96485.3365
  scanRate = 0.06 # as.numeric(readline("Enter Scan Rate (Vs-1): "))
  elecArea = 0.07 # as.numeric(readline("Enter Electrode Area (cm2): "))
  temp = 278 # as.numeric(readline("Enter Temperature (K): "))
  nElec = 1 # as.numeric(readline("Enter n(s): "))
  nAppElec = 1 # as.numeric(readline("Enter n(app): "))
  avgRedPot = 0
  elecCoverage = c(650, 150)
  exportParameters = c(gasConstant, faradConstant, scanRate, elecArea, temp, nElec, nAppElec, avgRedPot)
  finalPlotList = list()

  toContinue = FALSE
  while (toContinue == FALSE) {
    plot(smooth.spline(passedData), type = "l", main = "Modelled Peaks", xlab = allXAxisLabel, ylab = allYAxisLabel, frame.plot = FALSE)
    for (i in 1:length(peaks[, 1])) {
      # elecCoverage[i] = as.numeric(readline("Enter starting electrode coverage - less than baseline subtraction (pmol cm-2): "))
      avgError = generateCurrent(passedData, peaks, exportParameters, elecCoverage, TRUE, TRUE, "red", "blue", FALSE)
    }

    if (readline("Try more parameters? (y/n): ") == "n") {
      if (readline("Optimise Coverage? (y/n): ") == "n") {
        toContinue = TRUE
        message("Electrode coverage:")
        for (i in 1:length(peaks[, 1])) {
          print(paste("Peak ", i, ": ", elecCoverage[i], sep = ""))
        }
      }
      else {
        exportParameters = c(gasConstant, faradConstant, scanRate, elecArea, temp, nElec, nAppElec, avgRedPot, sum(avgError))
        optCoverage = optimiseCoverage(elecCoverage, passedData, exportParameters, peaks)

        print(optCoverage)
        for (i in 1:length(peaks[, 1])) {
          print(paste("Peak ", i, ": ", elecCoverage[i], sep = ""))
        }

        if (readline("Use optimised coverages? (y/n): ") == "y") {
          elecCoverage = optCoverage
          message("Electrode coverage:")
          for (i in 1:length(peaks[, 1])) {
            print(paste("Peak ", i, ": ", elecCoverage[i], sep = ""))
          }
          toContinue = TRUE
        }

        else {
          toContinue = TRUE
          message("Electrode coverage:")
          for (i in 1:length(peaks[, 1])) {
            print(paste("Peak ", i, ": ", elecCoverage[i], sep = ""))
          }
        }
      }
    }
  }
  generateCurrent(passedData, peaks, exportParameters, elecCoverage, TRUE, TRUE, "red", "blue", FALSE)
  #  peakWidthAtHalfHeight(passedData, peaks[,1:2])
  #  peakWidthAtHalfHeight(passedData, peaks[,3:4])

  toContinue = FALSE
  while (toContinue == FALSE) {
    if (readline("Continue? (y/n) ") == "y") {
      toContinue = TRUE
      finalPlotList = list(passedData, peaks, exportParameters, elecCoverage)
      return(finalPlotList)
    }
  }
}

# Read in initial CV file and plots
inputCV = read.csv("Test dcV.csv")

plot(x = inputCV[, 1], y = inputCV[, 2], main = "Loaded dcV", type = "l", xlab = allXAxisLabel, ylab = allYAxisLabel, frame.plot = FALSE)

redoxPairs(inputCV)
