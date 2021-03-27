import sys
import math
import numpy as np
import scipy
from scipy.signal import savgol_filter
import sympy
import matplotlib.pyplot as plt
from matplotlib.widgets import SpanSelector

allXAxisLabel = "Voltage / V"
allYAxisLabel = "Current / \u03BCA"

def runLog(message, newRun=False): # Saves any posts to terminal into log.out file instead
    if(newRun == True):
        f = open(".out/log.out", 'w')
    else:
        f = open(".out/log.out", 'a')
    f.write(str(message) + '\n')
    f.close()
    

def separatePeaks(inputdcV):   # Separates peaks by splitting dataset into 2 sections at the half-way point
    runLog("Separating forward and reverse scans...")
    if(inputdcV[2][1] < 0.001):
        inputdcV[:, 1] = inputdcV[:, 1] * 1000000
    nDataPoints = len(inputdcV)
    if(nDataPoints % 2 == 1):
        inputdcV = inputdcV[0:(nDataPoints-1)]
    midPoint = (nDataPoints//2)
    forwardScan = inputdcV[:midPoint]
    reverseScan = inputdcV[midPoint:]
    allData = np.column_stack((forwardScan, reverseScan))
    runLog("Separated forward and reverse scans")

    #plt.plot(allData[:, 0], allData[:, 1], 'r', label="Forward Scan")
    #plt.plot(allData[:, 2], allData[:, 3], 'b', label="Reverse Scan")
    #plt.xlabel(allXAxisLabel)
    #plt.ylabel(allYAxisLabel)
    #plt.legend(loc='upper left')
    #plt.show(block=False)

    # toContinue = True
    # while(toContinue == False):
    #    if(input("Continue? (y/n): ") == 'y'):
    #        toContinue = True
    return(allData)


def fitPolynomialBaseline(fullDataset, echemRegion):    # Removes Echem region from raw data to then fit polynomials
    np.savetxt(".out/fulldataset.out", np.c_[fullDataset])
    np.savetxt(".out/echem.out", np.c_[echemRegion])
    fittingData = np.zeros((len(fullDataset), 2))

    eWindowStart = echemRegion[0][0]
    eWindowEnd = echemRegion[len(echemRegion)-1][0]
    runLog("Removing region " + eWindowStart.astype('str') +
          " - " + eWindowEnd.astype('str') + " V")

    for i in range(len(fullDataset)):
        if(fullDataset[i, 0].astype('str') <= eWindowStart.astype('str')):
            fittingData[i, :] = fullDataset[i, :]
        elif(fullDataset[i, 0].astype('str') >= eWindowEnd.astype('str')):
            fittingData[i, :] = fullDataset[i, :]

    fittingData = np.ma.masked_equal(fittingData, 0)
    start = int(round((0.15*len(fittingData)), 1))
    stop = int(round((0.85*len(fittingData)), 1))
    slice_set = slice(start, stop)
    fittingData = fittingData[slice_set]
    
    np.savetxt(".out/fittingData.out", np.c_[fittingData])
    return(fittingData)


def plotPolys(plt, fullDataset, fittingData, echemRegion, polyQ=False):     # Generates polynomials and plots them to a given graph
    polyArray = np.zeros(len(fullDataset[:, 0]))
    for i in range(3, 11):
        coefficients = np.polyfit(
            x=fittingData[:, 0], y=fittingData[:, 1], deg=i)
        poly = np.poly1d(coefficients)
        baseline_y = poly(fullDataset[:, 0])
        plt.plot(fullDataset[:, 0], baseline_y,
                 alpha=0.5, label=i, linestyle="--")
        polyArray = np.c_[polyArray, baseline_y]
        plt.plot(echemRegion[:, 0], echemRegion[:, 1], 'r-')
        plt.legend()
        plt.set_xlim((0.95*min(echemRegion[:,0])), (1.05*max(echemRegion[:,0])))
        plt.set_ylim((0.95*min(echemRegion[:,1])), (1.05*max(echemRegion[:,1])))

    if(polyQ == True):
        polyArray = np.delete(polyArray, 0, 1)
        np.savetxt(".out/polys.out", polyArray)


def selectCroppingRegion(xmin, xmax):   # Selects region of interest and calls plotPolys to plot polynomial baselines
    global fig
    global ax
    global ax2
    global voltageUISelect
    global currentUISelect
    global zoomedLine

    try:

        indmin, indmax = np.searchsorted(voltageUISelect, (xmin, xmax))

        thisx = voltageUISelect[indmin: indmax]
        thisy = currentUISelect[indmin: indmax]
        zoomedLine.set_data(thisx, thisy)
        ax2.set_xlim(thisx.min()-0.1, thisx.max()+0.1)
        ax2.set_ylim(thisy.min()-1, thisy.max()+1)
        ax2.cla()

        echemData = np.c_[voltageUISelect, currentUISelect]
        echemWindow = np.c_[thisx, thisy]

        plotPolys(ax2, echemData, fittingData=fitPolynomialBaseline(
            echemData, echemWindow), echemRegion=echemWindow, polyQ=True)
        fig.canvas.draw_idle()

        baseline = np.genfromtxt(".out/polys.out")
        np.savetxt(".out/window.out", np.c_[thisx, thisy])
        np.savetxt(".out/baseline.out", baseline)
        runLog("updated window file")
    except Exception as e:
        runLog(e)


def baselineSubtraction(echemData, baseline):   # Subtracts one data set from another
    subtractedData = np.zeros((len(echemData[:, 0]), 2))
    for i in range(0, len(echemData)):
        subtractedData[i, 1] = echemData[i, 1] - baseline[i, 1]
        subtractedData[i, 0] = echemData[i, 0]
    return(subtractedData)


def focusData(baselineSub, echemWindow):    # Shortens baseline subtracted dataset to echem region window
    baselineReturnData = np.zeros((len(echemWindow), 2))
    echemWindowReturnData = np.zeros((len(echemWindow), 2))
    for i in range(0,len(echemWindow)):
        for j in range(0,len(baselineSub)):
            if(baselineSub[j,0] == echemWindow[i,0]):
                baselineReturnData[i,:] = baselineSub[j,:]
                echemWindowReturnData[i, :] = echemWindow[i, :]
    return(baselineReturnData, echemWindowReturnData)


def findPeaks(dataset, echemData, baseline):    # Identifies peaks of interest based on prominence
    if((sum(dataset[:,1])/len(dataset[:,1])) < 0):
        dataset[:,1] = dataset[:,1] * -1
        echemData[:,1] = echemData[:,1] * -1
        baseline[:,1] = baseline[:,1] * -1
    peaks, peakInfo = scipy.signal.find_peaks(dataset[:, 1], height=1, prominence=0.01, width=1)
    returnPeaks = np.zeros((len(peaks), 2))
    returnPeaks[:,0] = dataset[peaks, 0]
    returnPeaks[:,1] = dataset[peaks, 1]
    plt.plot(dataset[peaks, 0], dataset[peaks, 1], "x", markersize="10", label="Peaks")
    plt.plot(dataset[:, 0], dataset[:, 1], label="Subtraction")
    plt.plot(echemData[:,0], echemData[:,1], color='r', label="Raw Data")
    plt.xlim(0.9*min(dataset[:,0]), 1.1*max(dataset[:,0]))
    plt.legend()
    plt.show(block=False)
    
    return(returnPeaks)


def genCurrent(fParams, inputData):
    returnCurrent = np.zeros((len(inputData), 2))
    constants = (fParams["nElec"]*fParams["nAppElec"]*(fParams["faradConstant"]**2)*fParams["scanRate"]*fParams["elecArea"]*fParams["elecCoverage"])/(fParams["gasConstant"]*fParams["temp"])
    for i in range(0,len(inputData[:,0])):
        deltaE = inputData[i,0] - fParams["potential"]
        topHalfOfFraction = math.e**(fParams["nAppElec"] * fParams["faradConstant"] * deltaE/ (fParams["gasConstant"] * fParams["temp"]))
        bottomHalfOfFraction = (1 + topHalfOfFraction)**2
        modelledCurrent = constants*(topHalfOfFraction/bottomHalfOfFraction) * 1000000
        returnCurrent[i,0] = inputData[i,0]
        returnCurrent[i,1] = modelledCurrent
    plt.plot(returnCurrent[:, 0], returnCurrent[:,1], label="Model")
    plt.show()

def selectionUI(echemData, echemWindow=np.empty([]), fittingData=np.empty([])): # Front end of selection
    global fig
    global ax
    global ax2
    global voltageUISelect
    global currentUISelect
    global zoomedLine

    faradaicParameters = {
        "gasConstant" : 8.3144621,
        "faradConstant" : 96485.3365,
        "scanRate" : 0.06, 
        "elecArea" : 0.07, 
        "temp" : 278, 
        "nElec" : 1, 
        "nAppElec" : 1, 
        "elecCoverage" : 650E-12,
        "potential" : 0
    }

    runLog("Selecting activity region")

    fig = plt.figure()
    ax = fig.add_subplot(121)
    plt.ylabel(allYAxisLabel)
    plt.xlabel(allXAxisLabel)
    ax2 = fig.add_subplot(122)
    plt.xlabel(allXAxisLabel)

    voltageUISelect = echemData[:, 0]
    currentUISelect = echemData[:, 1]

    if((currentUISelect[10] - currentUISelect[15]) > 0):
        voltageUISelect = voltageUISelect[::-1]
        currentUISelect = currentUISelect[::-1]

    fig.suptitle(
        'Select region of interest by left clicking and dragging \n Close window when done')
    ax.plot(voltageUISelect, currentUISelect, 'b-')

    zoomedLine, = ax2.plot(voltageUISelect, currentUISelect,
                           'r-', label="Window")
    ax2.plot(voltageUISelect, currentUISelect, 'b--', alpha=0.3)

    toContinue = False
    while(toContinue == False):
        try:
            interestRegion = SpanSelector(ax, selectCroppingRegion, 'horizontal', useblit=True,
                                      rectprops=dict(alpha=0.5, facecolor='red'), span_stays=True)
            plt.legend()
            plt.show(block=False)
        except Exception as e:
            runLog(e)

        echemWindow = np.genfromtxt('.out/window.out', delimiter=" ")
        runLog("read window file")
        nPoly = int(input("Select Polynomial Order: "))
        plt.close()

        polyBaselines = np.genfromtxt('.out/polys.out', delimiter=" ")
        baseline = np.zeros((len(echemData), 2))
        (baseline[:, 0], baseline[:, 1]) = (
            echemData[:, 0], polyBaselines[:, (nPoly-3)])

        subtractedData = baselineSubtraction(echemData, baseline)
        (regionOfInterest, shorterEchem) = focusData(subtractedData, echemWindow)
        identifiedPeaks = findPeaks(regionOfInterest, shorterEchem, baseline)
        runLog("Identified Peaks at:" )
        runLog(identifiedPeaks)

        for i in range(0,len(identifiedPeaks - 1)):
            faradaicParameters["potential"] = identifiedPeaks[i,0]
            genCurrent(faradaicParameters, echemData)

        # if(input("Continue? (y/n): ") == "y"):
        toContinue = True
        return(echemData, echemWindow, baseline)



runLog("Running", True)
# loads in input dcV file and separates into anodic and cathodic data
inputdcV = np.loadtxt('Test dcV.csv', delimiter=',', skiprows=1)
allData = separatePeaks(inputdcV)
anodicData = allData[:, 0: 2]
cathodicData = allData[:, 2:]

(anodicData, anodicWindow, anodicBaseline) = selectionUI(anodicData)
#(cathodicData, cathodicWindow, cathodicBaseline) = selectionUI(cathodicData)

