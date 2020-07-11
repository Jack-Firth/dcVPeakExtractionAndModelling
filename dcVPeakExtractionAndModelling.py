import numpy as np
import matplotlib.pyplot as plt

allXAxisLabel = "Voltage / V"
allYAxisLabel = "Current / \u03BCA"


def separatePeaks(inputdcV):
    print("Separating forward and reverse scans...")
    if(inputdcV[2][1] < 0.001):
        inputdcV[:, 1] = inputdcV[:, 1]*1000000

    nDataPoints = len(inputdcV)
    midPoint = nDataPoints//2
    forwardScan = inputdcV[:midPoint]
    reverseScan = inputdcV[midPoint:]
    print("Separated forward and reverse scans")

    plt.plot(forwardScan[:, 0], forwardScan[:, 1], 'r', label="Forward Scan")
    plt.plot(reverseScan[:, 0], reverseScan[:, 1], 'b', label="Reverse Scan")
    plt.xlabel(allXAxisLabel)
    plt.ylabel(allYAxisLabel)
    plt.legend(loc='upper left')

    plt.show()

    toContinue = False
    while(toContinue == False):
        if (input("Continue? (y/n) ") == "y"):
            toContinue = True
    return(forwardScan, reverseScan)


def selectCroppingRegion(inputdcV):
    print("Selecting activity region")
    plt.plot(inputdcV[:, 0], inputdcV[:, 1])
    plt.show()


inputdcV = np.genfromtxt('Test dcV.csv', delimiter=',')
allData = separatePeaks(inputdcV)
anodicData = allData[0]
cathodicData = allData[1]

selectCroppingRegion(anodicData)
