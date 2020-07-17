import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import SpanSelector

allXAxisLabel = "Voltage / V"
allYAxisLabel = "Current / \u03BCA"


def continueQuestion(contCriteria):
    toContinue = False
    while(toContinue == False):
        if (input("Continue? (y/n) ") == "y"):
            toContinue = True
            return(contCriteria)


def separatePeaks(inputdcV):
    print("Separating forward and reverse scans...")
    if(inputdcV[2][1] < 0.001):
        inputdcV[:, 1] = inputdcV[:, 1]*1000000

    nDataPoints = len(inputdcV)
    midPoint = (nDataPoints//2)
    allData = inputdcV[:midPoint]
    reverseScan = inputdcV[midPoint:]
    allData = np.append(allData, reverseScan, axis=1)
    print("Separated forward and reverse scans")

    plt.plot(allData[:, 0], allData[:, 1], 'r', label="Forward Scan")
    plt.plot(allData[:, 2], allData[:, 3], 'b', label="Reverse Scan")
    plt.xlabel(allXAxisLabel)
    plt.ylabel(allYAxisLabel)
    plt.legend(loc='upper left')
    plt.show()
    return(continueQuestion(allData))


def selectCroppingRegion(xmin, xmax):
    indmin, indmax = np.searchsorted(x, (xmin, xmax))

    thisx = x[indmin: indmax]
    thisy = y[indmin: indmax]
    zoomedLine.set_data(thisx, thisy)
    ax2.set_xlim(thisx.min(), thisx.max())
    ax2.set_ylim(thisy.min(), thisy.max())
    fig.canvas.draw_idle()

    np.savetxt("text.out", np.c_[thisx, thisy])


# loads in input dcV file and separates into anodic and cathodic data
inputdcV = np.genfromtxt('Test dcV.csv', delimiter=',')
allData = separatePeaks(inputdcV)
anodicData = allData[:, 0: 2]
cathodicData = allData[:, 2:]

# selects anodic window and temporarily saves it in text.out file
print("Selecting anodic activity region")
fig = plt.figure()
ax = fig.add_subplot(211)
x = anodicData[:, 0]
y = anodicData[:, 1]
ax.plot(x, y, '-')
ax.set_title(
    'Select region of interest by left clicking and dragging \n Close window when done')
ax2 = fig.add_subplot(212)
zoomedLine, = ax2.plot(x, y, '-')
interestRegion = SpanSelector(ax, selectCroppingRegion, 'horizontal', useblit=True,
                              rectprops=dict(alpha=0.5, facecolor='red'), span_stays=True)
plt.show()
anodicWindow = np.genfromtxt('text.out', delimiter=" ")

# selects cathodic window and temporarily saves it in text.out file
print("Selecting cathodic activity region")
fig = plt.figure()
ax = fig.add_subplot(211)
x = cathodicData[:, 0]
y = cathodicData[:, 1]
x = x[::-1]
y = y[::-1]
ax.plot(x, y, '-')
ax.set_title(
    'Select region of interest by left clicking and dragging \n Close window when done')
ax2 = fig.add_subplot(212)
zoomedLine, = ax2.plot(x, y, '-')
interestRegion = SpanSelector(ax, selectCroppingRegion, 'horizontal', useblit=True,
                              rectprops=dict(alpha=0.5, facecolor='red'), span_stays=True)
plt.show()
cathodicWindow = np.genfromtxt('text.out', delimiter=" ")

plt.plot(anodicData[:, 0], anodicData[:, 1], 'b--')
plt.plot(cathodicData[:, 0], cathodicData[:, 1], 'r--')
plt.plot(anodicWindow[:, 0], anodicWindow[:, 1], 'b-')
plt.plot(cathodicWindow[:, 0], cathodicWindow[:, 1], 'r-')
plt.show()
