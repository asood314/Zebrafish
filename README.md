# Zebrafish
# Example of extracting neuron dF/F traces from movie

# setup
import sys
sys.path.append(pathToCalBlitz)
sys.path.append(pathToZebrafishTools)
from ZebrafishTools import ZebrafishTools

# initialize ZebrafishTools
zft = ZebrafishTools([MovieFile1,MovieFile2,...],nchan=2)
# choose movie to analyze and motion correct
zft.motionCorrect(index=movieIndex)
# OR
zft.motionCorrect(filename=movieFile)
# extract ROIs
zft.findNeuronsFromStDev()
# compute neuron positions, sizes, dF/F
zft.computeNeuronVars()
# compare ROIs to standard deviation of movie frames
zft.plotROIsOverlay()
# plot neuron dF/F
import pylab as plt
plt.plot(zft.neuronDFFs[neuronIndex])