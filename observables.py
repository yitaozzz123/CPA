import numpy as np
from forces import rMIC

import matplotlib.pyplot as plt


"""
Calculates an array of all pairwise distances r_ij [nParticles*(nParticles-1)]
pos is an array of all positions within the box [nParticles, nDimensions]
boxDimensions is the x,y,z size array of the box. [nDimensions]
"""
def calculatePairwiseDistances(pos, boxDimensions):
    nParticles = len(pos)
    pairwiseDistances=[]
    for i in range(nParticles):
        # 1. take the difference in position between particle i and each other particle
        # 2. remove the zero vector corresponding to self interaction
        # 3. convert seperations into MIC nearest clone seperations.
        deltaPos = rMIC(np.delete(pos[i]-pos, i, 0), boxDimensions)
        deltaPosNorm = np.linalg.norm(deltaPos, axis = 1)
        pairwiseDistances.extend(deltaPosNorm)
    return np.array(pairwiseDistances)




"""
Calculates total pressure (float)
pos is an array of all positions within the box [nParticles, nDimensions]
T is dimensionless temperature (float)
numberDensity is the number density of the particles (float)
boxDimensions is the x,y,z size array of the box. [nDimensions]
"""
def calculatePressure(pos, T, boxDimensions): 
    # number of particles and particle number density
    nParticles = len(pos)
    numberDensity = nParticles/np.prod(boxDimensions)
    # Get all pairwise distances from the calculatePairwiseDistances function
    pairwiseDistances = calculatePairwiseDistances(pos,boxDimensions)
    # sum over all contributions to pressure for each pairwise distance
    sum = 0
    for i in range(len(pairwiseDistances)):
        pairwiseDistance = pairwiseDistances[i]
        sum += -48*pairwiseDistance**-12 + 24*pairwiseDistance**-6
    # Put the sum into the formula for pressure (Verlet)
    P = T*numberDensity*(1-sum/(12*nParticles*T))
    return P
    



"""
Calculates the radial correlation function density g(r) and corresponding distance r ([nBins], [nBins]).
pos is an array of all positions within the box [nParticles, nDimensions]
boxDimensions is the x,y,z size array of the box. [nDimensions]
Only works in 3D!
"""
def calculateCorrelationFunction(pos, boxDimensions, nBins):
    # number of particles and volume
    nParticles = len(pos)
    volume = np.prod(boxDimensions)

    # Maximum valid distance in the histogram is L/2, beyond that g(r) is unphysical
    maxDistance = np.min(boxDimensions)/2

    # get all pairwise distances between particles
    pairwiseDistances = calculatePairwiseDistances(pos,boxDimensions)

    # calculate absolute correlation histogram
    absoluteCorrelation, binEdges = np.histogram(pairwiseDistances, bins=nBins, range=[0,maxDistance])
    # convert bin edges into bin centers
    rBins = (binEdges[:-1] + binEdges[1:])/2
    drBin = rBins[1]-rBins[0]               # distance spacing between bins

    # Convert absolute correlation to the proper radial correlation function
    if len(boxDimensions) == 1:    # 1D
        radialCorrelationDensities = 2*volume*absoluteCorrelation/(nParticles*(nParticles-1)*2*drBin)

    elif len(boxDimensions) == 2:   # 2D
        radialCorrelationDensities = 2*volume*absoluteCorrelation/(nParticles*(nParticles-1)*2*np.pi*rBins*drBin)

    elif len(boxDimensions) == 3:     # 3D
        radialCorrelationDensities = 2*volume*absoluteCorrelation/(nParticles*(nParticles-1)*4*np.pi*rBins**2*drBin)

    else:
        print("Not appropriate dimensionality. 1D, 2D or 3D only.")
    return radialCorrelationDensities, rBins


"""
# JUST SOME TESTING. YOU CAN IGNORE
L = 10
boxDimensions = np.array([L,L,L])
T = 1
xs = np.arange(0,L,1)
pos = []
for i in range(0,L,1):
    for j in range(0,L,1):
        for k in range(0,L,1):
            pos.append([i,j,k])
pos = np.array(pos)
nBins = 50

ys, xs = calculateCorrelationFunction(pos, boxDimensions, nBins)
plt.plot(xs,ys)
plt.show()
print(calculatePressure(pos, T, boxDimensions))
"""

