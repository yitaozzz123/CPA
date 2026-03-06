import numpy as np
from forces import pairwiseForce
from forces import rMIC

"""
Calculates total pressure (float)
pos is an array of all positions within the box [nParticles, nDimensions]
T is dimensionless temperature (float)
numberDensity is the number density of the particles (float)
boxDimensions is the x,y,z size array of the box. [nDimensions]
nDims is the number of dimensions (int)
"""
def calculatePressure(pos, T, boxDimensions, nDims): 
    # number of particles and particle number density
    nParticles = len(pos)
    numberDensity = nParticles/np.prod(boxDimensions)
    # loop through each particle i
    for i in range(nParticles):
        # 1. take the difference in position between particle i and each other particle
        # 2. remove the zero vector corresponding to self interaction
        # 3. convert seperations into MIC nearest clone seperations.
        deltaPos = rMIC(np.delete(pos[i]-pos, i, 0), boxDimensions)
        # calculate net force on particle i via Lennard Jones potential
        sum = 0
        for j in range(len(deltaPos)):
            sum -= np.dot(pairwiseForce(deltaPos[j]),deltaPos[j])
        P = T*numberDensity*(1-sum/(12*nParticles*T))
    return P
