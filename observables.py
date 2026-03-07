import numpy as np
from forces import pairwiseForce
from forces import rMIC

import matplotlib.pyplot as plt


"""
Calculates an array of all pairwise distances r_ij [nParticles*(nParticles-1)]
pos is an array of all positions within the box [nParticles, nDimensions]
boxDimensions is the x,y,z size array of the box. [nDimensions]
Can be used to get the correlation function by plotting into a histogram
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
    return pairwiseDistances




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
    


