import numpy as np
import matplotlib.pyplot as plt

nDim = 3
sigma = 3e-10
epsilon = 1



# returns force from LJ potential with rVec: target -> interacting particle
def singleInteractionForce(rVec):           # rVec is the seperation from target particle to surrounding particles
    r = np.sqrt(np.dot(rVec, rVec))         # norm of rVec
    F = epsilon*(48*(sigma**12)*(r**-14) - 24*(sigma**6)*(r**-8))*rVec
    return F
    # single interaction force. For total, sum up all interactions

# takes a list of vectors from target to interacting particles, returns total force
def totalInteractionForce(rs):
    totalForce = np.zeros(nDim)
    for i in range(len(rs)):
        totalForce += singleInteractionForce(rs[i])
    return totalForce


# Finds nearest copy for minimum interaction convention
def MICneighbour(rVec):
    return np.mod(rVec + 0.5*boxDimensions, boxDimensions) - 0.5*boxDimensions  # boxDimenions is the x,y,z vector of the size of the box. E.g. if its a square box of sizeN then boxDimensions = L*np.ones(nDim)
    # returns the conversion to nearest periodic clone






