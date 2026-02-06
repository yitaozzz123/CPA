import numpy
import matplotlib.pyplot as plt

sigma = 3e-10
epsilon = 1

def lennardJones(r):

def singleInteractionForce(rVec):
    r = np.sqrt(np.dot(rVec, rVec))
    F = epsilon*(48*(sigma**12)*(r**-14) - 24*(sigma**6)*(r**-8))*rVec
    return F



