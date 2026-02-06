import numpy
import matplotlib.pyplot as plt

sigma = 3e-10
epsilon = 1

def lennardJones(r):

# returns force from LJ potential with rVec: target -> interacting particle
def singleInteractionForce(rVec):           # rVec is the seperation from target particle to surrounding particles
    r = np.sqrt(np.dot(rVec, rVec))         # norm of rVec
    F = epsilon*(48*(sigma**12)*(r**-14) - 24*(sigma**6)*(r**-8))*rVec
    return F
    # single interaction force. For total, sum up all interactions

def 


