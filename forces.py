

import numpy as np
import matplotlib.pyplot as plt



# physical parameters corresponding to Argon
sigma = 1
epsilon = 1
mass = 1

"""
All functions, quantities are given in natural, units, i.e. dimensionless quantities.
"""


"""
Calculates interaction force between two particles. [nDimensions]
deltaR is the vector from target -> interacting particle [nDimensions]
returns the force experienced from the interacting particle
"""
def pairwiseForce(deltaR):
    # calculate the norm of rVec = r        
    rNorm = np.sqrt(np.dot(deltaR, deltaR))
    # calculate the force via F = -nabla U using the Lennard-Jones potential
    # force = epsilon*(48*(sigma**12)*(rNorm**-14) - 24*(sigma**6)*(rNorm**-8))*r
    force = (48*(rNorm**-14) - 24*(rNorm**-8))*deltaR
    return force



"""
Calculates net force of one particle experienced from all other particles [nDimensions]
deltaPos is an array of vectors from target -> interacting particles [nParticles, nDimensions]
nDims is the number of dimensions. (int)
returns the net force experienced from all the interactions on that particle
"""
def netForce(deltaPos, nDims):
    # loop through each particle and sum up its pairwise force contribution
    force = np.zeros(nDims)
    for i in range(len(deltaPos)):
        force += pairwiseForce(deltaPos[i])
    return force



"""
Converts a seperation vector between particles inside a box, to the seperation in the Minimum Image Convention (MIC) clone [nDimensions]
deltaR is the seperation vector between to particles [nDimensions] e.g. (delta x, delta y, delta z)
boxDimensions is the x,y,z size array of the box. E.g. a cubic box has (L,L,L)
"""
def rMIC(deltaR, boxDimensions):
    # shift the vector r, then take the modulus, which returns the right periodic image.
    # Then re-shift r back to origin.
    return np.mod(deltaPos + 0.5*boxDimensions, boxDimensions) - 0.5*boxDimensions




"""
THIS IS THE FUNCTION YOU USE
Calculates all forces of all particles [nParticles, nDimensions]
pos is an array of all positions within the box [nParticles, nDimensions]
boxDimensions is the x,y,z size array of the box. [nDimensions]
nDims is the number of dimensions. (int)
Returns an array of forces
"""
def calculateForces(pos, boxDimensions, nDims): 
    # loop through each particle i
    # fs is the an array of net-force vectors for each particle 
    fs = np.zeros((len(pos),nDims))
    for i in range(len(pos)):
        # 1. take the difference in position between particle i and each other particle
        # 2. remove the zero vector corresponding to self interaction
        # 3. convert seperations into MIC nearest clone seperations.
        deltaPos = rMIC(np.delete(pos[i]-pos, i, 0), boxDimensions)
        # calculate net force on particle i via Lennard Jones potential
        fs[i] = netForce(deltaPos, nDims)
    return fs




"""
# Messing around with some simulation. You can ignore

nTimesteps = 1000
nParticles = 5
nDims = 2
dt = 1e-15
L = sigma*1e1
boxDimensions = L*np.ones(nDims)
v0 = 0.01*sigma/dt
rss = np.zeros((nTimesteps, nParticles, nDims))
vss = np.zeros((nTimesteps, nParticles, nDims))

r0s = L*(np.random.rand(nParticles, nDims))
v0s = v0*(np.random.rand(nParticles, nDims)-0.5)
rss[0,:,:] = r0s
vss[0,:,:] = v0s

for i in range(nTimesteps-1):
    Fs = calculateForces(rss[i,:,:], boxDimensions=boxDimensions, nDims=nDims)
    rss[i+1,:,:] = rss[i,:,:] + dt * vss[i,:,:]
    vss[i+1,:,:] = vss[i,:,:] + dt * Fs/mass
    rss[i+1,:,:] = np.mod(rss[i+1,:,:],boxDimensions)



for j in range(nParticles):
    plt.scatter(rss[:,j,0], rss[:,j,1], alpha = np.linspace(0,0.1,nTimesteps))
    plt.scatter(rss[0,j,0], rss[0,j,1], s = 200, facecolors = "none", edgecolors = "red")
    plt.scatter(rss[-1,j,0], rss[-1,j,1], s = 200, facecolors = "none", edgecolors = "black")
plt.show()


"""