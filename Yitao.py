

import numpy as np
import matplotlib.pyplot as plt



# physical parameters corresponding to Argon
sigma = 3.5e-10
epsilon = 120*1.38e-23
mass = 40e-3/(6.022e-23)



"""
Calculates interaction force between two particles.
r is the vector from target -> interacting particle.
returns the force experienced from the interacting particle
"""
def pairwiseForce(r):
    # calculate the norm of rVec = r        
    rNorm = np.sqrt(np.dot(r, r))
    # calculate the force via F = -nabla U using the Lennard-Jones potential
    F = epsilon*(48*(sigma**12)*(rNorm**-14) - 24*(sigma**6)*(rNorm**-8))*r
    return F



"""
Calculates net force of one particle experienced from all other particles
rs is an array of vectors from target -> interacting particles
nDims is the number of dimensions.
returns the net force experienced from all the interactions on that particle
"""
def netForce(rs, nDims)
    # loop through each particle and and up its pairwise force contribution
    totalForce = np.zeros(nDims)
    for i in range(len(rs)):
        totalForce += pairwiseForce(rs[i])
    return totalForce

"""
Converts a seperation vector between particles inside a box, to the seperation in the Minimum Image Convention (MIC) clone
r is the seperation vector between to particles
boxDimensions is the x,y,z size array of the box. E.g. a cubic box has (L,L,L)
"""
def rMIC(r, boxDimensions):
    # shift the vector r, then take the modulus, which returns the right periodic image.
    # Then re-shift r back to origin.
    return np.mod(r + 0.5*boxDimensions, boxDimensions) - 0.5*boxDimensions




"""
THIS IS THE FUNCTION YOU USE
Calculates all forces of all particles
rs is an array of all positions within the box
boxDimensions is the x,y,z size array of the box. E.g. a cubic box has (L,L,L)
Returns an array of forces
"""
def calculateForces(rs, boxDimensions, nDims): 
    # loop through each particle i
    fs = np.zeros(len(rs),nDims)
    for i in range(len(rs)):
        # 1. take the difference in position between particle i and each other particle
        # 2. remove the zero vector corresponding to self interaction
        # 3. convert seperations into MIC nearest clone seperations.
        deltaRs = rMIC(np.delete(rs[i]-rs, i, 0), boxDimensions)
        # calculate net force on particle i via Lennard Jones potential
        fs[i] = netForce(deltaRs, nDims)
    return fs






# Messing around with some simulation. You can ignore

nTimesteps = 1000
nParticles = 6
nDims = 2
dt = 0.1
L = sigma
boxDimensions = L*np.ones(nDims)
v0 = L*0.001
rss = np.zeros((nTimesteps, nParticles, nDims))
vss = np.zeros((nTimesteps, nParticles, nDims))

r0s = L*(np.random.rand(nParticles, nDims))
v0s = 2*v0*(np.random.rand(nParticles, nDims)-0.5)
rss[0,:,:] = r0s
vss[0,:,:] = v0s

for i in range(nTimesteps-1):
    Fs = np.zeros((nParticles,nDims))
    for j in range(nParticles):
        deltaRList = np.ndarray.tolist(rss[i,:,:]-rss[i,j,:])
        del deltaRList[j]
        deltaRs = np.array(deltaRList)
        deltaRs = MICneighbour(deltaRs)
        Fs[j,:] = totalInteractionForce(deltaRs)
    rss[i+1,:,:] = rss[i,:,:] + dt * vss[i,:,:]
    vss[i+1,:,:] = vss[i,:,:] + dt * Fs/mass
    rss[i+1,:,:] = np.mod(rss[i+1,:,:],boxDimensions)



for j in range(nParticles):
    plt.scatter(rss[:,j,0], rss[:,j,1], alpha = np.linspace(0,1,nTimesteps))
plt.show()


