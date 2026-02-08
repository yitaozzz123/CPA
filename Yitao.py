import numpy as np
import matplotlib.pyplot as plt

sigma = 3.5e-10
epsilon = 120*1.38e-23*200
mass = 40/(6.022e-23)


# returns force from LJ potential with rVec: target -> interacting particle
def singleInteractionForce(rVec):           # rVec is the seperation from target particle to surrounding particles
    r = np.sqrt(np.dot(rVec, rVec))         # norm of rVec
    F = epsilon*(48*(sigma**12)*(r**-14) - 24*(sigma**6)*(r**-8))*rVec
    return F
    # single interaction force. For total, sum up all interactions

# takes a list of vectors from target to interacting particles, returns total force
def totalInteractionForce(rs):
    totalForce = np.zeros(nDims)
    for i in range(len(rs)):
        totalForce += singleInteractionForce(rs[i])
    return totalForce


# Finds nearest copy for minimum interaction convention
def MICneighbour(rVec):
    return np.mod(rVec + 0.5*boxDimensions, boxDimensions) - 0.5*boxDimensions  # boxDimenions is the x,y,z vector of the size of the box. E.g. if its a square box of sizeN then boxDimensions = L*np.ones(nDims)
    # returns the conversion to nearest periodic clone









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


