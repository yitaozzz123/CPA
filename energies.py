import numpy as np
from forces import rMIC

"""
All functions, quantities, use dimensionless quantities
"""

"""
Calculate kinetic energy (float)
vel is an array of velocity vectors [nParticles, nDimensions].
"""
def calculateKineticEnergy(vel):
    # elementwise multiplication to get the square of the velocity
    # sum to get total kinetic energy T
    kineticEnergy = 0.5*np.sum(np.multiply(vel,vel))
    return kineticEnergy

"""
Calculates pairwise Lennard-Jones (LJ) potential energy (float)
deltaR is a single distance vector [nDimensions]
"""
def pairwisePotential(deltaR):
    # calculate norm of r vector
    rNorm = np.sqrt(np.dot(deltaR, deltaR))
    # Calculate potential using LJ
    U = 4*(rNorm**-12 - rNorm**-6)
    return U

"""
Calculates net Lennard-Jones potential energy of a single particle (float)
deltaPos is an array of distance vectors [nParticles , nDimensions].
"""
def netPotential(deltaPos):
    netU = 0
    for i in range(len(deltaPos)):
        netU += pairwisePotential(deltaPos[i]) 
    return netU

"""
Calculate total potential energy of all particles (float)
pos is an array of all particle positions [nParticles, nDimensions].
"""
def calculatePotentialEnergy(pos, boxDimensions):
    # loop through each particle i
    potentialEnergy = 0
    for i in range(len(pos)):
        # 1. take the difference in position between particle i and each other particle
        # 2. remove the zero vector corresponding to self interaction
        # 3. convert seperations into MIC nearest clone seperations.
        deltaPos = rMIC(np.delete(pos[i]-pos, i, 0), boxDimensions)
        # calculate net potential energy of particle i via Lennard-Jones potential
        potentialEnergy += netPotential(deltaPos)
    return potentialEnergy

"""
Calculates total energy of all particles (float)
pos and vel are arrays of all particle positions and velocities respectively [nParticles, nDimensions]
"""
def calculateTotalEnergy(pos, vel, boxDimensions):
    potentialEnergy = calculatePotentialEnergy(pos, boxDimensions)
    kineticEnergy = calculateKineticEnergy(vel)
    totalEnergy = potentialEnergy + kineticEnergy
    return (totalEnergy)

"""
Calculates potential, kinetic and total energy of all particles (float, float, float)
pos and vel are arrays of all particle positions and velocities respectively [nParticles, nDimensions]
"""
def array_of_energies(pos, vel, boxDimensions):
    potentialEnergy = calculatePotentialEnergy(pos, boxDimensions)/2
    kineticEnergy = calculateKineticEnergy(vel)
    totalEnergy = potentialEnergy + kineticEnergy
    return potentialEnergy, kineticEnergy, totalEnergy