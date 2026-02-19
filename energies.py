import numpy as np
from forces import rMIC

"""
All functions, quantities, use dimensionless quantities
"""

"""
Calculate kinetic energy
vs is an array of velocity vectors (nParticles x nDims matrix).
"""
def calculateKineticEnergy(vs):
    # elementwise multiplication to get the square of the velocity
    # sum to get total kinetic energy T
    kineticEnergy = 0.5*np.sum(np.multiply(vs,vs))
    return kineticEnergy

"""
Calculates pairwise Lennard-Jones (LJ) potential energy
r is a single distance vector
"""
def pairwisePotential(r):
    # calculate norm of r vector
    rNorm = np.sqrt(np.dot(r, r))
    # Calculate potential using LJ
    U = 4*(rNorm**-12 - rNorm**-6)
    return U

"""
Calculates net Lennard-Jones potential energy of a single particle
rs is an array of distance vectors (nParticles x nDims matrix).
"""
def netPotential(rs):
    netU = 0
    for i in range(len(rs)):
        netU += pairwisePotential(rs[i]) 
    return netU

"""
Calculate total potential energy of all particles
rs is an array of all particle positions (nParticles x nDims matrix).
"""
def calculatePotentialEnergy(rs, boxDimensions):
    # loop through each particle i
    potentialEnergy = 0
    for i in range(len(rs)):
        # 1. take the difference in position between particle i and each other particle
        # 2. remove the zero vector corresponding to self interaction
        # 3. convert seperations into MIC nearest clone seperations.
        deltaRs = rMIC(np.delete(rs[i]-rs, i, 0), boxDimensions)
        # calculate net potential energy of particle i via Lennard-Jones potential
        potentialEnergy += netPotential(deltaRs)
    return potentialEnergy

"""
Calculates total energy of all particles
rs and vs are arrays of all particle positions and velocities respectively (nParticles x nDims matrices).
"""
def calculateTotalEnergy(rs, vs, boxDimensions):
    potentialEnergy = calculatePotentialEnergy(rs, boxDimensions)
    kineticEnergy = calculateKineticEnergy(vs)
    totalEnergy = potentialEnergy + kineticEnergy
    return (totalEnergy)


def array_of_energies(rs, vs, boxDimensions):
    potentialEnergy = calculatePotentialEnergy(rs, boxDimensions)/2
    kineticEnergy = calculateKineticEnergy(vs)
    totalEnergy = potentialEnergy + kineticEnergy
    return potentialEnergy, kineticEnergy, totalEnergy