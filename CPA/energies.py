"""
This code provides functions to calculate the kinetic, potential and total energies

All functions, quantities are given in dimensionless units
"""


# IMPORTING
import numpy as np
from forces import delta_pos_MIC





def calculate_kinetic_energy(vel):
    """
    Calculates the kinetic energy of all particles

    Arguments:
        vel: np.ndarray (n_particles, n_dimensions), dtype = float
            array of velocity vectors of all particles
        
    Returns:
        kinetic_energy: float
            total kinetic energy
    """
    # elementwise multiplication to get the square of the velocity
    # sum to get total kinetic energy T
    return 0.5 * np.sum(np.multiply(vel, vel))





def calculate_pairwise_potential(delta_pos):
    """
    Calculates the pairwise Lennard-Jones (LJ) potential energy

    Arguments:
        delta_pos: np.ndarray (n_dimensions), dtype = float
            separation vector between target and interacting particles
    
    Returns: 
        pairwise_potential_energy: float
            pairwise potential energy
    """
    # calculate norm of r vector
    delta_pos_norm = np.sqrt(np.dot(delta_pos, delta_pos))
    # Calculate potential using LJ
    pairwise_potential_energy = 4 * (delta_pos_norm**-12 - delta_pos_norm**-6)
    return pairwise_potential_energy




def calculate_net_potential(delta_pos):
    """
    Calculates the net LJ potential energy experienced by a single particle from all interacting particles
    
    Arguments:
        delta_pos: np.ndarray (n_particles - 1, n_dimensions), dtype = float
            array of separation vectors between target particle and interacting particles
    
    Returns:
        net_potential_energy: float
            potential energy experienced by a single particle
    """
    # sum over all pairwise potential contributions
    net_potential_energy = 0
    for i in range(len(delta_pos)):
        net_potential_energy += calculate_pairwise_potential(delta_pos[i])
    return net_potential_energy




def calculate_potential_energy(pos, box_dimensions):
    """
    Calculates total potential energy experienced by all particles
    
    Arguments:
        pos: np.ndarray (n_particles, n_dimensions), dtype = float
            array of position vectors of all particles
        box_dimensions: np.ndarray (n_dimensions), dtype = float
            size of the box x, y, z
            
    Returns:
        sum_potential_energy: float
            total potential energy of all particles
    """
    # loop through each particle i
    sum_potential_energy = 0
    for i in range(len(pos)):
        # the following steps calculates separation vectors between particle i and all other particles
        # 1. take the difference in position between particle i and each other particle
        # 2. remove the zero vector corresponding to self interaction
        # 3. convert seperations into MIC nearest clone seperations delta_pos.
        delta_pos = delta_pos_MIC(np.delete(pos[i] - pos, i, 0), box_dimensions)
        # calculate net potential energy of particle i via Lennard-Jones potential
        sum_potential_energy += 0.5*calculate_net_potential(delta_pos)
    return sum_potential_energy



def total_energy(pos, vel, box_dimensions):
    """
    Calculates total energy of all particles
    
    Arguments:
        pos: np.ndarray (n_particles, n_dimensions), dtype = float
            array of position vectors of all particles
        vel: np.ndarray (n_particles, n_dimensions), dtype = float
            array of velocity vectors of all particles
        box_dimensions: np.ndarray (n_dimensions), dtype = float
            size of the box x, y, z

    Returns:
        total_energy: float
            total energy of all particles
    """
    return calculate_potential_energy(pos, box_dimensions) + calculate_kinetic_energy(vel)




def array_of_energies(pos, vel, box_dimensions):
    """
    Calculates potential, kinetic and total energy of all particles. Separate from total_energy function for efficiency.
    
    Arguments:
        pos: np.ndarray (n_particles, n_dimensions), dtype = float
            array of position vectors of all particles
        vel: np.ndarray (n_particles, n_dimensions), dtype = float
            array of velocity vectors of all particles
        box_dimensions: np.ndarray (n_dimensions), dtype = float
            size of the box x, y, z
    
    Returns:
        potential_energy: float
            total potential energy
        kinetic_energy: float
            total kinetic energy
        total_energy: float
            total energy
    """
    U = calculate_potential_energy(pos, box_dimensions)
    T = calculate_kinetic_energy(vel)
    E = U + T
    return U, T, E


