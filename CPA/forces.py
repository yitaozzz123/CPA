"""
This code contains all the functions for calculating the total force experienced by particles

The final function used is the all_forces function

All quantities are given in natural units
"""


# IMPORTING
import numpy as np


# physical parameters corresponding to Argon. These are 1 in dimensionless units.
sigma = 1
epsilon = 1
mass = 1



def pairwise_force(delta_pos):
    """
    Calculates the pairwise force experienced by a target particle from an interacting particle

    Arguments:
        delta_pos: np.ndarray (n_dimensions), dtype = float
            vector from target to interacting particle
    
    Returns:
        force: np.ndarray (n_dimensions), dtype = float
            force vector acting on target particle
    """
    # calculate the norm of delta_pos vector    
    delta_pos_norm = np.sqrt(np.dot(delta_pos, delta_pos))
    # calculate the force via F = -nabla U using the Lennard-Jones potential
    force = (48*(delta_pos_norm**-14) - 24*(delta_pos_norm**-8))*delta_pos
    return force



def net_force(delta_pos, n_dims, external_field):
    """
    Calculates the net force experienced by a target particle from all other particles

    Arguments:
        delta_pos: np.ndarray (n_particles - 1, n_dims), dtype = float
            array of vectors from target to interacting particles
        n_dims: int
            number of dimensions
        external_field: np.ndarray (n_dims), dtype = float
            vector of external field
    
    Returns:
        force: np.ndarray (n_dims), dtype = float
            net force experienced by target particle from all interactions and external field
    """
    # loop through each particle and sum up its pairwise force contribution
    force = np.zeros(n_dims)
    for i in range(len(delta_pos)):
        force += pairwise_force(delta_pos[i])
    force += external_field
    return force




def delta_pos_MIC(delta_pos, box_dimensions):
    """
    Converts a delta_pos separation vector between particles into the minimum image convention (MIC) separation vector

    Arguments:
        delta_pos: np.ndarray (n_dims), dtype = float
            vector between to particles without MIC
        box_dimensions: np.ndarray (n_dims), dtype = float
            size of the box x, y, z
    """
    # shift the vector r, then take the modulus, which returns the right periodic image.
    # Then re-shift r back to origin.
    return np.mod(delta_pos + 0.5*box_dimensions, box_dimensions) - 0.5*box_dimensions





def calculate_forces(pos, box_dimensions, n_dims, external_field = 0): 
    """
    Calculates the total force experienced by all particles
    This is the function used in simulation
    
    Arguments:
        pos: np.ndarray (n_particles, n_dims), dtype = float
            array of position vectors of all particles
        box_dimensions: np.ndarray (n_dims), dtype = float
            size of the box x, y, z
        n_dims: int
            number of dimensions
        external_field: np.ndarray (n_dims), dtype = float
            vector of external field

    Returns:
        all_forces: np.ndarray (n_particles, n_dims), dtype = float
            array of all force vectors for all particles
    """
    # loop through each particle i
    # fs is the an array of net-force vectors for each particle 
    fs = np.zeros((len(pos),n_dims))
    for i in range(len(pos)):
        # the following steps calculates separation vectors between particle i and all other particles
        # 1. take the difference in position between particle i and each other particle
        # 2. remove the zero vector corresponding to self interaction
        # 3. convert seperations into MIC nearest clone seperations delta_pos.
        delta_pos = delta_pos_MIC(np.delete(pos[i]-pos, i, 0), box_dimensions)
        # calculate net force on particle i via Lennard Jones potential
        fs[i] = net_force(delta_pos, n_dims, external_field)
    return fs



