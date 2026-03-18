"""
This code provides functions used to calculate the observables of the simulation.
These are the correlation function and the pressure

All quantities are given in natural units
"""

# IMPORTING
import numpy as np
from forces import delta_pos_MIC




def calculate_pairwise_distances(pos, box_dimensions):
    """
    Calculates all pairwise distances between all particles

    Arguments:
        pos: np.ndarray (n_particles, n_dimensions), dtype = float
            array of all particle position vectors
        box_dimensions: np.ndarray(n_dimensions), dypte = float
            size of the box x, y, z
    
    Returns:
        all_delta_pos_norm: np.ndarray (n_particles*(n_particles-1)), dtype = float
            array of all particle separation distances
    """
    n_particles = len(pos)
    all_distances=[]
    for i in range(n_particles):
        # 1. take the difference in position between particle i and each other particle
        # 2. remove the zero vector corresponding to self interaction
        # 3. convert seperations into MIC nearest clone seperations.
        delta_pos = delta_pos_MIC(np.delete(pos[i]-pos, i, 0), box_dimensions)
        delta_pos_norm = np.linalg.norm(delta_pos, axis = 1)
        all_distances.extend(delta_pos_norm)
    return np.array(all_distances)




def calculate_pressure(pos, temperature, box_dimensions): 
    """
    Calculates total pressure
    
    Arguments:
        pos: np.ndarray(n_particles, n_dimensions), dtype = float
            array of position vectors of all particles
        temperature: float
        box_dimensions: np.ndarray(n_dimensions), dtype = float
            size of the box x, y, z
        """
    # number of particles and particle number density
    n_particles = len(pos)
    number_density = n_particles/np.prod(box_dimensions)
    # Get all pairwise distances from the calculatePairwiseDistances function
    all_distances = calculate_pairwise_distances(pos,box_dimensions)
    # sum over all contributions to pressure for each pairwise distance
    sum = 0
    for i in range(len(all_distances)):
        pairwiseDistance = all_distances[i]
        sum += -48*pairwiseDistance**-12 + 24*pairwiseDistance**-6
    # Put the sum into the formula for pressure (Verlet)
    P = temperature*number_density*(1-sum/(12*n_particles*temperature))
    return P
    



"""
Calculates the radial correlation function density g(r) and corresponding distance r ([nBins], [nBins]).
pos is an array of all positions within the box [nParticles, nDimensions]
boxDimensions is the x,y,z size array of the box. [nDimensions]
"""
def calculate_correlation_function(pos, box_dimensions, n_bins):
    """
    Calculates the radial correlation function based on the dimensionality of the box
    
    Arguments:
        pos: np.ndarray (n_particles, n_dimensions), dtype = float
            array of position vectors of all particles
        box_dimensions: np.ndarray (n_dimensions), dtype = float
            size of the box x, y, z
        n_bins: int
            number of bins for the histogram = number of distance sampling points

    Returns:
        radial_correlation_densities: np.ndarray (n_bins), dtype = float
            radial correlation function, normalised
        bin_distances: np.ndarray (nbins), dtype = float
            distance centers of the bins
    """
    # number of particles and volume
    n_particles = len(pos)
    volume = np.prod(box_dimensions)

    # Maximum valid distance in the histogram is L/2, beyond that g(r) is unphysical
    max_distance = np.min(box_dimensions)/2

    # get all pairwise distances between particles
    all_distances = calculate_pairwise_distances(pos,box_dimensions)

    # calculate absolute correlation histogram
    absolute_correlation_function, bin_edges = np.histogram(all_distances, bins=n_bins, range=[0,max_distance])
    # convert bin edges into bin centers
    bin_distances = (bin_edges[:-1] + bin_edges[1:])/2
    bin_separations = bin_distances[1]-bin_distances[0]               # distance spacing between bins

    # Convert absolute correlation to the proper radial correlation function
    if len(box_dimensions) == 1:    # 1D
        radial_correlation_densities = 2*volume*absolute_correlation_function/(n_particles*(n_particles-1)*2*bin_separations)

    elif len(box_dimensions) == 2:   # 2D
        radial_correlation_densities = 2*volume*absolute_correlation_function/(n_particles*(n_particles-1)*2*np.pi*bin_distances*bin_separations)

    elif len(box_dimensions) == 3:     # 3D
        radial_correlation_densities = 2*volume*absolute_correlation_function/(n_particles*(n_particles-1)*4*np.pi*bin_distances**2*bin_separations)
    else:
        print("Not appropriate dimensionality. 1D, 2D or 3D only.")
    return radial_correlation_densities, bin_distances


