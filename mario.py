import numpy as np
import scipy as sp
from scipy import sparse
from scipy.spatial import cKDTree
import matplotlib.pyplot as plt
import matplotlib.animation as animation

box_size=10
n_dim=2
num_particles=10
mass=0.1
timestep=0.01


rss = np.zeros((num_particles, n_dim))
pos=np.random.uniform(0, box_size, size=(num_particles, n_dim) )
vel=np.random.uniform(-1, 1, size=(num_particles, n_dim) )

def update(pos, vel, time, timestep=timestep):
    #computation of the force
    F=np.zeros((num_particles, n_dim))

    acceleration=F/mass
    pos+=vel*timestep
    vel+=acceleration*timestep

    pos[pos>box_size] -= box_size
    pos[pos<0] += box_size



