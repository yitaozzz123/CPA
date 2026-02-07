import numpy as np
import scipy as sp
from scipy import sparse
from scipy.spatial import cKDTree
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from collections import deque

#Numerical setup
box_size=10
n_dim=2
num_particles=100
mass=0.1
timestep=0.01
tail_lenght=6
num_iterations=100

#initialization of position, velocity 
rss = np.zeros((num_particles, n_dim))
pos=np.random.uniform(0, box_size, size=(num_particles, n_dim) )
vel=np.random.uniform(-10, 10, size=(num_particles, n_dim) )

#tail initialization
tail=deque(maxlen=tail_lenght)
tail.append(pos.copy())

for i in range(num_iterations):

    #computation of the force and acceleration
    F=np.zeros((num_particles, n_dim))    
    acceleration=F/mass

    #update of positions and velocity for every interaction
    pos+=vel*timestep
    vel+=acceleration*timestep

    #application of the periodic boundary conditions
    pos[pos>box_size] -= box_size
    pos[pos<0] += box_size

    #tail update
    tail.append(pos.copy())

    #from deque to numpy to have a functioning plot
    tail_numpy = np.stack(tail, axis=0)

    plt.clf()

    #we do not want to plot the tail of a particle if the periodic boundary condition happened 
    d = np.diff(tail_numpy, axis=0)                          
    wrapped = np.any(np.abs(d) > box_size/2, axis=2) 
    yes_tail_index = ~np.any(wrapped, axis=0)                


    if len(tail)==tail_lenght and n_dim==2:
        plt.scatter(pos[:,0],pos[:,1],marker='o')
            
        plottable_tail=tail_numpy[:,yes_tail_index,:]
        
        for plottable_particle in range(np.shape(plottable_tail)[1]):
            plt.plot(plottable_tail[:,plottable_particle,0],plottable_tail[:,plottable_particle,1])
    
        plt.xlim(0, box_size)
        plt.ylim(0, box_size)
        plt.pause(0.01)





