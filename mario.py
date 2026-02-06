import numpy as np
import scipy as sp
from scipy import sparse
from scipy.spatial import cKDTree
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from collections import deque


box_size=10
n_dim=2
num_particles=100
mass=0.1
timestep=0.01


rss = np.zeros((num_particles, n_dim))
pos=np.random.uniform(0, box_size, size=(num_particles, n_dim) )
vel=np.random.uniform(-10, 10, size=(num_particles, n_dim) )
tail=deque(maxlen=6)
tail.append(pos.copy())

for i in range(1000):

    #computation of the force
    F=np.zeros((num_particles, n_dim))

    acceleration=F/mass
    pos+=vel*timestep
    vel+=acceleration*timestep

    pos[pos>box_size] -= box_size
    pos[pos<0] += box_size

    tail.append(pos.copy())

    hist = np.stack(tail, axis=0)

    plt.clf()

    if len(tail)==6 and n_dim==2:
        plt.scatter(pos[:,0],pos[:,1],marker='o')

        for p in range(num_particles):
            plt.plot(hist[:,p,0],hist[:,p,1])
    
        plt.xlim(0, box_size)
        plt.ylim(0, box_size)
        plt.pause(0.01)





