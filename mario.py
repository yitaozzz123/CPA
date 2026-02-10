import numpy as np
import scipy as sp
from scipy import sparse
from scipy.spatial import cKDTree
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from collections import deque
from Yitao import calculateForces


"""
In this code, the numerical simulation of the molecule dinamic is animated and saved

"""

"""
This first part contains the numerical setup of the simulation: 
- box_size, num_particles, num_dim
    mass, are the physical settings of the experiment.
- timestep is used to compute the update of the position and the velocity
- tail_lenght determine how long is the tail, that is the component the makes the
    direction of the particles visible
-num_iteractions determine how long will be the simulation
- save (boolean) is a switch used to save the data of position, velocities, tails 
    and the final plot
"""


box_size=1e-9
n_dim=2
num_particles=20 #200
mass = 40/(6.022e+23)
timestep=1e-18
tail_lenght=50
num_iterations=500
save=False
fps=60

#initialization of position, velocity 
pos=np.random.uniform(0, box_size, size=(num_particles, n_dim))
vel=np.random.uniform(-box_size/(1000*timestep), box_size/(timestep*1000), size=(num_particles, n_dim))

box=box_size * np.ones(n_dim)

#tail initialization
tail=deque(maxlen=tail_lenght)
tail.append(pos.copy())

#setup for 3d
if n_dim==3:
    import matplotlib.pyplot as plt

    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(111, projection="3d")


for i in range(num_iterations):

    #computation of the force and acceleration
    F=calculateForces(rs=pos, boxDimensions=box, nDims=n_dim)  
    acceleration=F/mass

    #update of positions and velocity for every interaction
    pos+=vel*timestep
    vel+=acceleration*timestep

    #application of the periodic boundary conditions
    
    #first implementation of periodic boundary
    #pos[pos>box_size] -= box_size
    #pos[pos<0] += box_size

    #second implementation
    #pos_int=pos//box_size
    #pos-=pos_int*box_size

    #third implementation
    pos%=box_size

    #tail update
    tail.append(pos.copy())

    #from deque to numpy to have a functioning plot
    tail_numpy = np.stack(tail, axis=0)

    #we do not want to plot the tail of a particle if the periodic boundary condition happened 
    d = np.diff(tail_numpy, axis=0)                          
    wrapped = np.any(np.abs(d) > box_size/2, axis=2) 
    yes_tail_index = ~np.any(wrapped, axis=0)                

    #selection of the suitable tails
    plottable_tail=tail_numpy[:,yes_tail_index,:]
    """
        mod_vel=np.zeros((vel.shape[0]))
        mod_vel=np.dot(vel,vel)
        max_vel=np.max(mod_vel)
        print(max_vel*timestep/sigma)"""

    #animation in 2d or 3d and save of the last plot
    if len(tail)==tail_lenght and n_dim==2:

        plt.clf()

        plt.scatter(pos[:,0],pos[:,1],marker='o')
        
        for plottable_particle in range(np.shape(plottable_tail)[1]):
            plt.plot(plottable_tail[:,plottable_particle,0],plottable_tail[:,plottable_particle,1])
    
        plt.xlim(0, box_size)
        plt.ylim(0, box_size)
        plt.pause(1/fps)

        if save==True:
            np.save("tail.npy", plottable_tail)
            np.save("pos.npy", pos)
            np.save("vel.npy",vel)
            plt.savefig("figure.png", dpi=150, bbox_inches="tight")

    if len(tail) == tail_lenght and n_dim == 3:
        ax.cla()

        ax.scatter(pos[:, 0], pos[:, 1], pos[:, 2], marker='o')

        for plottable_particle in range(plottable_tail.shape[1]):
            ax.plot(plottable_tail[:, plottable_particle, 0],
                    plottable_tail[:, plottable_particle, 1],
                    plottable_tail[:, plottable_particle, 2])

        ax.set_xlim(0, box_size)
        ax.set_ylim(0, box_size)
        ax.set_zlim(0, box_size)

        plt.draw()
        plt.pause(1/fps)

        if save == True:
            np.save("tail.npy", plottable_tail)
            np.save("pos.npy", pos)
            np.save("vel.npy", vel)
            fig.savefig("figure.png", dpi=150, bbox_inches="tight")

