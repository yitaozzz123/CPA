import numpy as np
import scipy as sp
from scipy import sparse
from scipy.spatial import cKDTree
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from collections import deque
from forces import calculateForces
from pos_and_vel import box_array, position, velocity
import time
from energies import array_of_energies

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

simulation_time=10
mass = 1
timestep=1e-3
tail_lenght=50  
num_iterations=500
save=True
fps=60
periodic_tail=True
save_data=False

n_particles_1d=5
n_dim=3

empty_space=0.3

#ratio=0.5
#L=n_particles_1d * (ratio * empty_space)
L=10
n_particles=n_particles_1d**n_dim

box_size=L
#initialization of position, velocity 
pos=position(n_dim, empty_space, n_particles_1d, L)

vel=velocity(n_particles_1d, n_dim)

box = box_array(n_dim, L)

#tail initialization
if periodic_tail==False:
    tail=deque(maxlen=tail_lenght)
else:
    tail=[]

tail.append(pos.copy())     

#setup for 3d animation & plot
if n_dim==3:
    import matplotlib.pyplot as plt

    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(111, projection="3d")

#initialization arrays related to energy
potential=[]
kinetic=[]
total=[]

#first computation of energy
energy=array_of_energies(rs=pos, vs=vel, boxDimensions=box)
potential.append(energy[0])
kinetic.append(energy[1])
total.append(energy[2])


#starting simulation time
start_time=time.time()

while True:

    if simulation_time < (time.time()-start_time):
        break

    #computation of the force and acceleration
    F=calculateForces(rs=pos, boxDimensions=box, nDims=n_dim)  
    acceleration=F/mass

    #update of positions and velocity for every interaction
    pos+=vel*timestep
    vel+=acceleration*timestep

    #application of the periodic boundary conditions
    pos%=box_size

    #tail update
    tail.append(pos.copy())

    #from deque to numpy to have a functioning plot
    tail_numpy = np.stack(tail, axis=0)

    #we do not want to plot the tail of a particle if the periodic boundary condition happened 
    if periodic_tail==False:
        d = np.diff(tail_numpy, axis=0)                          
        wrapped = np.any(np.abs(d) > box_size/2, axis=2) 
        yes_tail_index = ~np.any(wrapped, axis=0) 
        plottable_tail=tail_numpy[:,yes_tail_index,:]
    else:
        plottable_tail=tail_numpy
        
    #computation of energy
    energy=array_of_energies(rs=pos, vs=vel, boxDimensions=box)
    potential.append(energy[0])
    kinetic.append(energy[1])
    total.append(energy[2])


    #animation in 2d or 3d and save of the last plot
    if n_dim==2:

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
            plt.savefig("2D_plot.png", dpi=150, bbox_inches="tight")

    if n_dim == 3:
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

        if save_data == True:
            np.save("tail.npy", plottable_tail)
            np.save("pos.npy", pos)
            np.save("vel.npy", vel)

        if save==True:
            fig.savefig("£D_plot.png", dpi=150, bbox_inches="tight")

x = np.arange(len(kinetic))

plt.close('all')

plt.figure(figsize=(8,5))

plt.plot(x, kinetic, label="Kinetic energy")
plt.plot(x, potential, label="Potential energy")
plt.plot(x, total, label="Total energy")

plt.xlabel("iteration number")
plt.ylabel("energy")
plt.title(f"Energy plot with t = {timestep}")

plt.legend()
plt.grid(True)
plt.savefig(f"Energy_{n_dim}D.png", dpi=150, bbox_inches="tight")
plt.show()
