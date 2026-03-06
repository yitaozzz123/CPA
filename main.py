import numpy as np
import scipy as sp
from scipy import sparse
from scipy.spatial import cKDTree
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from collections import deque
from forces import calculateForces
from pos_and_vel import box_array, position, velocity, toy_position, toy_velocity, renormalization, stable, FCC_pos, stable
import time
from energies import array_of_energies
from pressure import calculatePressure

"""
In this code, the numerical simulation of the molecule dinamic is animated and saved
"""


"""
This first part contains the numerical setup of the simulation: 
- L, num_particles, num_dim
    mass, are the physical settings of the experiment.
- timestep is used to compute the update of the position and the velocity
- tail_lenght determine how long is the tail, that is the component the makes the
    direction of the particles visible
-num_iteractions determine how long will be the simulation
- save (boolean) is a switch used to save the data of position, velocities, tails 
    and the final plot
"""
##################################################
#SIMULATION PHYSICAL PARAMETERS

n_dim=3
number_density = 1.2 # Dimensionless units!
d_less_T=0.5        #dless_T=T/120K

#################################################
# TIME
timestep = 1e-2
fps = 120
max_simulation_time = 60

tot_internal_time = 2
num_iterations = int(tot_internal_time / timestep)


##################################################
# BOX PARAMETERS 
""" NOT NEEDED ANYMORE!
empty_space = 1
ratio = 1  # ratio>=1
L = n_particles_1d * (ratio * empty_space)
# L=20
box = box_array(n_dim, L)
"""

####################################################
# FEATURES
tail_lenght = 20
full_tail = False

save = True
save_data = 1
plot_en_fluct = 0
toy_model = 0
############################################
"""
#PROPER INTIALIZATION OF THE SIMULATION PARAMETERS
#POSITION
#VELOCITY
#TAIL
#ENERGIES
"""
def simulation():
    # initialization of position, velocity depending on toy model switch
    if toy_model:
        pos = toy_position(n_dim, L=5)
        vel = toy_velocity(n_dim, pos)
    else:
        pos, box =FCC_pos(number_density)
        L = box[0]
        n_particles = len(pos)
        vel=velocity(n_particles, n_dim, mean=0, std=np.sqrt(d_less_T))

    # tail initialization
    if full_tail == False:
        tail = deque(maxlen=tail_lenght)
    else:
        tail = []

    tail.append(pos.copy())

    # setup for 3d animation & plot
    if n_dim == 3:
        import matplotlib.pyplot as plt

        fig = plt.figure(figsize=(8, 8))
        ax = fig.add_subplot(111, projection="3d")

    # initialization arrays related to energies
    potential = []
    kinetic = []
    total = []

    # first computation of energies
    energy = array_of_energies(pos=pos, vel=vel, boxDimensions=box)
    potential.append(energy[0])
    kinetic.append(energy[1])
    total.append(energy[2])

    ##################################################################
    # SIMULATION
    count=0
    # starting simulation time
    start_time = time.time()

    for i in range(num_iterations):
        # If we exceed with the simulation time, the simulation will stop
        if max_simulation_time < (time.time() - start_time):
            break

        internal_time = timestep * i

        #############################
        # pos-pos'<1.1*sigma
        ##################################################################
        # UPDATE OF THE PARAMETERS

        F = calculateForces(pos=pos, boxDimensions=box, nDims=n_dim)

        pos += vel * timestep + (timestep**2) * F / 2

        F_2 = calculateForces(pos=pos, boxDimensions=box, nDims=n_dim)

        vel += timestep * (F_2 + F) / 2 #industrial freezer effect

        # application of the periodic boundary conditions
        pos %= L

        # computation of energy
        energy = array_of_energies(pos=pos, vel=vel, boxDimensions=box)
        potential.append(energy[0])
        kinetic.append(energy[1])
        total.append(energy[2])

        # tail update
        tail.append(pos.copy())

        if stable(kinetic):
            if count<50:
                factor, kin_target=renormalization(d_less_T,energy[1],number_density, L*L*L)
                vel*=factor
            if count==100:
                pressure=calculatePressure(pos, d_less_T, box, n_dim)  
                #corr_function= calculate_corr_function(*variables) 
                break      
            count+=1

        ####################################################################
        # TAIL FIXES ON ARRAY HANDLING

        # from deque to numpy to have a functioning plot
        tail_numpy = np.stack(tail, axis=0)

        # we do not want to plot the tail of a particle if the periodic boundary condition happened
        if full_tail == False:
            d = np.diff(tail_numpy, axis=0)
            wrapped = np.any(np.abs(d) > L / 2, axis=2)
            yes_tail_index = ~np.any(wrapped, axis=0)
            plottable_tail = tail_numpy[:, yes_tail_index, :]
        else:
            plottable_tail = tail_numpy

        #########################################################################
        # STARTING ANIMATIONS

        # animation in 2d or 3d and save of the last plot
        if n_dim == 2:
            plt.clf()

            plt.scatter(pos[:, 0], pos[:, 1], marker="o")

            for plottable_particle in range(np.shape(plottable_tail)[1]):
                plt.plot(
                    plottable_tail[:, plottable_particle, 0],
                    plottable_tail[:, plottable_particle, 1],
                )

            plt.xlim(0, L)
            plt.ylim(0, L)
            plt.title(f"simulation time t={internal_time:.3f}")
            plt.pause(1 / fps)

            if save == True:
                np.save("tail.npy", plottable_tail)
                np.save("pos.npy", pos)
                np.save("vel.npy", vel)
                plt.savefig("2D_plot.png", dpi=150, bbox_inches="tight")

        if n_dim == 3:
            ax.cla()

            ax.scatter(pos[:, 0], pos[:, 1], pos[:, 2], marker="o")

            for plottable_particle in range(plottable_tail.shape[1]):
                ax.plot(
                    plottable_tail[:, plottable_particle, 0],
                    plottable_tail[:, plottable_particle, 1],
                    plottable_tail[:, plottable_particle, 2],
                )

            ax.set_xlim(0, L)
            ax.set_ylim(0, L)
            ax.set_zlim(0, L)

            plt.title(f"simulation time t={internal_time:.3f}")
            plt.draw()
            plt.pause(1 / fps)

            if save_data == True:
                np.save("tail.npy", plottable_tail)
                np.save("pos.npy", pos)
                np.save("vel.npy", vel)

            if save == True:
                fig.savefig("3D_plot.png", dpi=150, bbox_inches="tight")


    ########################################################################
    # PLOT OF ENERGIES
    if plot_en_fluct == True:
        x = np.arange(len(kinetic) - 1)
        total = np.diff(total)

        plt.close("all")

        plt.figure(figsize=(8, 5))

        # plt.plot(x, kinetic, label="Kinetic energy")
        # plt.plot(x, potential, label="Potential energy")
        plt.plot(x, total, label="Total energy")

        plt.xlabel("Iteration number")
        plt.ylabel("Energy")
        plt.title(f"Energy plot with t = {timestep}")

        plt.legend()
        plt.grid(True)
        plt.savefig(f"Energy_fluctuation_{n_dim}D.png", dpi=150, bbox_inches="tight")
        #plt.show()

    x = np.arange(len(kinetic))

    plt.close("all")

    plt.figure(figsize=(8, 5))

    plt.plot(x, kinetic, label="Kinetic energy")
    plt.plot(x, potential, label="Potential energy")
    plt.plot(x, total, label="Total energy")
    plt.axhline(y=float(kin_target), label="Energy target")

    plt.xlabel("Iteration number")
    plt.ylabel("Energy")
    plt.title(f"Energy plot with t = {timestep}")

    plt.legend()
    plt.grid(True)
    plt.savefig(f"Energy_fluctuation_{n_dim}D.png", dpi=150, bbox_inches="tight")
    #plt.show()

    return pressure #, corr_function


##################################################################################################

pressures, corr_funcs=[],[]

for i in range(5):
    pressure=simulation()
    pressures.append(pressure)
    #corr_funcs.append(corr_funct)

pressure_mean=np.mean(pressures)
pressure_std=np.std(pressures)

x = np.arange(len(pressures))

plt.close("all")
plt.figure()
plt.scatter(x,pressures, label="Pressure on single simulation")
plt.axhline(pressure_mean, color="red", label=f"Mean: {pressure_mean:.3f}")

plt.fill_between(x,pressure_mean-pressure_std, pressure_mean+pressure_std,
                 alpha=0.2, label=f"std deviation: {pressure_std:.3f}")

plt.xlabel("Simulation number")
plt.ylabel("Pressure")
plt.legend()
plt.title(f"Pressures on {len(x)} simulations")
plt.savefig(f"Pressures on {len(x)} simulations.png")


"""
corr_mean=np.mean(corr_funcs)
corr_std=np.std(corr_funcs)

x = np.arange(len(corr_funcs))

plt.close("all")
plt.figure()
plt.scatter(x,corr_funcs, label="Correlation function on single simulation")
plt.axhline(corr_mean, color="red", label=f"Mean: {corr_mean:.3f}")

plt.fill_between(x,corr_mean-corr_std, corr_mean+corr_std,
                 alpha=0.2, label=f"std deviation: {corr_std:.3f}")

plt.xlabel("Simulation number")
plt.ylabel("Correlation function")
plt.legend()
plt.title(f"Correlation functions on {x} simulations")
plt.savefig(f"Correlations on {x} simulations.png")
"""