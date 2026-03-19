"""
This code provides the backbone of all the project: it simulates the moldecule dynamic
using most of the tool to: 
    - initialize the correct initial conditions*
    - compute internal quantities* e.g. forces, energies
    - update variables
    - apply renormalization of velocity*
    - prepare observables for the final analysis*

* using importied function
"""


import numpy as np
import matplotlib.pyplot as plt
from collections import deque

from forces import calculate_forces
from pos_and_vel import FCC_pos, velocity, renormalization, stable
from energies import array_of_energies
from observables import calculate_pressure, calculate_correlation_function
import matplotlib.pyplot as plt
from animation import animation
from energies_plot import plot_energies_fluctuations, plot_energies




def simulation(number_density, d_less_T, num_iterations,
               timestep, field, n_counts, field_module, animate=False,
               plot_fluctuations=False, fps=60, tail_lenght=20, save=False, field_study_mode=False, show=False):
    """
    Run one molecular-dynamics simulation

    The function initializes the system, equilibrates it, optionally applies
    an external field, and stores pressure, radial correlation function,
    their bins, and measurement times

    Arguments: (all in dimension-less units)
        number_density: number density of the system
        d_less_T: temperature
        num_iterations: number of maximum steps in the for loop
        timestep: timestep
        field: if True, apply an external field after equilibration
        n_counts: number of measurements to collect in time-study mode
        field_module: strength of the external field
        animate: if True, show the animation
        plot_fluctuations: if True, plot energy fluctuations
        fps: frames per second for the animation (not that usefull)
        tail_lenght: tail length used in the animation (fixes some problems with the
          animation of periodic boundary conditions)
        save: if True, save plots and data
        field_study_mode: if True, return a single post-field measurement for
            field-scan studies

    Returns:
        pressures, radialCorrelationDensitiess, rBinss, measure_time
    """

    ##############################################################
    #INITIALIZATION PHASE
    #############################################################
    equilibrium=False
    
    pressures, radialCorrelationDensitiess, rBinss, measure_time = [],[],[],[]
    
    # initialization arrays related to energies
    potential, kinetic, total = [], [], []

    # initialization of position, velocity depending on toy model switch
    pos, box =FCC_pos(number_density)
    L = box[0]
    n_particles = len(pos)
    vel=velocity(n_particles, n_dim=3, mean=0, std=np.sqrt(d_less_T))

    # tail initialization
    tail = deque(maxlen=tail_lenght)

    tail.append(pos.copy())

    # setup for 3d animation & plot
    if animate == True:
        fig = plt.figure(figsize=(8, 8))
        ax = fig.add_subplot(111, projection="3d")
    else:
        fig = None
        ax = None

    # first computation of energies
    energy = array_of_energies(pos=pos, vel=vel, box_dimensions=box)
    potential.append(energy[0])
    kinetic.append(energy[1])
    total.append(energy[2])

    count=0
    field_count=0
    # starting simulation time
    field_start_time=0
    
    ##################################################################
    # SIMULATION LOOP
    #################################################################
    for i in range(num_iterations):

        internal_time = timestep * i
        ######################################################################
        # UPDATE OF THE PARAMETERS
        #####################################################################
        if equilibrium==True and field==True:
            F = calculate_forces(pos=pos, box_dimensions=box, n_dims=3, external_field=np.array([0,0,field_module]))
        else:
            F = calculate_forces(pos=pos, box_dimensions=box, n_dims=3)

        pos += vel * timestep + (timestep**2) * F / 2

        F_2 = calculate_forces(pos=pos, box_dimensions=box, n_dims=3)

        vel += timestep * (F_2 + F) / 2 #industrial freezer effect

        # application of the periodic boundary conditions
        pos %= L

        # computation of energy
        energy = array_of_energies(pos=pos, vel=vel, box_dimensions=box)
        potential.append(energy[0])
        kinetic.append(energy[1])
        total.append(energy[2])

        # tail update
        tail.append(pos.copy())

        ###################################################
        #APPLICATION OF RENORMALIZATION AND FIELD
        ####################################################
        #Going towards equilibrium: 
        # -  if kinetic energy is stable
        # -  application of renormalization 50 times
        # -  equilibrium is considered reached and eventual field can start
        # -  save pressure and correlation funct

        if stable(kinetic):
            if count<50:
                factor, kin_target=renormalization(d_less_T,energy[1],number_density, L*L*L)
                vel*=factor
            elif count==50:
                equilibrium=True
                field_start_time=i
                print("Equilibrium reached")
                if field==True:
                    print("Field started")
            elif count==100 and field==False:
                pressure=calculate_pressure(pos, d_less_T, box)  
                radialCorrelationDensities, rBins = calculate_correlation_function(pos, box_dimensions=box, n_bins=50) 
                pressures.append(pressure)
                radialCorrelationDensitiess.append(radialCorrelationDensities)
                rBinss.append(rBins)
                measure_time.append(internal_time)
                break      
            count+=1

        #with field on:
        #initialize field_time properly when field start to be applied
        # -  save every 20 time step pressure and correlation function
        field_time = timestep * (i-field_start_time)

        if (equilibrium == True) and (field == True) and field_study_mode==False:
            if field_count%20==0:
                pressure=calculate_pressure(pos, d_less_T, box)  
                radialCorrelationDensities, rBins = calculate_correlation_function(pos, box_dimensions=box, n_bins=50) 
                
                pressures.append(pressure)
                radialCorrelationDensitiess.append(radialCorrelationDensities)
                rBinss.append(rBins)
                measure_time.append(field_time)

                print(f"Sample number {len(pressures)}")

            if len(pressures)==n_counts:
                break

            field_count+=1

        if field_study_mode==True and equilibrium==True and field==True:
                if field_time>0.5:
                    pressure=calculate_pressure(pos, d_less_T, box)  
                    radialCorrelationDensities, rBins = calculate_correlation_function(pos, box_dimensions=box, n_bins=50) 
                    pressures.append(pressure)
                    radialCorrelationDensitiess.append(radialCorrelationDensities)
                    rBinss.append(rBins)
                    measure_time.append(field_time)
                    break   
                
        ####################################################################
        # TAIL LENGHT INDICES HANDLING FOR ANIMATION
        ##################################################################
        # from deque to numpy to have a functioning plot
        # we do not want to plot the tail of a particle if the periodic boundary condition happened
        tail_numpy = np.stack(tail, axis=0)
        d = np.diff(tail_numpy, axis=0)
        wrapped = np.any(np.abs(d) > L / 2, axis=2)
        yes_tail_index = ~np.any(wrapped, axis=0)
        plottable_tail = tail_numpy[:, yes_tail_index, :]

        if animate==True:
            animation(plt, fig, ax, plottable_tail, L, internal_time, fps, pos, vel, save_data=True, save=True)    


    ##############################################################################
    #PLOT OF ENERGIES
    #############################################################################
    if plot_fluctuations == True:
        plot_energies_fluctuations(timestep, total, field, show=show, save=save)

    plot_energies(kinetic, potential, total, kin_target, timestep, field, save=save, show=show)


    return pressures, radialCorrelationDensitiess, rBinss, measure_time