import numpy as np
"""
in this code we will do all the calculation to setup
-box (array of shaoe=(N_dim,) of lenght L)
-position (assuring that the atoms will be equally distributed and not overlapping)
-velocity (2d gaussian distribution with < ||v||^2 > = 1 )
"""
n_particles_1d=3
n_dim=2

L=n_particles_1d * np.sqrt(n_dim)
n_particles=n_particles_1d*n_dim

def box_size(n_dim, L):
    return L * np.ones(n_dim)


def position(n_dim, box_size, epty_space=1):
    return 0

def velocity(n_particles, n_dim):
    velocity_module=np.random.normal(1, 1, size=n_particles)
    theta=np.random.uniform(0, 2*np.pi, size=n_particles)
    velocity=np.zeros((n_particles, n_dim))
    velocity[:,0]=velocity_module*np.cos(theta)
    velocity[:,1]=velocity_module*np.sin(theta)
    return velocity

velocit=velocity(n_particles, n_dim)
mod_vel=np.sqrt(velocit[:,0]**2+velocit[:,1]**2)
print(mod_vel)