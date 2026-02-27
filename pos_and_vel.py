import numpy as np

"""
in this code we will do all the calculation to setup
-box (array of shaoe=(N_dim,) of lenght L)
-position (assuring that the atoms will be equally distributed and not overlapping)
-velocity (2d gaussian distribution with < ||v||^2 > = 1 )
"""
n_particles_1d = 1
n_dim = 2

empty_space = 1

L = n_particles_1d * (np.sqrt(n_dim) + empty_space)
n_particles = n_particles_1d * n_dim


def box_array(n_dim, L):
    return L * np.ones(n_dim)


# this function computes the position of all particles in a matrix (n_particles_1d,n_dim)
# the particles will spawn in suboxes with some empty space between each other to prevent overlap
def position(n_dim, empty_space, n_particles_1d, L):

    n_particles = n_particles_1d**n_dim

    subbox_size_1d = (L - (n_particles_1d) * empty_space) / n_particles_1d

    if n_dim not in (2, 3):
        raise ValueError("n_dim must be 2 or 3")

    if subbox_size_1d < 0:
        raise ValueError(
            "the spawn point is <0 so the particles will not spawn in the correct point"
        )

    pitch = subbox_size_1d + empty_space

    # set where all the particles will spawn inside its own sub-box
    pos = np.random.uniform(0, subbox_size_1d, size=(n_particles, n_dim))

    # computation of the index for each particles 0 1 2       0 0 0
    #                                            0 1 2  and  1 1 1
    #                    in 2d:                  0 1 2       2 2 2

    if n_dim == 2:
        i_x, i_y = np.meshgrid(np.arange(n_particles_1d), np.arange(n_particles_1d))
        offset = np.column_stack([i_x.ravel(), i_y.ravel()]) * pitch
    elif n_dim == 3:
        ix, iy, iz = np.meshgrid(
            np.arange(n_particles_1d),
            np.arange(n_particles_1d),
            np.arange(n_particles_1d),
        )
        offset = np.column_stack([ix.ravel(), iy.ravel(), iz.ravel()]) * pitch

    pos += offset + empty_space / 2

    return pos


def velocity(n_particles_1d, n_dim, mean=0, std=1):
    n_particles = n_particles_1d**n_dim

    velocity_module = np.random.normal(mean, std, size=(n_particles,n_dim))
    """theta = np.random.uniform(0, 2 * np.pi, size=n_particles)
    velocity = np.zeros((n_particles, n_dim))
    velocity[:, 0] = velocity_module * np.cos(theta)
    velocity[:, 1] = velocity_module * np.sin(theta)"""
    return velocity_module


def toy_position(n_dim, empty_space, n_particles_1d, L):

    pos = np.zeros((2, n_dim))
    pos[0] += np.ones(n_dim) * L / 3
    pos[1] += np.ones(n_dim) * 2 * L / 3

    return pos


def toy_velocity(n_particles_1d, n_dim, pos):

    vel = np.zeros((2, n_dim))
    vel[0] += pos[1] - pos[0]
    vel[1] -= pos[1] - pos[0]

    return vel


def renormalization( T, kinetic):
    T_dim_less=T/100       #temperature in dimensionless unit   
    kin_target=3 * (n_particles - 1) * T_dim_less   
    factor=np.sqrt(kin_target/kinetic)
    return factor
    

def stable(kinetics):
    if kinetics.shape[0]<30:
        return False
    elif np.abs(np.mean(np.diff(kinetics)))<0.1:
        return True
    else:
        return False