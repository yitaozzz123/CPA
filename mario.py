import numpy as np
import scipy as sp
from scipy import sparse
from scipy.spatial import cKDTree
import matplotlib.pyplot as plt
import matplotlib.animation as animation

box_size=10
n_dim=2
N=10

pos = np.random.uniform(0,box_size,size=(N,n_dim))
orient = np.random.uniform(-np.pi, np.pi,size=N)

fig, ax= plt.subplots(figsize=(6,6))

qv = ax.quiver(pos[:,0], pos[:,1], np.cos(orient[0]), np.sin(orient), orient, clim=[-np.pi, np.pi])