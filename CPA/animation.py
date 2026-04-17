import matplotlib.pyplot as plt
import numpy as np

def animation(plt, fig, ax, plottable_tail, L, internal_time, fps, pos, vel, save_data=True, save=True):

    # animation in 2d or 3d and save of the last plot
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

    return None