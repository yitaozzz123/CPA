import matplotlib.pyplot as plt
import numpy as np

########################################################################
    # PLOT OF ENERGIES
def plot_energies_fluctuations(timestep, total, field, show=False, save=False):

    x = np.arange(len(total) - 1)
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

    if save==True:
        plt.savefig(f"Energy_fluctuation_field_{field}.png", dpi=150, bbox_inches="tight")
    
    if show==True:
        plt.show()
    
    plt.close() 

    return None

def plot_energies(kinetic, potential, total, kin_target, timestep, field, show=False, save=False):

    x = np.arange(len(kinetic))

    plt.close("all")

    plt.figure(figsize=(8, 5))

    plt.plot(x, kinetic, label="Kinetic energy")
    plt.plot(x, potential, label="Potential energy")
    plt.plot(x, total, label="Total energy")
    plt.axhline(y=float(kin_target),linestyle="--", label="Energy target")

    plt.xlabel("Iteration number")
    plt.ylabel("Energy")
    plt.title(f"Energy plot with t = {timestep}")

    plt.legend()
    plt.grid(True)
    if save==True:
        plt.savefig(f"Energy_field_{field}.png", dpi=150, bbox_inches="tight")
    if show==True:
        plt.show()
    plt.close()