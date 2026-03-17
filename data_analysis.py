import matplotlib.pyplot as plt
import numpy as np


###############################################################
"""
RADIAL CORRELATION FUNCTION STATISTICS
Convert list of g(r) arrays into a numpy array
Shape: (num_runs, nBins)
"""
def radial_corr_stats(radialCorrelationDensitiess, rBinss, measure_times,
                      num_runs, number_density, d_less_T, field,
                      save=False, mean_std=True, show=True):

    corr_array = np.array(radialCorrelationDensitiess)

    # Use the r bins from the first simulation (they should all match)
    r_bins = np.array(rBinss[0])

    # Compute the mean value of g(r) for each bin across simulations
    corr_mean = np.mean(corr_array, axis=0)

    # Compute the standard deviation across simulations
    corr_std = np.std(corr_array, axis=0)

    # Compute the standard error of the mean
    if mean_std == True:
        corr_std = corr_std / np.sqrt(num_runs)

    ###############################################################
    # PLOT: individual simulations + average

    plt.figure(figsize=(8, 5))

    # Plot each simulation in light gray
    for i in range(num_runs):
        plt.plot(rBinss[i], radialCorrelationDensitiess[i],
                 color="gray", alpha=0.4)

    # Plot the averaged correlation function
    plt.plot(r_bins, corr_mean,
             color="red",
             linewidth=2,
             label="Average g(r)")

    # Plot the statistical uncertainty band
    plt.fill_between(r_bins,
                     corr_mean - corr_std,
                     corr_mean + corr_std,
                     alpha=0.2,
                     label="±1 std")

    plt.xlabel("r")
    plt.ylabel("g(r)")
    plt.title(f"Average radial correlation function ({num_runs} simulations)")
    plt.legend()
    plt.grid(True)

    # Save the plot
    if save:
        plt.savefig(
            f"Average_radial_correlation_function_rho_{number_density}_T_{d_less_T}_field_{field}.png",
            dpi=150,
            bbox_inches="tight"
        )
    if show:
        plt.show()

    plt.close()

    ###############################################################
    # SAVE NUMERICAL DATA
    # Columns: r, mean g(r), std
    if save == True:
        output = np.column_stack((r_bins, corr_mean, corr_std))

        np.savetxt(
            f"Average_radial_correlation_function_rho_{number_density}_T_{d_less_T}_field_{field}.txt",
            output,
            header="r_bins corr_mean corr_std"
        )

    return corr_mean, corr_std


###############################################################
"""
PRESSURE STATISTICS
"""
def press_stats(pressures, measure_times, num_runs, number_density,
                d_less_T, field, save=False, mean_std=True, show=True):

    # Convert pressure list to numpy array
    pressures_array = np.array(pressures)

    # Compute mean pressure
    pressure_mean = np.mean(pressures_array)

    # Compute standard deviation
    pressure_std = np.std(pressures_array)

    # Compute standard error of the mean
    if mean_std == True:
        pressure_std = pressure_std / np.sqrt(num_runs)

    ###############################################################
    # PLOT: pressure values from each simulation

    plt.figure(figsize=(8, 5))

    x = np.arange(num_runs)

    # Scatter plot of pressure values from each run
    plt.scatter(x, pressures_array,
                label="Pressure from individual simulations")

    # Plot the mean pressure as a horizontal dashed line
    plt.axhline(pressure_mean,
                color="red",
                linestyle="--",
                label=f"Mean pressure = {pressure_mean:.3f}")

    # Plot the uncertainty band (± standard deviation)
    plt.fill_between(x,
                     pressure_mean - pressure_std,
                     pressure_mean + pressure_std,
                     alpha=0.2,
                     label=f"Std deviation = {pressure_std:.3f}")

    plt.xlabel("Simulation index")
    plt.ylabel("Pressure")
    plt.title(f"Pressure measurements over {num_runs} simulations")
    plt.legend()
    plt.grid(True)

    # Save the plot
    if save:
        plt.savefig(
            f"Pressure_statistics_rho_{number_density}_T_{d_less_T}_field_{field}.png",
            dpi=150,
            bbox_inches="tight"
        )
    if show:
        plt.show()

    plt.close()

    ###############################################################
    # SAVE PRESSURE DATA
    pressure_output = np.column_stack((x, pressures_array))

    if save == True:
        np.savetxt(
            f"Pressure_values_rho_{number_density}_T_{d_less_T}_field_{field}.txt",
            pressure_output,
            header="simulation_index pressure"
        )

    return pressure_mean, pressure_std


###############################################################
# PLOT: pressure vs time or magnetic field
def pressure_vs_x_analysis(mean_pressure, std_pressure, x_values,
                           field_as_x=False, save=False, show=True):

    mean_pressure = np.asarray(mean_pressure).flatten()
    std_pressure = np.asarray(std_pressure).flatten()
    x_values = np.asarray(x_values).flatten()

    # Compute global mean and std over x
    pressure_mean_global = np.mean(mean_pressure)
    pressure_std_global = np.std(mean_pressure)

    if field_as_x == True:
        x_label = "magnetic field"
        title = "Average pressure as a function of magnetic field"
        file_tag = "field"
        header_x = "magnetic_field"
    else:
        x_label = "time"
        title = "Average pressure as a function of time"
        file_tag = "time"
        header_x = "measure_times"

    plt.figure(figsize=(8, 5))

    # Plot the averaged pressure vs x
    plt.plot(x_values, mean_pressure,
             color="red",
             linewidth=2,
             marker='o',
             markersize=4,
             label="Average pressure")

    # Local uncertainty band
    plt.fill_between(x_values,
                     mean_pressure - std_pressure,
                     mean_pressure + std_pressure,
                     alpha=0.2,
                     label="±1 std (local)")

    # Global mean line
    plt.axhline(pressure_mean_global,
                color="black",
                linestyle="--",
                label=f"Mean pressure = {pressure_mean_global:.3f}")

    # Global std band
    plt.fill_between(x_values,
                     pressure_mean_global - pressure_std_global,
                     pressure_mean_global + pressure_std_global,
                     alpha=0.1,
                     label=f"Global std = {pressure_std_global:.3f}")

    plt.xlabel(x_label)
    plt.ylabel("pressure")
    plt.title(title)
    plt.legend()
    plt.grid(True)

    # Save the plot
    if save:
        plt.savefig(
            f"Average_pressure_vs_{file_tag}.png",
            dpi=150,
            bbox_inches="tight"
        )
    if show:
        plt.show()

    plt.close()

    ###############################################################
    # SAVE NUMERICAL DATA
    if save == True:
        output = np.column_stack((x_values, mean_pressure, std_pressure))

        np.savetxt(
            f"Average_pressure_vs_{file_tag}.txt",
            output,
            header=f"{header_x} mean_pressure std_pressure"
        )

    return pressure_mean_global, pressure_std_global


###############################################################
# PLOT: peak position of g(r) vs time or magnetic field
def corr_vs_x_analysis(rad_corr_dens, rBinss, x_values,
                       num_runs, number_density, d_less_T, field,
                       field_as_x=False, save=False, mean_std=True, show=True):

    x_values = np.asarray(x_values).flatten()

    # Find the r_bin corresponding to the maximum of g(r) for each measurement
    peak_positions = []

    for g_r, r_bins in zip(rad_corr_dens, rBinss):
        g_r = np.asarray(g_r).flatten()
        r_bins = np.asarray(r_bins).flatten()

        max_idx = np.argmax(g_r)
        peak_positions.append(r_bins[max_idx])

    peak_positions = np.asarray(peak_positions)

    ###############################################################
    # Simple statistics on peak positions
    peak_mean = np.mean(peak_positions)
    peak_std = np.std(peak_positions)

    # Compute standard error of the mean
    if mean_std == True:
        peak_std = peak_std / np.sqrt(num_runs)

    if field_as_x == True:
        x_label = "magnetic field"
        title = "Position of the maximum of g(r) as a function of magnetic field"
        file_tag = "field"
        header_x = "magnetic_field"
    else:
        x_label = "time"
        title = "Position of the maximum of g(r) as a function of time"
        file_tag = "time"
        header_x = "measure_times"

    ###############################################################
    # PLOT: peak position vs x
    plt.figure(figsize=(8, 5))

    plt.plot(x_values, peak_positions,
             color="red",
             linewidth=2,
             marker='o',
             markersize=4,
             label="Peak position")

    plt.axhline(peak_mean,
                color="black",
                linestyle="--",
                label=f"Mean peak position = {peak_mean:.3f}")

    plt.fill_between(x_values,
                     peak_mean - peak_std,
                     peak_mean + peak_std,
                     alpha=0.2,
                     label=f"Std = {peak_std:.3f}")

    plt.xlabel(x_label)
    plt.ylabel("r at max g(r)")
    plt.title(title)
    plt.legend()
    plt.grid(True)

    # Save the plot
    if save:
        plt.savefig(
            f"Peak_position_vs_{file_tag}_rho_{number_density}_T_{d_less_T}_field_{field}.png",
            dpi=150,
            bbox_inches="tight"
        )

    if show:
        plt.show()

    plt.close()

    ###############################################################
    # SAVE NUMERICAL DATA
    if save == True:
        output = np.column_stack((x_values, peak_positions))

        np.savetxt(
            f"Peak_position_vs_{file_tag}_rho_{number_density}_T_{d_less_T}_field_{field}.txt",
            output,
            header=f"{header_x} peak_positions"
        )

    return peak_positions, peak_mean, peak_std