import matplotlib.pyplot as plt
import numpy as np


###############################################################
"""
RADIAL CORRELATION FUNCTION STATISTICS
Convert list of g(r) arrays into a numpy array
Shape: (num_runs, nBins)
"""

###############################################################
# RADIAL CORRELATION FUNCTION STATISTICS VS FIELD
def radial_corr_stats(radialCorrelationDensitiess, rBinss,
                      num_runs, number_density, d_less_T, field,
                      save=False, show=True):
    """
    Compute mean and standard deviation of g(r) over independent runs
    at a fixed field value.

    Parameters
    ----------
    radialCorrelationDensitiess : list of arrays
        List of g(r) arrays from different simulation runs.
        Expected shape: (num_runs, nBins)

    rBinss : list of arrays
        Radial bins corresponding to each run.
        Assumed identical for all runs.

    num_runs : int
        Number of independent simulations.

    Returns
    -------
    r_bins : ndarray
    corr_mean : ndarray
        Mean radial correlation function.
    corr_std : ndarray
        Standard deviation across runs (NOT standard error).
    """

    # Convert input list to numpy array
    corr_array = np.array(radialCorrelationDensitiess)

    # Use r bins from first run (assumed identical)
    r_bins = np.array(rBinss[0])

    # Compute mean and standard deviation across runs
    corr_mean = np.mean(corr_array, axis=0)
    corr_std = np.std(corr_array, axis=0)

    ###############################################################
    # PLOT: mean g(r) + standard deviation band
    plt.figure(figsize=(3.0, 2.2))

    plt.plot(
        r_bins, corr_mean,
        linewidth=2,
        label="Mean g(r)"
    )

    plt.fill_between(
        r_bins,
        corr_mean - corr_std,
        corr_mean + corr_std,
        alpha=0.25,
        label="±1 std"
    )

    plt.xlabel("r [σ]", fontsize=11)
    plt.ylabel("g(r)", fontsize=11)
    #plt.title(f"Radial correlation function ({num_runs} runs), field={field}")
    plt.legend(fontsize=9)
    plt.grid(True, linewidth=0.4, alpha=0.6)
    
    if save:
        plt.savefig(
            f"Radial_correlation_mean_std_rho_{number_density}_T_{d_less_T}_field_{field}.pdf",
            dpi=150,
            bbox_inches="tight"
        )

    if show:
        plt.show()

    plt.close()

    ###############################################################
    # SAVE NUMERICAL DATA
    if save:
        output = np.column_stack((r_bins, corr_mean, corr_std))
        np.savetxt(
            f"Radial_correlation_mean_std_rho_{number_density}_T_{d_less_T}_field_{field}.txt",
            output,
            header="r_bins corr_mean corr_std"
        )

    return corr_mean, corr_std


###############################################################
# PRESSURE STATISTICS VS FIELD
def press_stats(pressures, num_runs, number_density, d_less_T, field,
                save=False, show=True):
    """
    Compute mean and standard deviation of pressure across runs
    at a fixed field.

    Parameters
    ----------
    pressures : list or array
        Pressure values from independent runs.

    num_runs : int
        Number of simulations.

    Returns
    -------
    pressure_mean : float
    pressure_std : float
        Standard deviation across runs.
    """

    pressures_array = np.array(pressures)

    # Compute statistics
    pressure_mean = np.mean(pressures_array)
    pressure_std = np.std(pressures_array)

    ###############################################################
    # PLOT: individual runs + mean + std band
    x=np.arange(num_runs)
    plt.figure(figsize=(3.0, 2.2))

    plt.scatter(x, pressures_array, s=25, label="Runs")

    plt.axhline(
        pressure_mean,
        linestyle="--",
        linewidth=1.5,
        label=f"Mean = {pressure_mean:.4f}"
    )

    plt.fill_between(
        x,
        pressure_mean - pressure_std,
        pressure_mean + pressure_std,
        alpha=0.3,
        label=f"±1 std = {pressure_std:.4f}"
    )

    plt.xlabel("Simulation index", fontsize=11)
    plt.ylabel("Pressure", fontsize=11)

    plt.legend(fontsize=9)
    plt.grid(True, linewidth=0.5, alpha=0.6)

    if save:
        plt.savefig(
            f"Pressure_statistics_rho_{number_density}_T_{d_less_T}_field_{field}.pdf",
            dpi=150,
            bbox_inches="tight"
        )

    if show:
        plt.show()

    plt.close()

    ###############################################################
    # SAVE NUMERICAL DATA
    if save:
        # Save individual pressure values
        pressure_values_output = np.column_stack((x, pressures_array))
        np.savetxt(
            f"Pressure_values_rho_{number_density}_T_{d_less_T}_field_{field}.txt",
            pressure_values_output,
            header="simulation_index pressure"
        )

        # Save mean and std
        pressure_stats_output = np.array([[pressure_mean, pressure_std]])
        np.savetxt(
            f"Pressure_mean_std_rho_{number_density}_T_{d_less_T}_field_{field}.txt",
            pressure_stats_output,
            header="pressure_mean pressure_std"
        )

    return pressure_mean, pressure_std


###############################################################
# PRESSURE VS FIELD
def pressure_vs_field_analysis(mean_pressure, std_pressure, field_values,
                               save=False, show=True):
    """
    Plot mean pressure as a function of field with standard deviation.

    Parameters
    ----------
    mean_pressure : array-like
        Mean pressure for each field.

    std_pressure : array-like
        Standard deviation for each field.

    field_values : array-like
        Field values.
    """

    mean_pressure = np.asarray(mean_pressure).flatten()
    std_pressure = np.asarray(std_pressure).flatten()
    field_values = np.asarray(field_values).flatten()

    plt.figure(figsize=(3.0, 2.2))

    plt.plot(
        field_values, mean_pressure,
        linewidth=1.5,
        marker='o',
        markersize=4,
        label="Mean pressure"
    )

    plt.fill_between(
        field_values,
        mean_pressure - std_pressure,
        mean_pressure + std_pressure,
        alpha=0.3,
        label="±1σ"
    )

    plt.xlabel("Field", fontsize=11)
    plt.ylabel("Pressure", fontsize=11)

    plt.legend( fontsize=9)
    plt.grid(True, linewidth=0.5, alpha=0.6)

    plt.tight_layout()

    if save:
        plt.savefig(
            "Pressure_vs_field.pdf",
            dpi=150,
            bbox_inches="tight"
        )

    if show:
        plt.show()

    plt.close()

    ###############################################################
    # SAVE NUMERICAL DATA
    if save:
        output = np.column_stack((field_values, mean_pressure, std_pressure))
        np.savetxt(
            "Pressure_vs_field.txt",
            output,
            header="field mean_pressure std_pressure"
        )

    return mean_pressure, std_pressure


###############################################################
# CORRELATION FUNCTION VS FIELD
def corr_vs_field_analysis(corr_means, r_binss, field_values,
                           save=False, show=True):
    """
    Plot mean g(r) for different field values on the same graph.

    No standard deviation is shown here to allow a clearer comparison
    of structural changes.

    Parameters
    ----------
    corr_means : list of arrays
        Mean g(r) for each field.

    r_binss : list of arrays
        Corresponding r bins.

    field_values : array-like
        Field values.
    """

    field_values = np.asarray(field_values).flatten()

    plt.figure(figsize=(3.0, 2.2))

    for corr_mean, r_bins, field in zip(corr_means, r_binss, field_values):
        corr_mean = np.asarray(corr_mean).flatten()
        r_bins = np.asarray(r_bins).flatten()

        plt.plot(
            r_bins, corr_mean,
            linewidth=2,
            label=f"field = {field}"
        )

    plt.xlabel("r [σ]",fontsize=11)
    plt.ylabel("g(r)",fontsize=11)
    #plt.title("Mean radial correlation function vs field")
    plt.legend(fontsize=9)
    plt.grid(True)

    if save:
        plt.savefig(
            "Radial_correlation_vs_field.pdf",
            dpi=150,
            bbox_inches="tight"
        )

    if show:
        plt.show()

    plt.close()

    ###############################################################
    # SAVE DATA (using npz for multiple curves)
    if save:
        save_dict = {"field_values": field_values}

        for i, (r_bins, corr_mean) in enumerate(zip(r_binss, corr_means)):
            save_dict[f"r_bins_{i}"] = np.asarray(r_bins)
            save_dict[f"corr_mean_{i}"] = np.asarray(corr_mean)

        np.savez("Radial_correlation_vs_field.npz", **save_dict)

    return corr_means, r_binss, field_values


###############################################################
# PLOT: pressure vs time or field (just time)
def pressure_vs_x_analysis(mean_pressure, std_pressure, x_values,
                           field_as_x=False, save=False, show=True):

    mean_pressure = np.asarray(mean_pressure).flatten()
    std_pressure = np.asarray(std_pressure).flatten()
    x_values = np.asarray(x_values).flatten()

    # Compute global mean and std over x
    pressure_mean_global = np.mean(mean_pressure)
    pressure_std_global = np.std(mean_pressure)

    if field_as_x == True:
        x_label = "field"
        title = "Average pressure as a function of field"
        file_tag = "field"
        header_x = "field"
    else:
        x_label = "time"
        title = "Average pressure as a function of time"
        file_tag = "time"
        header_x = "measure_times"

    plt.figure(figsize=(3.0, 2.2))

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

    plt.xlabel(x_label,fontsize=11)
    plt.ylabel("pressure",fontsize=11)
    #plt.title(title)
    plt.legend(fontsize=9)
    plt.grid(True)

    # Save the plot
    if save:
        plt.savefig(
            f"Average_pressure_vs_{file_tag}.pdf",
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
# PLOT: peak position of g(r) vs time or field (just for time)
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
        x_label = "field"
        title = "Position of the maximum of g(r) as a function of field"
        file_tag = "field"
        header_x = "field"
    else:
        x_label = "time"
        title = "Position of the maximum of g(r) as a function of time"
        file_tag = "time"
        header_x = "measure_times"

    ###############################################################
    # PLOT: peak position vs x
    plt.figure(figsize=(3.0, 2.2))

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

    plt.xlabel(x_label,fontsize=11)
    plt.ylabel("r at max g(r)",fontsize=11)
    #plt.title(title)
    plt.legend(fontsize=9)
    plt.grid(True)

    # Save the plot
    if save:
        plt.savefig(
            f"Peak_position_vs_{file_tag}_rho_{number_density}_T_{d_less_T}_field_{field}.pdf",
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

