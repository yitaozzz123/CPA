"""Main script to run molecular dynamics studies and related data analysis.

This file provides three study modes:
- no external field;
- fixed external field with measurements over time;
- scan over different field values.

The simulation itself is implemented in `simulation.py`, while the plotting
and statistics are handled in `data_analysis.py`.
"""

import numpy as np

from data_analysis import (
    radial_corr_stats,
    press_stats,
    pressure_vs_x_analysis,
    corr_vs_x_analysis,
)
from simulation import simulation


##################################################
# SIMULATION PHYSICAL PARAMETERS + STATISTICS
#
# num_runs:
#     number of independent simulations used to estimate averages and
#     uncertainties on pressure and radial correlation function.
# number_density:
#     dimensionless particle density.
# d_less_T:
#     reduced temperature, defined here as T / 120 K.
# field:
#     switch for turning the external field on or off.
# field_study:
#     if False, study the evolution at fixed field over time;
#     if True, scan different field strengths.
# field_module:
#     fixed field strength used in the time-dependent study.
# n_counts:
#     number of measurements collected during a time study.
# field_max, n_field_values:
#     maximum field and number of sampled field values in the field scan.
##################################################

num_runs = 3
number_density = 1.2
d_less_T = 0.5

field = True
field_study = True
field_module = 50
n_counts = 5
field_max = 100
n_field_values = 5


#################################################
# TIME PARAMETERS
#
# timestep:
#     integration timestep for the equations of motion.
# tot_internal_time:
#     total simulated time.
# num_iterations:
#     total number of integration steps.
#################################################

timestep = 1e-2
tot_internal_time = 5
num_iterations = int(tot_internal_time / timestep)


####################################################
# EXECUTION / OUTPUT FEATURES
#
# save:
#     save plots and data to disk if True.
# plot_fluctuations:
#     save fluctuation plots if enabled in simulation.py.
# animate:
#     run the animation.
# show:
#     display plots interactively.
####################################################

save = False
plot_fluctuations = False
animate = False
show = True


def main_no_field(num_runs):
    """Run several simulations without external field and analyse final observables.

    Parameters
    ----------
    num_runs : int
        Number of independent simulations.

    Returns
    -------
    int
        Zero on successful completion.
    """
    pressures, radialCorrelationDensitiess, rBinss, measure_times = [], [], [], []

    for i in range(num_runs):
        print(f"Simulation {i} started")

        pressure, radialCorrelationDensities, rBins, measure_time = simulation(
            number_density,
            d_less_T,
            num_iterations,
            timestep,
            field,
            n_counts,
            field_module,
            animate=animate,
            plot_fluctuations=plot_fluctuations,
            save=save,
            field_study_mode=False,
        )

        pressures.append(pressure[0])
        radialCorrelationDensitiess.append(radialCorrelationDensities[0])
        rBinss.append(rBins[0])
        measure_times.append(measure_time[0])

    press_stats(
        pressures,
        measure_times,
        num_runs,
        number_density,
        d_less_T,
        field,
        save=save,
        show=show,
    )

    radial_corr_stats(
        radialCorrelationDensitiess,
        rBinss,
        measure_times,
        num_runs,
        number_density,
        d_less_T,
        field,
        save=save,
        show=show,
    )

    return 0


def main_time(num_runs, n_counts):
    """Run several simulations at fixed field and analyse observables versus time.

    Parameters
    ----------
    num_runs : int
        Number of independent simulations.
    n_counts : int
        Number of measurements collected during each simulation.

    Returns
    -------
    int
        Zero on successful completion.
    """
    pressures, radialCorrelationDensitiess, rBinss = [], [], []

    for i in range(num_runs):
        print(f"Simulation {i} started")

        pressure, radialCorrelationDensities, rBins, measure_time = simulation(
            number_density,
            d_less_T,
            num_iterations,
            timestep,
            field,
            n_counts,
            field_module,
            animate=animate,
            plot_fluctuations=plot_fluctuations,
            save=save,
            field_study_mode=False,
        )

        pressures.append(pressure)
        radialCorrelationDensitiess.append(radialCorrelationDensities)
        rBinss.append(rBins)

    mean_pressures, std_pressures = [], []
    mean_radialCorrelationDensitiess, std_radialCorrelationDensitiess = [], []

    pressures = np.array(pressures)
    radialCorrelationDensitiess = np.array(radialCorrelationDensitiess)
    rBinss = np.array(rBinss)
    measure_times = np.array(measure_time.copy())

    for j in range(len(pressures[0])):
        mean_pressure, std_pressure = press_stats(
            np.array(pressures[:, j]),
            measure_times,
            num_runs,
            number_density,
            d_less_T,
            field,
            save=save,
            show=False,
        )

        mean_radialCorrelationDensities, std_radialCorrelationDensities = radial_corr_stats(
            radialCorrelationDensitiess[:, j],
            rBinss[:, j],
            measure_times,
            num_runs,
            number_density,
            d_less_T,
            field,
            save=save,
            show=False,
        )

        mean_pressures.append(mean_pressure)
        std_pressures.append(std_pressure)
        mean_radialCorrelationDensitiess.append(mean_radialCorrelationDensities)
        std_radialCorrelationDensitiess.append(std_radialCorrelationDensities)

        print(j)

    mean_pressures = np.array(mean_pressures)
    std_pressures = np.array(std_pressures)
    mean_radialCorrelationDensitiess = np.array(mean_radialCorrelationDensitiess)

    pressure_vs_x_analysis(
        mean_pressures,
        std_pressures,
        measure_times,
        field_as_x=False,
        save=save,
        show=show,
    )

    corr_vs_x_analysis(
        mean_radialCorrelationDensitiess,
        rBinss[0],
        measure_times,
        num_runs,
        number_density,
        d_less_T,
        field,
        field_as_x=False,
        save=save,
        show=show,
    )

    return 0


def main_field(num_runs, field_max, n_field_values):
    """Run several simulations for different field values and analyse observables.

    Parameters
    ----------
    num_runs : int
        Number of independent simulation sets.
    field_max : float
        Maximum field value included in the scan.
    n_field_values : int
        Number of sampled field values, including 0 and field_max.

    Returns
    -------
    int
        Zero on successful completion.
    """
    pressures, radialCorrelationDensitiess, rBinss = [], [], []

    field_values = np.linspace(0.0, field_max, n_field_values)

    for i in range(num_runs):
        print(f"Set of simulation {i} started")

        run_pressures = []
        run_radialCorrelationDensitiess = []
        run_rBinss = []

        for j, field_value in enumerate(field_values):
            print(f"Step {j} started, field = {field_value}")

            if field_value==0:
                pressure, radialCorrelationDensities, rBins, measure_time = simulation(
                    number_density,
                    d_less_T,
                    num_iterations,
                    timestep,
                    field=False,
                    n_counts=1,
                    field_module=field_value,
                    animate=animate,
                    plot_fluctuations=plot_fluctuations,
                    save=save,
                    field_study_mode=False,
                )
            else:
                pressure, radialCorrelationDensities, rBins, measure_time = simulation(
                    number_density,
                    d_less_T,
                    num_iterations,
                    timestep,
                    field=True,
                    n_counts=1,
                    field_module=field_value,
                    animate=animate,
                    plot_fluctuations=plot_fluctuations,
                    save=save,
                    field_study_mode=True,
                )

            run_pressures.append(pressure[0])
            run_radialCorrelationDensitiess.append(radialCorrelationDensities[0])
            run_rBinss.append(rBins[0])

        pressures.append(run_pressures)
        radialCorrelationDensitiess.append(run_radialCorrelationDensitiess)
        rBinss.append(run_rBinss)

    mean_pressures, std_pressures = [], []
    mean_radialCorrelationDensitiess, std_radialCorrelationDensitiess = [], []

    pressures = np.array(pressures)
    radialCorrelationDensitiess = np.array(radialCorrelationDensitiess)
    rBinss = np.array(rBinss)

    for j in range(len(field_values)):
        mean_pressure, std_pressure = press_stats(
            pressures[:, j],
            field_values,
            num_runs,
            number_density,
            d_less_T,
            field=True,
            save=save,
            show=False,
        )

        mean_radialCorrelationDensities, std_radialCorrelationDensities = radial_corr_stats(
            radialCorrelationDensitiess[:, j],
            rBinss[:, j],
            field_values,
            num_runs,
            number_density,
            d_less_T,
            True,
            save=save,
            show=False,
        )

        mean_pressures.append(mean_pressure)
        std_pressures.append(std_pressure)
        mean_radialCorrelationDensitiess.append(mean_radialCorrelationDensities)
        std_radialCorrelationDensitiess.append(std_radialCorrelationDensities)

        print(j)

    mean_pressures = np.array(mean_pressures)
    std_pressures = np.array(std_pressures)
    mean_radialCorrelationDensitiess = np.array(mean_radialCorrelationDensitiess)
    std_radialCorrelationDensitiess = np.array(std_radialCorrelationDensitiess)

    pressure_vs_x_analysis(
        mean_pressures,
        std_pressures,
        field_values,
        field_as_x=True,
        save=save,
        show=show,
    )

    corr_vs_x_analysis(
        mean_radialCorrelationDensitiess,
        rBinss[0],
        field_values,
        num_runs,
        number_density,
        d_less_T,
        True,
        field_as_x=True,
        save=save,
        show=show,
    )

    return 0


def main():
    """Select and run the requested study mode.

    Returns
    -------
    int
        Zero on successful completion.

    Raises
    ------
    ValueError
        Raised if the combination of switches is inconsistent.
    """
    if field is False:
        print("Running no-field study")
        main_no_field(num_runs)
        return 0

    if field is True and field_study is False:
        print("Running time-dependent field study")
        main_time(num_runs, n_counts)
        return 0

    if field is True and field_study is True:
        print("Running field-scan study")
        main_field(num_runs, field_max, n_field_values)
        return 0

    raise ValueError("Invalid combination of field and field_study.")


if __name__ == "__main__":
    main()