# Distributed under the MIT License.
# See LICENSE.txt for details.

# A script to convert binary-black-hole simulations from the SXS catalog
# to the LVC format described in https://arxiv.org/abs/1703.01076
# Based on an earlier script by Patricia Schmidt, with input from
# Geoffrey Lovelace and Alyssa Garcia.

import argparse
import h5py
import json
import numpy as np
import os
import sys
import time

# Update next line with path to romspline, if not in a standard location
sys.path.append("/Users/geoffrey/Codes/CatalogAnalysis")
import romspline

########################################################
# Logging
########################################################
# A global variable holding all of the messages printed.
# Will be saved in auxiliary-info of the conversion script.
history = ""

def log(string):
    global history
    history += string + "\n"
    print(string)

########################################################
# Functions to convert waveform quantities to LVC format
########################################################


def sxs_id_from_alt_names(alt_names):
    """Takes an array of alternative names from an SXS metadata.json file
    and returns the SXS ID of the simulation."""
    pattern = 'SXS'
    if not isinstance(alt_names, (list, tuple)):
        alt_names = [alt_names]
    sxs_id = str(next((ss for ss in alt_names if pattern in ss), None))
    return sxs_id


def first_index_after_time(times, target_time):
    """Returns the index of the first time in a list of times after
    time target_time."""
    return np.abs(times - target_time).argmin() + 1


def first_index_after_relaxation_time(times, metadata, offset=-2):
    """Returns the index of the first time in a list of times after the
    relaxation time, which is given as a key in metadata, i.e.
    metadata['relaxed_measurement_time'], except actually return an 
    index offset earlier."""
    relaxation_time = metadata['relaxed_measurement_time']
    return first_index_after_time(times, relaxation_time) - offset


def waveform_norm_squared(
        sxs_format_waveform,
        extrapolation_order="Extrapolated_N2"):
    """Takes an SXS-format waveform and returns the sum of the squared
    amplitude of each (l,m) mode of the wave."""
    extrap = str(extrapolation_order) + ".dir"
    times = sxs_format_waveform[extrap]['Y_l2_m2.dat'][:, 0]
    sum_amp_squared = 0.0 * times
    for key in sxs_format_waveform[extrap].keys():
        if "Y_" in key:
            for i in range(1, 3):
                sum_amp_squared += np.square(
                    sxs_format_waveform[extrap][key][:, i])
    return sum_amp_squared


def peak_time_from_sxs(
        sxs_format_waveform,
        metadata,
        extrapolation_order='Extrapolated_N2'):
    """Returns the time when the sum of the squared amplitudes of an
    SXS-format waveform is largest. Note: this is not necessarily the time of
    the peak of the l=m=2 mode."""
    extrap = extrapolation_order + ".dir"
    # All modes have the same time, so just look at the l=m=2 mode to get the
    # times
    times = sxs_format_waveform[extrapolation_order +
                                ".dir"]['Y_l2_m2.dat'][:, 0]
    start = first_index_after_relaxation_time(times, metadata)
    sum_amp_squared = waveform_norm_squared(
        sxs_format_waveform, extrapolation_order)
    index_peak = start + sum_amp_squared[start:].argmax()
    return sxs_format_waveform[extrap]['Y_l2_m2.dat'][index_peak][0]


def amp_phase_from_sxs(sxs_format_waveform, metadata, modes,
                       extrapolation_order="Extrapolated_N2"):
    """Returns amplitude and phase for an SXS-format waveform, for a list of
    Ylm modes. If modes='all', return all modes for l=2 through l=8,
    inclusive."""
    extrap = str(extrapolation_order) + ".dir"

    if modes == "all":
        modes = [[l,m] for l in range(2,9) for m in range(-l,l+1)]
        
    log("Modes: " + str(modes))
    amps = []
    phases = []
    times_list = []
    l_max = 0
    # All modes have the same time, so just look at the l=m=2 mode to get                                                                                                                                                                                                            
    # the times                                                                                                                                                                                                                                                                      
    times = sxs_format_waveform[extrap]['Y_l2_m2.dat'][:, 0]
    start = first_index_after_relaxation_time(times, metadata)
    peak = peak_time_from_sxs(
        sxs_format_waveform, metadata, extrapolation_order)
    for mode in modes:
        l = mode[0]
        m = mode[1]
        log("Computing mode: l = " + str(l) + ", m = " + str(m))
        mode = "Y_l" + str(l) + "_m" + str(m) + ".dat"
        hlm = sxs_format_waveform[extrap][mode]

        # CHECK ME: is the + sign correct here?
        h = hlm[start:, 1] + 1j * hlm[start:, 2]

        amps.append(np.abs(h))
        phases.append(np.unwrap(np.angle(h)))
        times_list.append(times[start:] - peak)

        if (l > l_max):
            l_max = l

    return modes, times_list, amps, phases, times[start], peak, l_max


def spline_amp_phase_from_sxs(sxs_format_waveform, metadata, modes,
                              extrapolation_order="Extrapolated_N2"):
    """Returns spline amplitude and phase for an SXS-format waveform, for a
    list of Ylm modes. If modes='all', return all modes for l=2 through l=8,
    inclusive."""
    modes, times, amps, phases, start_time, peak_time, l_max \
        = amp_phase_from_sxs(sxs_format_waveform, metadata, modes,
                             extrapolation_order)
    spline_amps = []
    spline_phases = []
    for i, mode in enumerate(modes):
        log("Computing spline for amplitude of mode " + str(mode))
        spline_amps.append(romspline.ReducedOrderSpline(times[i], amps[i]))
        log("Computing spline for phase of mode " + str(mode))
        spline_phases.append(romspline.ReducedOrderSpline(times[i], phases[i]))
    return modes, times, spline_amps, spline_phases, start_time, peak_time, \
        l_max


def write_splines_to_H5(
        out_filename,
        modes,
        spline_amps,
        spline_phases,
        times):
    """Writes spline amplitudes and phases to an HDF5 file
    named out_filename."""
    log("Writing waveform data to " + str(out_filename))
    out_file = h5py.File(out_filename, 'w')
    for i, mode in enumerate(modes):
        l = mode[0]
        m = mode[1]
        out_group_amp = out_file.create_group('amp_l%d_m%d' % (l, m))
        out_group_phase = out_file.create_group('phase_l%d_m%d' % (l, m))
        spline_amps[i].write(out_group_amp)
        spline_phases[i].write(out_group_phase)
        if (l == 2 and m == 2):
            out_file.create_dataset('NRtimes', data=times[i])
    out_file.close()

########################################################
# Functions to convert horizon quantities to LVC format
########################################################


def prepare_horizon_quantity(sxs_horizon_quantity, start_time, peak_time):
    """Returns times and values of an SXS-format horizon quantity, such as
    AhA.dir/ArealMass.dat. This function first truncates the horizon data,
    including only data after the relaxation time. Then, it shifts the time
    by the same amount as the waveforms (i.e., by the peak time). Then,
    return the truncated/shifted times and truncated values."""
    # First, figure out the correct time series
    times_raw_AH = sxs_horizon_quantity[:, 0]
    start_AH = np.argmin(np.abs(times_raw_AH - start_time))
    times_AH = times_raw_AH[start_AH:] - peak_time

    # Loop over remaining components, truncating each one to match times_AH
    quantity_AH_list = []
    for i in range(1, len(sxs_horizon_quantity[0])):
        quantity_AH_list.append(sxs_horizon_quantity[start_AH:, i])
    quantity_AH = np.array(quantity_AH_list)
    return times_AH, quantity_AH


def spline_horizon_quantity(sxs_horizon_quantity, start_time, peak_time):
    """Prepares sxs_horizon_quantity by passing it to
    prepare_horizon_quantity() and then returns a spline of the result."""
    times_AH, quantity_AH = prepare_horizon_quantity(sxs_horizon_quantity,
                                                     start_time, peak_time)
    spline_AH_list = []
    for i in range(0, len(sxs_horizon_quantity[0]) - 1):
        spline_AH_list.append(
            romspline.ReducedOrderSpline(times_AH, quantity_AH[i]))
    return np.array(spline_AH_list)


def insert_spline(sxs_horizons, spline_dictionary, spline_keys,
                  horizon_key, quantity_key, start_time, peak_time):
    """Inserts a spline for a quantity in Horizons.h5 into a dictionary of
    horizon splines. Note: spline_keys is a vector of key names (should be
    length 1 for scalars, length 3 for vectors), where each key is the name
    of a group that will be written in the LVC file, such as mass1-vs-time or
    spin1x-vs-time; horizon_key is AhA, AhB, or AhC; and quantity_key
    is e.g. ChristodoulouMass, the name of the quantity in Horizons.h5
    to be read and splined."""
    log("Computing " + str(spline_keys))
    ah = str(horizon_key) + ".dir"
    qty = str(quantity_key) + ".dat"
    quantity = sxs_horizons[ah][qty]
    spline = spline_horizon_quantity(quantity, start_time, peak_time)
    for i, spline_key in enumerate(spline_keys):
        spline_dictionary[spline_key] = spline[i]


def derived_horizon_quantities_from_sxs(sxs_horizons, start_time, peak_time):
    """Computes quantities derived from quanities in the SXS-format file
    Horizons.h5. Specifically, returns a tuple containing nhat, a unit vector
    from the secondary black hole to the primary black hole; omega_orbit,
    the orbital frequency; LNhat, a unit vector in the direction of the
    orbital angular momentum; and the horizon times of the primary,
    secondary, and remnant, truncated and shifted."""
    t_A, x_A = prepare_horizon_quantity(
        sxs_horizons['AhA.dir']['CoordCenterInertial.dat'], start_time,
        peak_time)
    t_B, x_B = prepare_horizon_quantity(
        sxs_horizons['AhB.dir']['CoordCenterInertial.dat'], start_time,
        peak_time)
    # This is used only for t_C
    t_C, x_C = prepare_horizon_quantity(
        sxs_horizons['AhC.dir']['CoordCenterInertial.dat'], start_time,
        peak_time)

    # n_vec is a unit vector pointing from the secondary (black hole B) to
    # the primary (black hole A). We use np.linalg.norm() to get the Euclidean
    # magnitude of n_vec.
    x_A = x_A.T
    x_B = x_B.T
    n_vec = x_A-x_B
    n_vec_norm = np.linalg.norm(n_vec,axis=-1)
    n_hat  = n_vec/n_vec_norm[:,None]
    
    # We compute dn_vec/dt to get a velocity vector, by computing differences
    # for dn_vec and dt.
    dn_vec = np.diff(n_vec, axis=0)
    dt = np.diff(t_A)
    dn_vec_dt = dn_vec/dt[:,None]

    # The orbital frequency is the magnitude of n_vec x dn_vec/dt
    # This is just from the Newtonian expression |r x v| = r^2 \omega.
    # Because taking a time derivative reduces the length of the array by 1,
    # drop the last time of n_vec so n_vec and dn_vec_dt have the same number
    # of points.
    r_cross_v = np.array([np.cross(n_vec[i], dn_vec_dt[i])
                          for i in range(len(dn_vec_dt))])
    r_cross_v_norm = np.linalg.norm(r_cross_v,axis=-1)
    omega_orbit = r_cross_v_norm / n_vec_norm[:-1]**2

    # Finally, LNhat is a unit vector in the direction of the orbital
    # angular momentum. That is, it is a unit vector in the direction of
    # r x p, which is the same direction as r x v.
    LN_hat = r_cross_v / r_cross_v_norm[:,None]
              
    # Horizons.h5 stores quantities as functions of time. Append time to the
    # derived quantities.
    n_hat_vs_time = np.c_[t_A,n_hat]
    omega_orbit_vs_time = np.c_[t_A[:-1],omega_orbit]
    LN_hat_vs_time = np.c_[t_A[:-1],LN_hat]
    return n_hat_vs_time, omega_orbit_vs_time, LN_hat_vs_time, t_A, t_B, t_C


def spline_horizon_quantity(sxs_horizon_quantity, start_time, peak_time):
    """Prepares sxs_horizon_quantity by passing it to
    prepare_horizon_quantity() and then returns a spline of the result."""
    times_AH, quantity_AH = prepare_horizon_quantity(sxs_horizon_quantity,
                                                     start_time, peak_time)
    spline_AH_list = []
    for i in range(0, len(sxs_horizon_quantity[0]) - 1):
        spline_AH_list.append(
            romspline.ReducedOrderSpline(times_AH, quantity_AH[i]))
    return np.array(spline_AH_list)


def insert_derived_spline(spline_dictionary, spline_keys, derived_quantity):
    """Inserts a spline into a dictionary of horizon splines for a quantity
    derived from Horizons.h5. Note: spline_keys is a vector of key names
    (should be length 1 for scalars, length 3 for vectors), where each key is
    the name of a group that will be written in the LVC file, such as
    Omega-vs-time or LNhatx-vs-time; and derived_quantity is the quantity to
    be splined. The derived_quantity should be computed using
    erived_horizon_quantities_from_sxs()."""
    # N.B. times already shifted, truncated to remove junk by
    # derived_horizon_quantities_from_sxs()
    log("Computing " + str(spline_keys))
    times_AH = derived_quantity[:, 0]
    quantity_AH = derived_quantity[:, 1:]
    spline_AH_list = []
    for a in range(0, len(derived_quantity[0]) - 1):
        spline_AH_list.append(romspline.ReducedOrderSpline(
            times_AH, quantity_AH[:, a]))
    spline = np.array(spline_AH_list)
    for i, spline_key in enumerate(spline_keys):
        spline_dictionary[spline_key] = spline[i]


def horizon_splines_from_sxs(horizons, start_time, peak_time):
    """Prepares a dictionary of the horizon-quantity splines that the LVC
    format expects, starting with an SXS-format Horizons.h5. start_time and
    peak_time are determined by amp_phase_from_sxs()."""
    horizon_splines = {}

    # Christodoulou mass
    insert_spline(horizons, horizon_splines, ['mass1-vs-time'],
                  'AhA', 'ChristodoulouMass', start_time, peak_time)
    insert_spline(horizons, horizon_splines, ['mass2-vs-time'],
                  'AhB', 'ChristodoulouMass', start_time, peak_time)
    insert_spline(horizons, horizon_splines, ['remnant-mass-vs-time'],
                  'AhC', 'ChristodoulouMass', start_time, peak_time)

    # Dimensionless spin
    insert_spline(horizons, horizon_splines,
                  ['spin1x-vs-time', 'spin1y-vs-time', 'spin1z-vs-time'],
                  'AhA', 'chiInertial', start_time, peak_time)
    insert_spline(horizons, horizon_splines,
                  ['spin2x-vs-time', 'spin2y-vs-time', 'spin2z-vs-time'],
                  'AhB', 'chiInertial', start_time, peak_time)
    insert_spline(horizons, horizon_splines,
                  ['remnant-spinx-vs-time', 'remnant-spiny-vs-time',
                      'remnant-spinz-vs-time'],
                  'AhC', 'chiInertial', start_time, peak_time)

    # Position
    insert_spline(horizons,
                  horizon_splines,
                  ['position1x-vs-time',
                   'position1y-vs-time',
                   'position1z-vs-time'],
                  'AhA',
                  'CoordCenterInertial',
                  start_time,
                  peak_time)
    insert_spline(horizons,
                  horizon_splines,
                  ['position2x-vs-time',
                   'position2y-vs-time',
                   'position2z-vs-time'],
                  'AhB',
                  'CoordCenterInertial',
                  start_time,
                  peak_time)
    insert_spline(horizons, horizon_splines,
                  ['remnant-positionx-vs-time', 'remnant-positiony-vs-time',
                   'remnant-positionz-vs-time'],
                  'AhC', 'CoordCenterInertial', start_time, peak_time)

    # Derived quantities: nhat, omega_orbit, LNhat
    n_hat, omega_orbit, LN_hat, t_A, t_B, t_C \
        = derived_horizon_quantities_from_sxs(horizons, start_time, peak_time)
    insert_derived_spline(
        horizon_splines, [
            'nhatx-vs-time', 'nhaty-vs-time', 'nhatz-vs-time'], n_hat)
    insert_derived_spline(horizon_splines, ['Omega-vs-time'], omega_orbit)
    insert_derived_spline(
        horizon_splines, [
            'LNhatx-vs-time', 'LNhaty-vs-time', 'LNhatz-vs-time'], LN_hat)

    return horizon_splines, t_A, t_B, t_C


def write_horizon_splines_from_sxs(
        out_filename,
        horizon_splines,
        primary_horizon_times,
        secondary_horizon_times,
        remnant_horizon_times):
    """Takes a dictionary of horizon splines, prepared with
    horizon_splines_from_sxs, and writes each spline into an HDF5 file. Also
    outputs the horizon times for the individual and remnant black holes,
    truncated to remove junk radiation and shifted."""
    log("Writing horizon data to " + str(out_filename))
    out_file = h5py.File(out_filename, 'a')
    for key in horizon_splines.keys():
        out_group = out_file.create_group(key)
        horizon_splines[key].write(out_group)
    out_file.create_dataset('HorizonATimes', data=primary_horizon_times)
    out_file.create_dataset('HorizonBTimes', data=secondary_horizon_times)
    out_file.create_dataset('CommonHorizonTimes', data=remnant_horizon_times)
    out_file.close()

########################################################
# Functions to convert metadata to the LVC format
########################################################


def simulation_type_from_spins(dimensionless_spin_1, dimensionless_spin_2):
    """Helper function that classifies a simulation with dimensionless_spin_1
    on the primary (larger) black hole and dimensionless_spin_2 on the
    secondary black hole as nonspinning (if no spin component > 0.01),
    aligned-spin (if only z components are > 0.01),
    or precessing (otherwise)."""
    spin_zero_threshold = 0.01  # treat spins smaller than this as zero
    # Types defined in arXiv:1703.01076
    nonspinning_type = "non-spinning"
    aligned_type = "aligned-spins"
    precessing_type = "precessing"

    simulation_type = nonspinning_type

    if (abs(dimensionless_spin_1[2]) > spin_zero_threshold or
            abs(dimensionless_spin_2[2]) > spin_zero_threshold):
        simulation_type = aligned_type

    if (abs(dimensionless_spin_1[0]) > spin_zero_threshold or
        abs(dimensionless_spin_2[0]) > spin_zero_threshold or
        abs(dimensionless_spin_1[1]) > spin_zero_threshold or
            abs(dimensionless_spin_2[1]) > spin_zero_threshold):
        simulation_type = precessing_type

    return simulation_type


def find_comparable_simulations(sxs_id, catalog, catalog_resolutions):
    """Given an SXD ID sxs_id, a dictionary catalog_data (whose keys
    are SXS ID numbers and whose values are dictionaries including
    masses and spins) and a dictionary catalog_resolutions
    (whose keys are SXS ID numbers and whose values are dictionaries
    containing lists of integers indicating what resolutions are
    available in the catalog for that SXS ID number), return an
    LVC-format filename that is a comparable simulation with
    multiple resolutions available.
    """
    mass1 = catalog[sxs_id]['initial_mass1']
    mass2 = catalog[sxs_id]['initial_mass2']
    spin1 = catalog[sxs_id]['initial_dimensionless_spin1']
    spin2 = catalog[sxs_id]['initial_dimensionless_spin2']
    mass_ratio = mass1 / mass2
    spin1_magnitude = np.linalg.norm(spin1)
    spin2_magnitude = np.linalg.norm(spin2)

    # Select SXS ID numbers with multiple resolutions available
    has_multiple_resolutions = []
    for key in catalog_resolutions:
        if len(catalog_resolutions[key]) > 1:
            has_multiple_resolutions.append(key)

    # Select SXS ID numbers with multiple resolutions and the same
    # initial data type and spec revision
    same_id_and_spec_revision = []
    for key in has_multiple_resolutions:
        if catalog[key]['spec_revisions'] == catalog[sxs_id]['spec_revisions'] \
            and catalog[key]['initial_data_type'] \
            == catalog[sxs_id]['initial_data_type']:
                same_id_and_spec_revision.append(key)

    if len(same_id_and_spec_revision) > 0:
        has_multiple_resolutions = same_id_and_spec_revision

    # Select an SXS ID number with multiple resolutions that is "similar"
    # to sxs_id. First, try to find a case with the same initial data type
    # and SpEC revision. If that is possible, then choose from among the
    # ones with the closest mass ratio and spin magnitudes.
    # Otherwise, just select a simulation with the closest mass ratios
    # and spin magnitudes.
    # Note: as in the previous script, here I use initial masses and spins,
    # not relaxed masses and spins.
    mass_spin_diff_best = 1.e100
    key_best = ""
    for key in has_multiple_resolutions:
        current_mass1 = catalog[key]['initial_mass1']
        current_mass2 = catalog[key]['initial_mass2']
        current_spin1 = catalog[key]['initial_dimensionless_spin1']
        current_spin2 = catalog[key]['initial_dimensionless_spin2']
        current_mass_ratio = current_mass1 / current_mass2
        current_spin1_magnitude = np.linalg.norm(current_spin1)
        current_spin2_magnitude = np.linalg.norm(current_spin2)
        mass_spin_diff = np.abs(mass_ratio - current_mass_ratio) \
                          + np.abs(spin1_magnitude - current_spin1_magnitude) \
                          + np.abs(spin2_magnitude - current_spin2_magnitude)
        if mass_spin_diff < mass_spin_diff_best:
            mass_spin_diff_best = mass_spin_diff
            key_best = key

    resolution_best = np.max(catalog_resolutions[key])
    return key_best.replace(':', '_') + "_Res" + str(resolution_best) + ".h5"

def write_metadata_from_sxs(out_filename, resolution, metadata, catalog,
                            catalog_resolutions, start_time, peak_time, l_max):
    """Writes metadata to an LVC-format file, including both adding
    auxiliary-info/metadata.json and setting attributes conforming to
    the format given by arXiv:1703.01076. Input arguments are the
    output filename out_filename, the resolution of the simulation (an int),
    metadata (read from metadata.json) for this simulation/resolution,
    the start_time, peak time, and l_max determined by amp_phase_from_sxs().
    The argument catalog_resolutions is a dictionary (read from a
    json file) whose keys are SXS ID numbers, and whose values are lists
    of integers corresponding to the resolutions available for that
    SXS ID number in the SXS catalog. The argument catalog can be read
    from e.g. arXiv:1904.04831's ancilliary file sxs_catalog.json; the
    keys are SXS ID numbers, and the values are metadata for the simulation,
    such as masses and spins.
    """
    log("Writing metadata")
    # Get the SXS ID for the current simulation
    names = metadata['alternative_names']
    if not isinstance(names, (list, tuple)):
        names = [names]
    sxs_id = sxs_id_from_alt_names(names)
    out_file = h5py.File(out_filename, 'a')

    # Put the metadata for this simulation into the auxiliary-info group
    aux_group = out_file.create_group('auxiliary-info')
    json_this_sim = json.dumps(metadata)
    aux_group.create_dataset('metadata.json', data=json_this_sim)

    # Note: all numerical quantities for LVC attributes
    # are at relaxation time unless otherwise indicated
    mass1 = metadata["relaxed_mass1"]
    mass2 = metadata["relaxed_mass2"]
    total_mass = mass1 + mass2
    eta = (mass1 * mass2) / np.square(total_mass)

    # Round slightly larger mass ratios to 1
    # CHECK ME: should we do this?
    if eta > 0.25 and eta < 0.2501:
        eta = 0.25

    dimensionless_spin1 = metadata["relaxed_dimensionless_spin1"]
    dimensionless_spin2 = metadata["relaxed_dimensionless_spin2"]

    position1 = np.array(metadata["relaxed_position1"])
    position2 = np.array(metadata["relaxed_position2"])
    separation = position1 - position2

    n_hat = separation / np.linalg.norm(separation)
    omega_orbit_vec = metadata["relaxed_orbital_frequency"]
    omega_orbit = np.linalg.norm(omega_orbit_vec)
    omega_grav_22 = 2.0 * omega_orbit

    # Unit vector in direction of orbital angular momentum, which
    # is the same as the direction of the orbital angular velocity
    # at the relaxation time
    ln_hat = omega_orbit_vec / omega_orbit

    # Compute 1 solar mass * G/c^3
    # Use the same value as in Patricia Schmidt's conversion script,
    # attributed to lal
    msun_seconds = 4.925491025543576e-06

    # Eccentricity quantities
    eccentricity = metadata['relaxed_eccentricity']
    if str(eccentricity)[0] == "<":
        eccentricity = float(str(eccentricity)[1:])
    if str(eccentricity)[0] == ">":
        eccentricity = float(str(eccentricity)[1:])
    elif str(eccentricity) == '[simulation too short]' \
                         or str(eccentricity) == '[unknown]':
        # CHECK ME: is this the right thing to do for cases where we can't
        # measure eccentricity?
        log("Warning: eccentricity not measured for this simulation")
        eccentricity = -1.0
    else:
        eccentricity = float(eccentricity)

    mean_anomaly = metadata['relaxed_mean_anomaly']
    if isinstance(mean_anomaly, str):
        if mean_anomaly == '[unknown]':
            mean_anomaly = -1.0
        else:
            mean_anomaly = float(mean_anomaly)

    # Set metadata attributes of the output file

    out_file.attrs['NR-group'] = "SXS"
    out_file.attrs['name'] = sxs_id
    out_file.attrs['alternative-names'] = ",".join(names)
    out_file.attrs['simulation-type'] \
        = simulation_type_from_spins(dimensionless_spin1, dimensionless_spin2)

    # our metdata uses lowercase for object types (bh, ns), but LVC wants
    # BH or NS, so convert to uppercase with upper()
    out_file.attrs['object1'] = metadata['object1'].upper()
    out_file.attrs['object2'] = metadata['object2'].upper()

    out_file.attrs['PN_approximant'] = "none, NR only"

    out_file.attrs['NR_reference_time'] = float()

    out_file.attrs['mass1'] = mass1
    out_file.attrs['mass2'] = mass2
    out_file.attrs['eta'] = eta

    # Metadata not in arXiv:1703.01076
    out_file.attrs['NR_reference_time'] = metadata['relaxed_measurement_time']
    out_file.attrs['NR_start_time'] = start_time
    out_file.attrs['NR_peak_time'] = peak_time  # not in Patricia's script
    out_file.attrs['NR_frame'] = 'inertial'

    out_file.attrs['f_lower_at_1MSUN'] = omega_grav_22 / \
        (2.0 * np.pi * msun_seconds)
    out_file.attrs['Omega'] = omega_orbit
    out_file.attrs['spin1x'] = dimensionless_spin1[0]
    out_file.attrs['spin1y'] = dimensionless_spin1[1]
    out_file.attrs['spin1z'] = dimensionless_spin1[2]
    out_file.attrs['spin2x'] = dimensionless_spin2[0]
    out_file.attrs['spin2y'] = dimensionless_spin2[1]
    out_file.attrs['spin2z'] = dimensionless_spin2[2]
    out_file.attrs['nhatx'] = n_hat[0]
    out_file.attrs['nhaty'] = n_hat[1]
    out_file.attrs['nhatz'] = n_hat[2]
    out_file.attrs['LNhatx'] = ln_hat[0]
    out_file.attrs['LNhaty'] = ln_hat[1]
    out_file.attrs['LNhatz'] = ln_hat[2]

    out_file.attrs['eccentricity'] = eccentricity
    out_file.attrs['mean_anomaly'] = mean_anomaly

    out_file.attrs['type'] = "NRinjection"
    out_file.attrs['Format'] = 3
    out_file.attrs['NR-code'] = "SpEC"
    out_file.attrs['modification-date'] = time.strftime("%Y-%m-%d")
    out_file.attrs['point-of-contact-email'] = 'questions@black-holes.org'
    # CHECK ME: is this always the right reference to cite?
    out_file.attrs['INSPIRE-bibtex-keys'] = "Boyle:2019kee"

    # CHECK-ME: should this be public?
    out_file.attrs['license'] = "LVC-internal"

    out_file.attrs['NR-techniques'] = "Quasi-Equilibrium-ID, GH, RWZ-h, " \
        + "Extrapolated-Waveform, ApproxKillingVectorSpin, Christodoulou-Mass"
    out_file.attrs['Lmax'] = l_max

    # Resolution-related functions
    # If highest resoultion available, 'production-run' = 1, otherwise
    # 'production-run' = 0. If multiple resolutions available, list their
    # file names for 'files-in-error-series'. Otherwise, list filenames of a
    # comparable simulation with multiple resolutions.
    # This can all be determined from metadata_all_simulations, a json file
    # containing all simulations in the catalog.

    resolutions = catalog_resolutions[sxs_id]
    if len(resolutions) > 1:
        error_name_base = out_filename.split('/')[-1].split("_Res")[0] + "_Res"
        error_series = [error_name_base + str(i) + ".h5" for i in resolutions]
        out_file.attrs['files-in-error-series'] = ",".join(error_series)
        out_file.attrs['comparable-simulation'] = ""
        if (resolution == np.max(resolutions)):
            out_file.attrs['production-run'] = 1
        else:
            out_file.attrs['production-run'] = 0
    else:
        out_file.attrs['production-run'] = 1
        out_file.attrs['files-in-error-series'] = ""
        comparable = find_comparable_simulations(sxs_id, catalog,
                                                 catalog_resolutions)
        out_file.attrs['comparable-simulation'] = comparable

    out_file.close()


def convert_simulation(sxs_data_path, resolution, sxs_catalog_metadata_path,
                       sxs_catalog_resolutions_path, modes, out_path):
    """Convert a simulation from the SXS BBH catalog into the LVC format.
    sxs_data_path is a path to a directory that must contain the following:
        i)   rhOverM_Asymptotic_GeometricUnits_CoM.h5
        ii)  Horizons.h5
        iii) metadata.json
    Additionally, the function requires paths to the following:
        iv)  sxs_catalog_metadata_path points to sxs_catalog.json, which
             contains information such as masses and spins. A file in this
             format is available from arXiv:1904.04831. Keys are
             SXS ID numbers, and values are dictionaries of metadata.
        v)   sxs_catalog_resolutions_path points to a file whose keys are
             SXS ID numbers and whose values are lists of integers, where each
             integer corresponds to a resolution available in the catalog.
    The option resolution is an integer labeling the resolution of the
    converted waveform. Modes is an array of the format
    [[l1, m1], [l2, m2], ...] listing the l,m modes to convert. This function
    outputs a file in LVC format named SXS_BBH_\#\#\#\#_Res\#.h5
    in out_path."""
    horizons = h5py.File(sxs_data_path + "/Horizons.h5", 'r')
    rhOverM = h5py.File(sxs_data_path
                         + "/rhOverM_Asymptotic_GeometricUnits_CoM.h5")
    with open(sxs_data_path + "/metadata.json", 'r') as file:
        metadata = json.load(file)
        file.close()
    with open(sxs_catalog_metadata_path, 'r') as file:
        sxs_catalog_metadata = json.load(file)
        file.close()
    with open(sxs_catalog_resolutions_path, 'r') as file:
        sxs_catalog_resolutions = json.load(file)
        file.close()

    sxs_id = sxs_id_from_alt_names(metadata['alternative_names'])
    log("Converting " + sxs_id)
    log("Running script " + os.path.realpath(__file__) + " in " + os.getcwd())
    log("convert_sxs_to_lvc.py called with the following parameters:")
    log("  sxs_data_path: " + sxs_data_path)
    log("  resolution: " + str(resolution))
    log("  sxs_catalog_metadata_path: " + sxs_catalog_metadata_path)
    log("  sxs_catalog_resolutions_path: " + sxs_catalog_resolutions_path)
    log("  modes: " + str(modes))
    log("  out_path: " + str(out_path))

    out_name = out_path + "/" + sxs_id.replace(':', '_') + "_Res" \
                        + str(resolution) + ".h5"
    log("Output filename is " + out_name)

    modes, times, spline_amps, spline_phases, start_time, peak_time, l_max \
        = spline_amp_phase_from_sxs(rhOverM, metadata, modes)
    write_splines_to_H5(out_name, modes, spline_amps, spline_phases, times)

    horizon_splines_to_write, t_A, t_B, t_C \
        = horizon_splines_from_sxs(horizons, start_time, peak_time)
    write_horizon_splines_from_sxs(out_name, horizon_splines_to_write,
                                   t_A, t_B, t_C)

    write_metadata_from_sxs(out_name, resolution, metadata,
                            sxs_catalog_metadata, sxs_catalog_resolutions,
                            start_time, peak_time, l_max)

    out_file = h5py.File(out_name, 'a')

    # Save a copy of this script in auxiliary-info
    with open(os.path.realpath(__file__), 'r') as file:
        out_file["auxiliary-info"].create_dataset('convert_sxs_to_lvc.py',
                                                   data=file.read())

    if "VersionHist.ver" in list(rhOverM.keys()):
        log("Writing VersionHist.ver")
        out_file["auxiliary-info"].create_dataset('VersionHist.ver', \
                                                data=rhOverM['VersionHist.ver'])
    else:
        log("No VersionHist.ver found. Data being converted is version 0.")
    log("Writing log")
    out_file["auxiliary-info"].create_dataset('ConversionLog.txt',
                                              data=history)
    out_file.close()

    horizons.close()
    rhOverM.close()

if __name__ == "__main__":
    p = argparse.ArgumentParser(description="Convert SXS data to LVC format")
    p.add_argument("--sxs_data",
                   help="Path containing rh*h5, Horizons.h5, metadata.json",
                   required=True)
    p.add_argument("--resolution", type=int,
                   help="Int giving the resolution of the data to convert",
                   required=True)
    p.add_argument("--sxs_catalog_metadata",
                   help="Path to sxs_catalog.json (e.g. via arXiv:1904.04831)",
                   required=True)
    p.add_argument("--sxs_catalog_resolutions",
                   help="Path to sxs_catalog_resolutions.json",
                   required=True)
    p.add_argument("--modes", help="'all', '22only', or int for max l")
    p.add_argument("--out_path", help="Path where LVC format file is output")
    args = p.parse_args()

    if args.modes == 'all':
        modes = [[l,m] for l in range(2,9) for m in range(-l,l+1)]
    elif args.modes == '22only':
        modes = [[2, 2], [2, -2]]
    else:
        l_max = int(args.modes)
        modes = [[l,m] for l in range(2,l_max+1) for m in range(-l,l+1)]
        

    convert_simulation(args.sxs_data, args.resolution, args.sxs_catalog_metadata,
                       args.sxs_catalog_resolutions, modes, args.out_path)
