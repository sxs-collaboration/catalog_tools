import argparse
import h5py
import json
import numpy as np
import sys
import time

# Update next line with path to romspline, if not in a standard location
sys.path.append("/Users/geoffrey/Codes/CatalogAnalysis")
import romspline

def sxs_id_from_alt_names(alt_names):
    """Takes an array of alternative names from an SXS metadata.json file
    and returns the SXS ID of the simulation."""
    pattern = 'SXS'
    sxs_id = str(next((ss for ss in alt_names if pattern in ss), None))
    return sxs_id

def compare_attribute(lvc_key, lvc_value, sxs_value):
    if (lvc_value == sxs_value):
        print("= " + lvc_key + " (" + str(lvc_value) + ")")
    else:
        print("x " + lvc_key + " (lvc: " + str(lvc_value) + ", sxs: " \
                   + str(sxs_value) + ")")

def compare_attributes(lvc, metadata):
    """Compares LVC attributes to SXS metadata"""
    print("Comparing metadata")
    sxs_id = sxs_id_from_alt_names(metadata["alternative_names"])
    compare_attribute("name", lvc.attrs["name"], sxs_id)

    mass1 =  metadata["relaxed_mass1"]
    mass2 =  metadata["relaxed_mass2"]
    compare_attribute("mass1", lvc.attrs["mass1"], mass1)
    compare_attribute("mass2", lvc.attrs["mass2"], mass2)

    eta = (mass1 * mass2) / np.square(mass1 + mass2)
    compare_attribute("eta", lvc.attrs["eta"], eta)

    compare_attribute("spin1x", lvc.attrs["spin1x"],
                      metadata["relaxed_dimensionless_spin1"][0])
    compare_attribute("spin1y", lvc.attrs["spin1y"],
                      metadata["relaxed_dimensionless_spin1"][1])
    compare_attribute("spin1z", lvc.attrs["spin1z"],
                      metadata["relaxed_dimensionless_spin1"][2])
    compare_attribute("spin2x", lvc.attrs["spin2x"],
                      metadata["relaxed_dimensionless_spin2"][0])
    compare_attribute("spin2y", lvc.attrs["spin2y"],
                      metadata["relaxed_dimensionless_spin2"][1])
    compare_attribute("spin2z", lvc.attrs["spin2z"],
                      metadata["relaxed_dimensionless_spin2"][2])

    omega_orbit_vec = metadata["relaxed_orbital_frequency"]
    omega_orbit = np.linalg.norm(omega_orbit_vec)
    compare_attribute("Omega", lvc.attrs["Omega"], omega_orbit)

    # Compute 1 solar mass * G/c^3
    # Use the same value as in Patricia Schmidt's conversion script,
    # attributed to lal
    msun_seconds = 4.925491025543576e-06
    f_lower_at_1_msun = 2.0 * omega_orbit / (2.0 * np.pi * msun_seconds)
    compare_attribute("f_lower_at_1MSUN", lvc.attrs["f_lower_at_1MSUN"],
                      f_lower_at_1_msun)

    eccentricity = metadata['relaxed_eccentricity']
    if str(eccentricity)[0] == "<":
        eccentricity = float(str(eccentricity)[1:])
    if str(eccentricity)[0] == ">":
        eccentricity = float(str(eccentricity)[1:])
    elif str(eccentricity) == '[simulation too short]' \
                         or str(eccentricity) == '[unknown]':
        eccentricity = -44444444.44444444
    else:
        eccentricity = float(eccentricity)
    compare_attribute("eccentricity", lvc.attrs["eccentricity"], eccentricity)

    mean_anomaly = metadata['relaxed_mean_anomaly']
    if isinstance(mean_anomaly, str):
        if mean_anomaly == '[unknown]':
            mean_anomaly = -1.0
        else:
            mean_anomaly = float(mean_anomaly)
    compare_attribute("mean_anomaly", lvc.attrs["mean_anomaly"], mean_anomaly)

    compare_attribute("Omega", lvc.attrs["Omega"], omega_orbit)


def compare_splines(lvc_filename, rhOverM, horizons, key,
                    extrap_order="Extrapolated_N2"):
    spline = romspline.ReducedOrderSpline()
    spline.read(lvc_filename, group=key)

    lvc_file = h5py.File(lvc_filename, 'r')
    start_time = lvc_file.attrs["NR_start_time"]
    peak_time = lvc_file.attrs["NR_peak_time"]

    key_split = key.split('_')
    if (key_split[0] == "amp"):
        sxs_key = "Y_" + key_split[1] + "_" + key_split[2] + ".dat"
        hlm = rhOverM[extrap_order + ".dir"][sxs_key]
        times_raw = hlm[:, 0]
        start_h = np.abs(times_raw - start_time).argmin() + 1
        h = hlm[start_h:, 1] + 1j * hlm[start_h:, 2]
        sxs_values = np.abs(h)
        sxs_times = hlm[start_h:, 0] - peak_time
    elif (key_split[0] == "phase"):
        sxs_key = "Y_" + key_split[1] + "_" + key_split[2] + ".dat"
        hlm = rhOverM[extrap_order + ".dir"][sxs_key]
        times_raw = hlm[:, 0]
        start_h = np.abs(times_raw - start_time).argmin() + 1
        h = hlm[start_h:, 1] + 1j * hlm[start_h:, 2]
        sxs_values = np.unwrap(np.angle(h))
        sxs_times = hlm[start_h:, 0] - peak_time
    else:
        if "remnant" in key:
            ah_key = "AhC.dir"
        elif "2" in key:
            ah_key = "AhB.dir"
        elif "1" in key:
            ah_key = "AhA.dir"

        if "mass" in key:
            quantity_key = "ChristodoulouMass"
            component_key = 1
        elif "spin" in key:
            quantity_key = "chiInertial"
            if 'x' in key:
                component_key = 1
            elif 'y' in key:
                component_key = 2
            elif 'z' in key:
                component_key = 3
        elif "position" in key:
            quantity_key = "CoordCenterInertial"
            if 'x' in key:
                component_key = 1
            elif 'y' in key:
                component_key = 2
            elif 'z' in key:
                component_key = 3
        quantity = horizons[ah_key][quantity_key + ".dat"]
        times_raw_ah = quantity[:, 0]
        start_ah = np.argmin(np.abs(times_raw_ah - start_time))
        sxs_values = quantity[start_ah:, component_key]
        sxs_times = quantity[start_ah:, 0] - peak_time

    # Linf norm of difference
    diff = np.max(np.abs(spline(sxs_times) - sxs_values))
    lvc_file.close()
    return diff


def compare_time_series(lvc, rhOverM, horizons, key, extrap="Extrapolated_N2"):
    lvc_times = np.array(lvc[key])
    start_time = lvc.attrs["NR_start_time"]
    peak_time = lvc.attrs["NR_peak_time"]
    if "HorizonA" in key:
        times_raw_ah = horizons["AhA.dir"]["ArealMass.dat"][:,0]
        start_ah = np.argmin(np.abs(times_raw_ah - start_time))
        sxs_times = times_raw_ah[start_ah:] - peak_time
    elif "HorizonB" in key:
        times_raw_ah = horizons["AhB.dir"]["ArealMass.dat"][:,0]
        start_ah = np.argmin(np.abs(times_raw_ah - start_time))
        sxs_times = times_raw_ah[start_ah:] - peak_time
    elif "CommonHorizon" in key:
        times_raw_ah = horizons["AhC.dir"]["ArealMass.dat"][:,0]
        start_ah = np.argmin(np.abs(times_raw_ah - start_time))
        sxs_times = times_raw_ah[start_ah:] - peak_time

    # Linf norm of difference
    diff = np.max(np.abs(lvc_times - sxs_times))
    return diff

def compare_wave_time_series(lvc, rhOverM, extrap="Extrapolated_N2"):
    lvc_wave_times = np.array(lvc["NRtimes"])
    start_time = lvc.attrs["NR_start_time"]
    peak_time = lvc.attrs["NR_peak_time"]
    # Ensure that all waveform modes were generated using the same times
    diff_all_times_identical = 0.0
    for times in lvc_wave_times:
        diff_all_times_identical += np.max(np.abs(times - lvc_wave_times[0]))
    if diff_all_times_identical > 0.0:
        print("x NRtimes not all identical")
    else:
        print("= NRtimes are all identical")

    # Check the first time series is as expected
    sxs_key = "Y_l2_m2.dat"
    hlm = rhOverM[extrap + ".dir"][sxs_key]
    times_raw = hlm[:, 0]
    start_h = np.abs(times_raw - start_time).argmin() + 1
    sxs_times = hlm[start_h:, 0] - peak_time

    # FIX ME: Why do the lvc NRtimes have one extra point?
    diff = np.max(np.abs(sxs_times - lvc_wave_times[0][1:]))
    if (diff > 0.0):
        print("x NRtimes differs from SXS l=m=2 times (diff = " \
              + str(diff) + ")")
    else:
        print("= NRtimes agrees with SXS l=m=2 times (diff = " \
              + str(diff) + ")")

if __name__ == "__main__":
    p = argparse.ArgumentParser(description="Compare SXS, LVC data")
    p.add_argument("--lvc_file",
                   help="Path to SXS_BBH_????_Res?.h5",
                   required=True)
    p.add_argument("--sxs_waveform",
                   help="Path to rhOverM_Asymptotic_GeometricUnits_Com.h5",
                   required=True)
    p.add_argument("--sxs_horizons",
                   help="Path to Horizons.h5",
                   required=True)
    p.add_argument("--sxs_json",
                   help="Path to metadata.json",
                   required=True)
    args = p.parse_args()
    lvc = h5py.File(args.lvc_file, 'r')
    rhOverM = h5py.File(args.sxs_waveform, 'r')
    horizons = h5py.File(args.sxs_horizons, 'r')
    with open(args.sxs_json, 'r') as file:
        metadata = json.load(file)

    compare_attributes(lvc, metadata)

    print("Comparing datasets")
    # Get the keys for each spline
    splines_to_check = []
    time_series_to_check = []
    other_keys = []
    for key in lvc:
        try:
            # If the data set is a spline, it should have a key 'deg'
            deg = lvc[key]['deg']
            splines_to_check.append(key)
        except:
            try:
                # If the file is a group, it should have its own subkeys
                subkeys = lvc[key].keys()
                other_keys.append(key)
            except:
                time_series_to_check.append(key)

    for key in splines_to_check:
        if "remnant" in key or "1" in key or "2" in key:
            diff = compare_splines(args.lvc_file, rhOverM, horizons, key)
            tol = float(np.array(lvc[key]['tol']))
            if np.abs(diff) < tol:
                print("= " + key + " (diff = " + str(diff) + ", tol = " \
                           + str(tol) + ")")
            else:
                print("x " + key + " (diff = " + str(diff)  + ", tol = " \
                           + str(tol) + ")")

    print("Comparing time series")
    eps = 1.e-15
    for key in time_series_to_check:
        if "Common" in key or "HorizonA" in key or "HorizonB" in key:
            try:
                diff = compare_time_series(lvc, rhOverM, horizons, key)
                if diff < eps:
                    print("= " + key + " (diff = " + str(diff) + ")")
                else:
                    print("x " + key + " (diff = " + str(diff) + ")")
            except:
                print("Cannot diff time series " + key)

    compare_wave_time_series(lvc, rhOverM)

    lvc.close()
    rhOverM.close()
    horizons.close()
    file.close()
