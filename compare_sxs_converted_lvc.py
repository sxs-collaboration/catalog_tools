import argparse
import h5py
import json
import numpy as np
import romspline
import sys
import time

if __name__ == "__main__":
    p = argparse.ArgumentParser(description="Compare two SXS LVC format files")
    p.add_argument("--file1",
                   help="Path to first file",
                   required=True)
    p.add_argument("--file2",
                   help="Path to second file",
                   required=True)
    args = p.parse_args()

def compare_attributes(file1, file2):
    print("Comparing attributes")
    for key in file1.attrs:
        try:
          if file1.attrs[key] == file2.attrs[key]:
              print("= " + key)
          else:
              print("x " + key)
              print("  file1: " + file1.attrs[key])
              print("  file2: " + file2.attrs[key])
        except:
          print("x " + key)
          try:
            print("  file1: " + str(file1.attrs[key]))
          except:
            print("  file1 does not contain key " + key)
          try:
            print("  file2: " + str(file2.attrs[key]))
          except:
            print("  file2 does not contain key " + key)

def compare_splines(filename1, filename2, key):
    spline1 = romspline.ReducedOrderSpline()
    spline1.read(filename1, group=key)
    spline2 = romspline.ReducedOrderSpline()
    spline2.read(filename2, group=key)
    x_min = np.max([spline1.X[0], spline2.X[0]])
    x_max = np.min([spline1.X[-1], spline2.X[-1]])
    x = np.arange(x_min, x_max, np.abs(spline1.X[1]-spline1.X[0]))
    diff = np.max(np.abs(spline1(x) - spline2(x)))
    return diff

if __name__ == "__main__":
    file1 = h5py.File(args.file1)
    file2 = h5py.File(args.file2)

    compare_attributes(file1, file2)

    print("Comparing datasets")
    # Get the keys for each spline
    splines_to_check = []
    time_series_to_check = []
    other_keys = []
    for key in file1:
        try:
            # If the data set is a spline, it should have a key 'deg'
            deg = file1[key]['deg']
            splines_to_check.append(key)
        except:
            try:
                # If the file is a group, it should have its own subkeys
                subkeys = file1[key].keys()
                other_keys.append(key)
            except:
                time_series_to_check.append(key)

    eps = 1.e-15
    for key in splines_to_check:
        diff = compare_splines(args.file1, args.file2, key)
        if np.abs(diff) < eps:
            print("= " + key + " (diff = " + str(diff) + ")")
        else:
            print("x " + key + " (diff = " + str(diff) + ")")

    print("Comparing time series")
    for key in time_series_to_check:
        try:
            diff = np.max(np.abs(np.array(file1[key]) - np.array(file2[key])))
            if diff < eps:
                print("= " + key + " (diff = " + str(diff) + ")")
            else:
                print("x " + key + " (diff = " + str(diff) + ")")
        except:
            print("Cannot diff time series " + key)
            print("File 1: " + str(np.array(file1[key]).shape))
            print("File 2: " + str(np.array(file2[key]).shape))

    file1.close()
    file2.close()

