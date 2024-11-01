import pickle
import matplotlib.pyplot as plt
from pyproj import Proj
import geopandas as gpd
from shapely.geometry import box
import numpy as np
import os
import sys
import argparse

import utils

current_dir = os.path.dirname(os.path.abspath(__file__))

def print_results(H, N, D, T_Max, T_CH, halt_times, max_flight_times, 
                  flight_times, c_pos, transitions, subtour_u, optimal_sol, gap, optimal_value, 
                  process_time, nconss, nvars):
    """
    Print detailed results about the solution.
    """
    print(f"Number of vertices: {N}")
    print(f"Number of hotels: {H}")
    print(f"Number of trips: {D}")
    print(f"Maximum tour time (T_Max): {T_Max} seconds")
    print(f"Max flight time on full charge (T_CH): {T_CH} seconds")

    print(f"Halt times: {halt_times}")
    print("Max flight times: ", max_flight_times)
    print("Flight times: ", flight_times)
    print("Subtour variables (u): ", subtour_u)

    print(f"Optimal solution found? {optimal_sol}")
    print("===========================")
    print("Gap: ", gap)
    print("Objective value: ", optimal_value)
    print("Total tour time in minutes: ", 1 / 60.0 * (sum(flight_times) + sum(halt_times)))
    print(f"Total charging time: {1 / 60.0 * sum(halt_times)} minutes")
    print(f"Process time: {process_time / 60.0} minutes")
    print("===========================")
    print("# Constraints: ", nconss)
    print("# Variables: ", nvars)

    print("Trip transition details:")
    for d in range(D):
        print("Trip #", d)
        for i in range(H+N):
            for j in range(H+N):
                if transitions[i][j][d] == 1:
             
                    if(i>=H+1 and j>=H+1):
                        print(f"({i}, {j}) <-> {subtour_u[i-H]}, {subtour_u[j-H]})")
                    else:
                        print(f"({i}, {j})")
                    [x1, y1] = c_pos[i]
                    [x2, y2] = c_pos[j]

    sys.stdout.flush()  # Ensure the print statements are flushed to the console

def main(mip_results_fp, cartesian_plot=False, map_plot=False):
    """
    Main function to load results, print them, and plot (Cartesian or geographic coordinates) if directed.
    
    Plotting on geographic coordinates with a river map in the background requires the `river_data/river_map.pkl` file 
    or corresponding shape files.
    
    Parameters:
    -----------
    mip_results_fp : str
        Absolute filepath of the results pickle file.
    disable_plot : bool
        If True, plotting is disabled.
    """
    (H, N, D, T_Max, T_CH, c_pos, Si, score, transitions, halt_times, 
     max_flight_times, flight_times, subtour_u, optimal_value, gap, nconss, 
     nvars, optimal_sol, process_time) = utils.load_results(mip_results_fp)

    # Print results
    print_results(H, N, D, T_Max, T_CH, halt_times, max_flight_times, flight_times, c_pos,
                  transitions, subtour_u, optimal_sol, gap, optimal_value, process_time, nconss, nvars)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Display results from a pickle object.')

    parser.add_argument('mip_results_fp', type=str, 
                        help='Absolute file path to the pickle file which contains the results of the optimization.')

    args = parser.parse_args()
    
    main(args.mip_results_fp)