import pickle
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pyproj import Proj, transform
import geopandas as gpd
from shapely.geometry import box

import os
import sys

import argparse

current_dir = os.path.dirname(os.path.abspath(__file__))

# Define the Lambert Conformal Conic projection parameters for Texas
lambert_conformal_conic = Proj(
    proj='lcc',
    lat_1=33.0,  # First standard parallel (near the northern boundary of Texas)
    lat_2=27.0,  # Second standard parallel (near the southern boundary of Texas)
    lat_0=31.0,  # Latitude of origin (central latitude of Texas)
    lon_0=-100.0,  # Central meridian (central longitude of Texas)
    x_0=0,
    y_0=0,
    datum='WGS84'
)
x_min = -328108.7153171253 # run the process_input_data_verX script for these values
y_min = 6603.913017557618


# Inverse transformation function
def inverse_project_coords(x, y):
    lon, lat = lambert_conformal_conic(x, y, inverse=True)
    return (lon, lat)

'''
# Load the river map shapefile
river_map = gpd.read_file(os.path.join(current_dir, 'river_data/pfaf_07_riv3sMERIT_sort/pfaf_07_riv_3sMERIT_sort.shp'))

# Save the GeoDataFrame to a pickle file
with open('river_map.pkl', 'wb') as f:
    pickle.dump(river_map, f)
'''
# Load the GeoDataFrame from the pickle file
with open('river_map.pkl', 'rb') as f:
    river_map = pickle.load(f)
# Define the bounding box for the region of interest (xmin, ymin, xmax, ymax)
bbox = box(minx=-103.5, miny=31, maxx=-103.2, maxy=31.25)
# Clip the river map to the bounding box
river_map = river_map[river_map.intersects(bbox)]

def main(results_fn, disable_plot=False):
    """ Example usage: 
            python print_results.py Fracking_25_scival_N22_H3_D2_Tmax28800_Tch1800_UAVs5_kch5_kdis1.pkl
    """
    
    results_dir = os.path.join(current_dir, 'results')

    # Load data objects from the file
    results_fl =  os.path.join(results_dir,results_fn)
    #results_fl = 'dataSet.pkl'
    #results_fl = './results/data_objects_SET1 2-3_64-45-2-3.pkl'
    with open(results_fl, 'rb') as file:
        H = pickle.load(file)
        N = pickle.load(file)
        D = pickle.load(file)
        T_Max = pickle.load(file)
        T_CH = pickle.load(file)
        c_pos = pickle.load(file)
        Si = pickle.load(file)
        score = pickle.load(file)
        transitions = pickle.load(file)
        halt_times = pickle.load(file)
        max_flight_times = pickle.load(file)
        flight_times = pickle.load(file)
        subtour_u = pickle.load(file)
        subtour_u = np.array(subtour_u) - min(np.array(subtour_u))
        optimal_value = pickle.load(file)
        gap = pickle.load(file)
        nconss = pickle.load(file)
        nvars = pickle.load(file)
        optimal_sol = pickle.load(file)
        process_time = pickle.load(file)


    print(f"number of vertices is N: {N}")
    print(f"number of hotels is H: {H}")
    print(f"number for trips is D: {D}")
    print(f"maximum tour time T_Max: {T_Max}")
    print(f"Max flight time on full charge T_CH: {T_CH}")

    print(f"halt times: {halt_times}")
            
    print("max flight times: ", max_flight_times)

    print("flight times: ", flight_times)

    print("u is: ", subtour_u)

    print(f"Optimal solution is found? {optimal_sol}")

    print("===========================")

    print("Gap:", gap)

    print("Objective value:", optimal_value)

    print("Total tour time in minutes: ", 1/60.0 * (sum(flight_times)+sum(halt_times)))

    print(f"Total charging time: {1/60.0 * sum(halt_times)} minutes")

    print(f"Process time {process_time/ 60.0} minutes")

    print("===========================")
    
    print("# constraints:", nconss)
    
    print("# Variables:", nvars)

    

    ### Plot the solution ###
    # Initialize the figure
    fig, ax = plt.subplots(figsize=(12, 8))

    # Plot the river map
    river_map.plot(ax=ax, color='blue', linewidth=0.5, alpha=0.7)

    # adjust the c_pos coordinates (re-center) so that the reverse projection shall work properly
    # Add the constants to each coordinate pair
    c_pos = [[x + x_min, y + y_min] for x, y in c_pos]

    # Separate the first H_count 'H' coordinates from the rest
    hotels = c_pos[:H]
    #print('hotels: ', hotels)
    vertices = c_pos[H:]
    #print('vertices: ', vertices)

    # Extract the 'x' and 'y' values for both 'H' and the rest
    x_h, y_h = zip(*hotels)
    x_n, y_n = zip(*vertices)

    Si = np.array(Si)
    # Create the scatter plot with different marker styles and colors
    # Apply the inverse transformation
    (lon_h, lat_h) = inverse_project_coords(x_h, y_h)
    plt.scatter(lon_h, lat_h, color='orange', s=50, marker='o', label='Hotel') 
    (lon_n, lat_n) = inverse_project_coords(x_n, y_n)
    plt.scatter(lon_n, lat_n, color='black', s=30*Si[H:], marker='s', label='Vertex')
    plt.xlabel("Longitude [degrees]", fontsize=14)
    plt.ylabel("Latitude [degrees]", fontsize=14)
    # Set x and y-axis tick label size
    plt.tick_params(axis='both', which='major', labelsize=10)
    plt.tick_params(axis='both', which='minor', labelsize=10)


    # Draw lines connecting cities
    colors = [
        'black', 'red', 'green', 'blue', 'brown', 'cyan', 'magenta', 'orange', 'purple', 
        'yellow', 'pink', 'lime', 'teal', 'navy', 'maroon', 'olive', 'grey', 'gold', 'silver' ]# corresponding to different trips
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
                    
                    # Apply the inverse transformation
                    (lon1, lat1) = inverse_project_coords(x1, y1)
                    (lon2, lat2) = inverse_project_coords(x2, y2)
                    plt.plot([lon1, lon2], [lat1, lat2], color=colors[d])
                    
                    # Adding labels to the points.
                    plt.text(lon1, lat1, f'{i}', fontsize=12, ha='right')
                    plt.text(lat2, lat2, f'{j}', fontsize=12, ha='right')


    # Add a colorbar
    #plt.colorbar()

    sys.stdout.flush()  # Ensure the print statements are flushed to the console


    # Show the plot
    if disable_plot is False:
        plt.show()
        


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Display results')

    parser.add_argument('results_fn', type=str, help='results pickle object name, e.g., "Fracking_25_scival_N24_H1_D10_Tmax28800_Tch1800.pkl"')
    parser.add_argument('--disable_plot', action='store_true', help='Enable plotting of diagram')
                       
    args = parser.parse_args()
    
    main(args.results_fn, args.disable_plot)