import pickle
import matplotlib.pyplot as plt
from pyproj import Proj
import geopandas as gpd
from shapely.geometry import box
import numpy as np
import os
import sys
import argparse

current_dir = os.path.dirname(os.path.abspath(__file__))

def load_results(mip_results_fp):
    """
    Load result data from a pickle file.
    
    Parameters:
    -----------
    mip_results_fp : str
        Absolute file path to the pickle file which contains the results of the optimization.
    
    Returns:
    --------
    tuple : All loaded results from the pickle file.
    """
    # Ensure file exists
    if not os.path.exists(mip_results_fp):
        raise FileNotFoundError(f"Results file '{mip_results_fp}' not found.")
    
    with open(mip_results_fp, 'rb') as file:
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

    return (H, N, D, T_Max, T_CH, c_pos, Si, score, transitions, halt_times, 
            max_flight_times, flight_times, subtour_u, optimal_value, gap, 
            nconss, nvars, optimal_sol, process_time)

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


def plot_on_cartesian(H, N, D, c_pos, Si, transitions):
    """
    Plot the solution on a 2D map with hotels and vertices.
    """
    # Initialize the figure
    fig, ax = plt.subplots(figsize=(12, 8))

    # Separate hotels and vertices
    hotels = c_pos[:H]
    vertices = c_pos[H:]

    # Extract x, y coordinates for hotels and vertices
    x_h, y_h = zip(*hotels)
    x_v, y_v = zip(*vertices)

    Si = np.array(Si)

    # Plot hotels and vertices
    plt.scatter(y_h, x_h, color='orange', s=50, marker='o', label='Hotel')
    plt.scatter(y_v, x_v, color='black', s=30 * Si[H:], marker='s', label='Vertex')

    plt.xlabel("X [meters]", fontsize=14)
    plt.ylabel("Y [meters]", fontsize=14)
    plt.tick_params(axis='both', which='major', labelsize=10)

    # Define trip colors
    colors = [
        'black', 'red', 'green', 'blue', 'brown', 'cyan', 'magenta', 'orange', 'purple', 
        'yellow', 'pink', 'lime', 'teal', 'navy', 'maroon', 'olive', 'grey', 'gold', 'silver'
    ]

    # Plot transitions between hotels and vertices for each trip
    for d in range(D):
        print(f"Trip #{d}")
        for i in range(H + N):
            for j in range(H + N):
                if transitions[i][j][d] == 1:
                    [x1, y1] = c_pos[i]
                    [x2, y2] = c_pos[j]
                    plt.plot([y1, y2], [x1, x2], color=colors[d])

                    # Label points
                    plt.text(y1, x1, f'{i}', fontsize=10, ha='right')
                    plt.text(y2, x2, f'{j}', fontsize=10, ha='right')

    plt.legend()
    plt.show()

def plot_on_map(H, N, D, c_pos, Si, transitions):
    """
    Plot results on a map using the Lambert Conformal Conic projection for Texas.
    This function plots the positions of hotels and vertices on a map, along with the 
    connections between them (transitions) for multiple trips.

    Plotting on geographic coordinates with a river map in the background requires the `river_data/river_map.pkl` file 
    or corresponding shape files.s
    """

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

    # Offset values (these need to be pre-calculated for your dataset)
    # run the process_input_data_verX script for these values
    x_min = -328108.7153171253 
    y_min = 6603.913017557618

    # Inverse transformation function
    def inverse_project_coords(x, y):
        lon, lat = lambert_conformal_conic(x, y, inverse=True)
        return (lon, lat)

    '''
    # Below snippet is to be executed if the river_map.pkl file is not present.
    # Load the river map shapefile
    river_map = gpd.read_file(os.path.join(current_dir, 'river_data/pfaf_07_riv3sMERIT_sort/pfaf_07_riv_3sMERIT_sort.shp'))

    # Save the GeoDataFrame to a pickle file
    with open('river_map.pkl', 'wb') as f:
        pickle.dump(river_map, f)
    '''
    # Load the river map (GeoDataFrame) from the pickle file
    with open('river_data/river_map.pkl', 'rb') as f:
        river_map = pickle.load(f)
        
    # Define the bounding box for the region of interest (xmin, ymin, xmax, ymax)
    bbox = box(minx=-103.5, miny=31, maxx=-103.2, maxy=31.25)
    river_map = river_map[river_map.intersects(bbox)] # Clip the river map to the bounding box

    ### Plotting the solution ###
    fig, ax = plt.subplots(figsize=(12, 8))

    # Plot the river map
    river_map.plot(ax=ax, color='blue', linewidth=0.5, alpha=0.7)

    # adjust the c_pos coordinates (re-center) so that the reverse projection works properly.
    c_pos = [[x + x_min, y + y_min] for x, y in c_pos] # Add the constants to each coordinate pair

    # Separate hotel and vertex coordinates
    # Separate the first H_count 'H' coordinates from the rest
    hotels = c_pos[:H]
    #print('hotels: ', hotels)
    vertices = c_pos[H:]
    #print('vertices: ', vertices)

    # Extract the 'x' and 'y' values for both 'H' and the rest
    x_h, y_h = zip(*hotels)
    x_n, y_n = zip(*vertices)

    Si = np.array(Si)
    
    # Convert coordinates to lat/lon for plotting
    # Apply the inverse transformation
    (lon_h, lat_h) = inverse_project_coords(x_h, y_h)
    (lon_n, lat_n) = inverse_project_coords(x_n, y_n)
    
    # Scatter plot for hotels and vertices
    plt.scatter(lon_h, lat_h, color='orange', s=50, marker='o', label='Hotel') 
    plt.scatter(lon_n, lat_n, color='black', s=30*Si[H:], marker='s', label='Vertex')
    
    plt.xlabel("Longitude [degrees]", fontsize=14)
    plt.ylabel("Latitude [degrees]", fontsize=14)
    # Set x and y-axis tick label size
    plt.tick_params(axis='both', which='major', labelsize=10)
    plt.tick_params(axis='both', which='minor', labelsize=10)

    # Draw lines between vertices, hotels according to transitions
    colors = [
        'black', 'red', 'green', 'blue', 'brown', 'cyan', 'magenta', 'orange', 'purple', 
        'yellow', 'pink', 'lime', 'teal', 'navy', 'maroon', 'olive', 'grey', 'gold', 'silver'
    ]  # Colors for different trips

    for d in range(D):
        for i in range(H+N):
            for j in range(H+N):
                if transitions[i][j][d] == 1:
                    [x1, y1] = c_pos[i]
                    [x2, y2] = c_pos[j]
                    # Apply the inverse transformation
                    (lon1, lat1) = inverse_project_coords(x1, y1)
                    (lon2, lat2) = inverse_project_coords(x2, y2)
                    plt.plot([lon1, lon2], [lat1, lat2], color=colors[d])
                    
                    # Adding labels to the points.
                    plt.text(lon1, lat1, f'{i}', fontsize=12, ha='right')
                    plt.text(lat2, lat2, f'{j}', fontsize=12, ha='right')
    # Show the plot                 
    plt.show()

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
     nvars, optimal_sol, process_time) = load_results(mip_results_fp)

    # Print results
    print_results(H, N, D, T_Max, T_CH, halt_times, max_flight_times, flight_times, c_pos,
                  transitions, subtour_u, optimal_sol, gap, optimal_value, process_time, nconss, nvars)

    # Plot the solution
    if cartesian_plot is True:
        plot_on_cartesian(H, N, D, c_pos, Si, transitions)
    
    if map_plot is True:
        plot_on_map(H, N, D, c_pos, Si, transitions)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Display results from a pickle object.')

    parser.add_argument('mip_results_fp', type=str, 
                        help='Absolute file path to the pickle file which contains the results of the optimization.')
    parser.add_argument('--cartesian_plot', action='store_true', help='Plot UAV route on a Cartesian map of rivers.')
    parser.add_argument('--map_plot', action='store_true', help='Plot UAV route on a geographic map of rivers.')

    args = parser.parse_args()
    
    main(args.mip_results_fp, args.cartesian_plot, args.map_plot)