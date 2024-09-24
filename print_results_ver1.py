import pickle
import matplotlib.pyplot as plt
import numpy as np
import os
import sys

import argparse

def main(results_fn, disable_plot=False):
    """ Example usage: 
            python print_results.py Fracking_25_scival_N22_H3_D2_Tmax28800_Tch1800_UAVs5_kch5_kdis1.pkl
    """
    current_dir = os.path.dirname(os.path.abspath(__file__))
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

    print("halt times: ", halt_times)
            
    print("max flight times: ", max_flight_times)

    print("flight times: ", flight_times)

    print("Total tour time: ", sum(flight_times)+sum(halt_times))

    print("u is: ", subtour_u)

    print("Objective value:", optimal_value)

    print(f"Optimal solution is found? {optimal_sol}")

    print("Gap:", gap)
    
    print("# constraints:", nconss)
    
    print("# Variables:", nvars)

    print("Process time:", process_time)

    ### Plot the solution ###
    # Initialize the figure
    plt.figure(figsize=(8, 6))

    # Separate the first H_count 'H' coordinates from the rest
    hotels = c_pos[:H]
    #print('hotels: ', hotels)
    vertices = c_pos[H:]
    #print('vertices: ', vertices)

    # Extract the 'x' and 'y' values for both 'H' and the rest
    x_h, y_h = zip(*hotels)
    x_v, y_v = zip(*vertices)

    Si = np.array(Si)
    # Create the scatter plot with different marker styles and colors
    plt.scatter(x_h, y_h, color='orange', s=50, marker='o', label='Hotel') 
    plt.scatter(x_v, y_v, color='black', s=30*Si[H:], marker='s', label='Vertex')

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
                    plt.plot([x1, x2], [y1, y2], color=colors[d])
                    # Adding labels to the points
                    plt.text(y1, x1, f'{i}', fontsize=8, ha='right')
                    plt.text(y2, x2, f'{j}', fontsize=8, ha='right')

    # Add a colorbar
    plt.colorbar()

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