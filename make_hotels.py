import numpy as np
from sklearn.cluster import KMeans
import argparse
import os
import matplotlib.pyplot as plt


def main(vertex_data_fp, H, hotel_fp, algo='kmeans', random_state=0):
    """
    This script strategically places H hotels using K-Means clustering based on input data of vertex locations (scores are not considered).

    Args:
        vertex_data_fp (str): Absolute file path to the file containing the vertex data (i.e. list of vertices and their scores).
        H (int): Number of hotels to place strategically.
        hotel_fp (str): Absolute file path where the results, i.e. strategic hotel locations are to be written.
        algo (str): Specify the algorithm to be implemented for generating the hotels. Can be either 'kmeans' or 'random_uniform'.
        random_state (int or None): Integer to control the seeding to the kmeans algorithm. When an integer value is supplied, it allows for deterministic results from the k-means.
    """

    def read_data_file(input_fp):
        """
        Reads input data file and extracts vertex locations and scores.
        Assumes the file contains 'x y Si' data where:
            x  = x coordinate
            y  = y coordinate
            Si = score

        Args:
            input_fp (str): Absolute path to the input file.

        Returns:
            location (list): List of (x, y) coordinates.
            Si (list): List of corresponding scores.
        """
        if not os.path.exists(input_fp):
            raise FileNotFoundError(f"Input file '{input_fp}' not found.")

        location = []
        Si = []

        with open(input_fp, 'r') as file:
            next(file)  # Skip the first line (header)
            for line in file:
                values = list(map(float, line.split()))
                x, y, score = values[-3], values[-2], values[-1]
                location.append([x, y])
                Si.append(score)

        return location, Si

    # read in the input data 
    c_pos, Si = read_data_file(vertex_data_fp)
    c_pos = np.array(c_pos)
    Si = np.array(Si)
    
    if algo == 'kmeans':
        # Apply K-means clustering with only the locations
        kmeans = KMeans(n_clusters=H, random_state=random_state).fit(c_pos)

        # Get the cluster centers (strategic hotel locations)
        strategic_points = kmeans.cluster_centers_

        # Prepare data for saving: hotel locations and zeroed Si scores
        hotel_Si = np.zeros(H,)
        
    elif algo == 'random_uniform':
        # Get the minimum and maximum values for each dimension (x and y)
        x_min, y_min = np.min(c_pos, axis=0)
        x_max, y_max = np.max(c_pos, axis=0)

        # Generate H random points within the bounds
        strategic_points = np.random.uniform(low=[x_min, y_min], high=[x_max, y_max], size=(H, 2))

        # Prepare data for saving: hotel locations and zeroed Si scores
        hotel_Si = np.zeros(H,)
    
    else: 
        raise RuntimeError('Unknown algorithm specified for generating hotels.')

    combined_array = np.column_stack((strategic_points, hotel_Si))

    # Save the strategic points to a tab-separated file
    np.savetxt(hotel_fp, combined_array, delimiter='\t', header='x\ty\tSi', comments='')

    print(f'{algo} was selected for the hotel generation. If kmeans, the random_state variable is {random_state}')

    # Plot the original points
    plt.scatter(c_pos[:, 0], c_pos[:, 1], c='blue', label='Original Points')

    # Plotting
    to_plt = False
    if to_plt is True:
        plt.scatter(c_pos[:, 0], c_pos[:, 1], c='blue', label='Original Points')
        plt.scatter(strategic_points[:, 0], strategic_points[:, 1], c='red', marker='x', label='Strategic Points')
        plt.legend()
        plt.show()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Find strategic hotel locations using K-Means clustering.')

    parser.add_argument('vertex_data_fp', type=str, 
                        help='Absolute file path to the file containing the vertex data (i.e. list of vertices and their scores).')
    parser.add_argument('H', type=int, 
                        help='Number of hotels to place strategically.')
    parser.add_argument('hotel_fp', type=str, 
                        help='Absolute file path where the results, i.e. strategic hotel locations are to be written.')

    args = parser.parse_args()
    main( args.vertex_data_fp, args.H, args.hotel_fp)