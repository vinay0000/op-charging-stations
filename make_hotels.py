import numpy as np
from sklearn.cluster import KMeans
import argparse
import os
import matplotlib.pyplot as plt


def main(sci_fn, H, hotel_fn):
    """
    This script strategically places H hotels using K-Means clustering based on 
    input data of vertex locations and scores.

    Args:
        sci_fn (str): Input science data filename (e.g., 'Fracking_25_scival.opc'). File must be in the data/science_data directory.
        H (int): Number of hotels to place strategically.
        hotel_fn (str): Output filename for hotel locations (e.g., 'strategic_hotels.txt'). File is written in the "data/science_data" directory.
    """

    # Define filepaths
    current_dir = os.path.dirname(os.path.abspath(__file__))
    input_data_dir = os.path.join(current_dir,'data/science_data')
    sci_fp = os.path.join(input_data_dir, sci_fn)
    hotel_fp = os.path.join(input_data_dir, hotel_fn)

    def read_data_file(input_fp):
        """
        Reads input data file and extracts vertex locations and scores.
        Assumes the file contains 'x y Si' data where:
            x  = x coordinate
            y  = y coordinate
            Si = score

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
    c_pos, Si = read_data_file(sci_fp)
    c_pos = np.array(c_pos)

    # Apply K-means clustering
    kmeans = KMeans(n_clusters=H, random_state=0).fit(c_pos)

    # Get the cluster centers (strategic hotel locations)
    strategic_points = kmeans.cluster_centers_
    print("Strategically placed points:", strategic_points)

    # Prepare data for saving: hotel locations and zeroed Si scores
    hotel_Si = np.zeros(H,)
    combined_array = np.column_stack((strategic_points, hotel_Si))
    
    # Save the strategic points to a tab-separated file
    np.savetxt(hotel_fp, combined_array, delimiter='\t', header='x\ty\tSi', comments='')

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

    parser.add_argument('sci_fn', type=str, 
                        help='Input science data filename. Must be in the "data/science_data" directory. E.g., "Fracking_25_scival.opc".')
    parser.add_argument('H', type=int, 
                        help='Number of hotels to place strategically.')
    parser.add_argument('hotel_fn', type=str, 
                        help='Output filename for strategic hotel locations. File is written in the "data/science_data" directory.')

    args = parser.parse_args()
    main(args.sci_fn, args.H, args.hotel_fn)