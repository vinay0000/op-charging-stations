import numpy as np
from sklearn.cluster import KMeans
import argparse

import os

import matplotlib.pyplot as plt

def main(sci_fn, H, hotel_fn):
        """
        """
        #### User inputs ####
        # sci_fn = 'Fracking_25_scival.ophs'
        # H = 5 # Number of Hotels to place strategically
        # 
        ################

        current_dir = os.path.dirname(os.path.abspath(__file__))
        input_data_dir = os.path.join(current_dir,'data/science_data')
        sci_fp = os.path.join(input_data_dir, sci_fn)
        hotel_fp = os.path.join(input_data_dir, hotel_fn)

        def read_data_file(input_fp):
                with open(input_fp, 'r') as file:
                        """
                                location and score:  x y Si
                                x	= x coordinate 
                                y	= y coordinate
                                Si	= score

                                the first line is the starting hotel
                                the second line is the ending hotel
                                the next "H" lines are the extra hotels
                                the remaining lines (N) are the vertices

                        """
                        # Extract location and score data
                        location = []
                        Si = []
                        next(file)  # Skip the first line (header)
                        for line in file:
                                values = list(map(float, line.split()))
                                x, y, _Si = values[-3], values[-2], values[-1]  # Assuming the first two values are x and y, and third is Si
                                location.append([x, y])
                                Si.append(_Si)

                return location, Si

        # read in the input data 
        c_pos, Si = read_data_file(sci_fp)
        c_pos = np.array(c_pos)

        # Apply K-means clustering
        kmeans = KMeans(n_clusters=H, random_state=0).fit(c_pos)

        # Get the cluster centers
        strategic_points = kmeans.cluster_centers_

        print("Strategically placed points:", strategic_points)

        # Combine the arrays
        hotel_Si = np.zeros(H,)
        combined_array = np.column_stack((strategic_points, hotel_Si))
        # Save the combined array to a tab-separated file
        np.savetxt(hotel_fp, combined_array, delimiter='\t', header='x\ty\tSi', comments='')

        # Plot the original points
        plt.scatter(c_pos[:, 0], c_pos[:, 1], c='blue', label='Original Points')

        # Plot the strategically placed points
        plt.scatter(strategic_points[:, 0], strategic_points[:, 1], c='red', label='Strategic Points', marker='x')

        plt.legend()
        #plt.show()

if __name__ == "__main__":
        parser = argparse.ArgumentParser(description='Find strategic hotel locations with k-means.')

        parser.add_argument('sci_fn', type=str, help='Input science data filename. File must be present in the data/science_data directory with .ophs extension. e.g., Fracking_25_scival.ophs') 
        parser.add_argument('H', type=int, help='Number of hotels')
        parser.add_argument('hotel_fn', type=str, help='File name where the resulting hotel data is stored. File is written in the data/science_data directory.')

        args = parser.parse_args()
        
        main(args.sci_fn, args.H, args.hotel_fn)