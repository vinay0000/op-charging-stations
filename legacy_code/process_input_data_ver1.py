""" Script to produce data files for the MIP process, from csv files given by Molly/ George.
"""

import pandas as pd

from pyproj import Proj
import matplotlib.pyplot as plt

import os
import argparse


def main(sci_fn, N, H, D, Tmax, T_CH):
    """ Example usage: 
            python process_input_data.py Fracking_25_scival.csv 24 1 10 28800 1800
        
        Input file must be in the 'data/science_data' directory
        
        The resultant file is stored in the 'data/science_data' directory with the same filename, and extension .ophs
    """

    #### User inputs ####
    #sci_fn = 'Fracking_25_scival.csv' # file must be in the data/science_data directory
    #N = 24 # number of vertices + 2
    #H = 1 # number of extra hotels
    #D = 10 # number of trips
    #Tmax = 28800 # [seconds] = 8 hrs total tour length 
    #T_CH = 1800 # [seconds] = 30 min maximum flight time on full-charge 

    #### Define filepaths ####

    result_fn = sci_fn.replace(".csv", "") + '_N' + str(N) + '_H' + str(H) + '_D' + str(D) + '_Tmax' + str(Tmax) + '_Tch' + str(T_CH) + '.ophs'


    current_dir = os.path.dirname(os.path.abspath(__file__))
    sci_file_dir= 'data/science_data'

    sci_fp = os.path.join(sci_file_dir, sci_fn)
    result_fp = os.path.join(sci_file_dir, result_fn)

    #### Read the CSV file into a pandas DataFrame ####
    df = pd.read_csv(os.path.join(current_dir, sci_fp))


    #### Project ####
    # Define the Lambert Conformal Conic projection parameters for Texas
    # TODO: Provide LCC parameters more appropriate to the input dataset
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

    # Function to project latitude and longitude to x, y
    def project_coords(row):
        x, y = lambert_conformal_conic(row['long'], row['lat'])
        return pd.Series({'x': x, 'y': y})

    # Apply the projection function to each row in the DataFrame
    df[['x', 'y']] = df.apply(project_coords, axis=1)
    #  adjust the coordinates so that the minimum x and y values start at zero.
    df['x'] = df['x'] - df['x'].min()
    df['y'] = df['y'] - df['y'].min()
    
    #### Save the results in the "ophs" format ####

    # Select specific columns from the DataFrame
    selected_columns = df[['x', 'y', 'yearly_avg']]

    headers = [[str(N) + '\t' + str(H) + '\t' + str(D)], [str(Tmax)], [str(T_CH)], []]

    # Write the headers to the file
    with open(result_fp, 'w') as f:
        for header_line in headers:
            f.write('\t'.join(header_line) + '\n')

    # Write the selected columns to a CSV file with custom headers
    selected_columns.to_csv(result_fp, header=False, mode='a', sep='\t', index=False)

    with open(result_fp, 'a') as f:
        f.write('---------------------------------------------------')

    #### Plotting the points weighted by the scores ####
    plt.figure(figsize=(10, 8))
    plt.scatter(df['x'], df['y'], c='blue', marker='o', edgecolor='k', label='Projected Points')
    plt.xlabel('X coordinate')
    plt.ylabel('Y coordinate')
    plt.title('Projected Coordinates (Lambert Conformal Conic)')
    plt.legend()
    plt.grid(True)
    #plt.show()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process user inputs for the tour planning script.')

    parser.add_argument('sci_fn', type=str, help='Scientific filename, File must be in the data/science_data directory. E.g., "Fracking_25_scival.csv"') 
    parser.add_argument('N', type=int, help='Number of vertices + 2')
    parser.add_argument('H', type=int, help='Number of extra hotels')
    parser.add_argument('D', type=int, help='Number of trips')
    parser.add_argument('Tmax', type=float, help='Total tour length in seconds (e.g., 28800 for 8 hours)')
    parser.add_argument('T_CH', type=float, help='Maximum flight time on full-charge in seconds (e.g., 1800 for 30 minutes)')

    args = parser.parse_args()
    
    main(args.sci_fn, args.N, args.H, args.D, args.Tmax, args.T_CH)