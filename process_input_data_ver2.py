""" Script to produce data files for the MIP process, from csv files given by Molly/ George.
"""

import pandas as pd

from pyproj import Proj
import matplotlib.pyplot as plt

import os
import argparse


def main(sci_raw_fn, result_fn):
    """ Example usage: 
            python process_input_data_ver2.py Fracking_25_scival.csv
        
        Input file must be in the 'data/science_data' directory
        
        The resultant file is stored in the 'data/science_data' directory with the same filename, and extension .ophs
    """

    #### User inputs ####
    #sci_fn = 'Fracking_25_scival.csv' # file must be in the data/science_data directory


    #### Define filepaths ####
    current_dir = os.path.dirname(os.path.abspath(__file__))
    sci_file_dir= 'data/science_data'

    sci_raw_fp = os.path.join(sci_file_dir, sci_raw_fn)
    result_fp = os.path.join(sci_file_dir, result_fn)

    #### Read the CSV file into a pandas DataFrame ####
    df = pd.read_csv(os.path.join(current_dir, sci_raw_fp))

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
    print("df['x'].min(): ", df['x'].min())
    print("df['y'].min(): ", df['y'].min())
    
    df['x'] = df['x'] - df['x'].min()
    df['y'] = df['y'] - df['y'].min()
    
    #### Save the results in the "ophs" format ####

    # Select specific columns from the DataFrame
    selected_columns = df[['x', 'y', 'yearly_avg']].rename(columns={
    'yearly_avg': 'Si'
    })
    # normalize the score
    selected_columns['Si'] = selected_columns['Si']/ selected_columns['Si'].max()
    # Write the selected columns to a CSV file with custom headers
    selected_columns.to_csv(result_fp, header=True, mode='w', sep='\t', index=False)


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

    parser.add_argument('sci_raw_fn', type=str, help='Raw science filename, File must be in the data/science_data directory. E.g., "Fracking_25_scival.csv"') 
    parser.add_argument('result_fn', type=str, help='Processed science filename E.g., "Fracking_25_scival.ophs"') 

    args = parser.parse_args()
    
    main(args.sci_raw_fn, args.result_fn)