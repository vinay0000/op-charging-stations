import pandas as pd
from pyproj import Proj
import matplotlib.pyplot as plt
import os
import argparse

def main(sci_raw_fn, result_fn):
    """
    Processes a CSV file containing science data (from Molly/ George) and projects coordinates using a Lambert Conformal Conic projection.
    Saves the output as a tab-delimited file with projected coordinates and normalized yearly averages as scores.

    Example:
        python process_input_data_ver2.py Fracking_25_scival.csv Fracking_25_scival.opc

    Args:
        sci_raw_fn (str): Name of the input CSV file (e.g., "Fracking_25_scival.csv"). File must be in the data/science_data directory.
        result_fn (str): Name of the output file (e.g., "Fracking_25_scival.opc").

    Notes:
        Input files must be located in the 'data/science_data' directory.
    """
    # Define filepaths
    current_dir = os.path.dirname(os.path.abspath(__file__))
    sci_file_dir= 'data/science_data'

    sci_raw_fp = os.path.join(sci_file_dir, sci_raw_fn)
    result_fp = os.path.join(sci_file_dir, result_fn)

    # Check if input file exists
    if not os.path.exists(sci_raw_fp):
        raise FileNotFoundError(f"Input file '{sci_raw_fp}' not found.")

    # Read the CSV file into a pandas DataFrame
    df = pd.read_csv(os.path.join(current_dir, sci_raw_fp))

    # Define Lambert Conformal Conic projection parameters (Texas)
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
        try:
            x, y = lambert_conformal_conic(row['long'], row['lat'])
            return pd.Series({'x': x, 'y': y})
        except Exception as e:
            raise ValueError(f"Error projecting coordinates: {e}")

    # Apply the projection function to each row
    df[['x', 'y']] = df.apply(project_coords, axis=1)

    # Adjust the coordinates so that the minimum x and y values start at zero.    
    df['x'] = df['x'] - df['x'].min()
    df['y'] = df['y'] - df['y'].min()
    
    #  Select columns and normalize the 'yearly_avg' field as score Si
    selected_columns = df[['x', 'y', 'yearly_avg']].rename(columns={
    'yearly_avg': 'Si'
    })
    selected_columns['Si'] = selected_columns['Si']/ selected_columns['Si'].max()
    
    # Save the output in a tab-delimited format
    selected_columns.to_csv(result_fp, header=True, mode='w', sep='\t', index=False)


    #### Plotting the points weighted by the scores ####
    to_plot = False
    if to_plot is True:
        plt.figure(figsize=(10, 8))
        plt.scatter(df['x'], df['y'], c='blue', marker='o', edgecolor='k', label='Projected Points')
        plt.xlabel('X coordinate')
        plt.ylabel('Y coordinate')
        plt.title('Projected Coordinates (Lambert Conformal Conic)')
        plt.legend()
        plt.grid(True)
        plt.show()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Processes CSV science data file from scientists and generates a projected output file with targets and scores.')
    parser.add_argument('sci_raw_fn', type=str, help='Raw science filename (e.g., "Fracking_25_scival.csv"). File must be in the data/science_data directory."') 
    parser.add_argument('result_fn', type=str, help='Output filename (e.g., "Fracking_25_scival.opc")."') 
    args = parser.parse_args()
    
    main(args.sci_raw_fn, args.result_fn)