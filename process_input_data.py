import pandas as pd
from pyproj import Proj
import matplotlib.pyplot as plt
import os
import argparse

def main(sci_raw_data_fp, result_fp, score_col_num=3):
    """
    Processes a CSV file containing science data (from Molly/ George) and projects coordinates using a Lambert Conformal Conic projection.
    The output is a file with lists of vertex locations and their corresponding scores.
    Saves the output as a tab-delimited file with projected coordinates and normalized values from the third column (column position = 2) as scores.

    Example:
        python process_input_data_ver2.py Fracking_25_scival.csv Fracking_25_scival.opc

    Args:
        sci_raw_data_fp (str): Absolute path to the input science file.
        result_fp (str): Abosulte path to the output file (e.g., "Fracking_25_scival.opc") containing the vertices and their scores.
        score_col_num (int): Column number of the data file which is to be considered for score calculation

    Notes:
        Input files must be located in the 'data/science_data' directory.
    """
    # Check if input file exists
    if not os.path.exists(sci_raw_data_fp):
        raise FileNotFoundError(f"Input file '{sci_raw_data_fp}' not found.")

    # Read the CSV file into a pandas DataFrame
    df = pd.read_csv(sci_raw_data_fp)

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
    print("df['x'].min(): ", df['x'].min())
    print("df['y'].min(): ", df['y'].min())
    df['x'] = df['x'] - df['x'].min()
    df['y'] = df['y'] - df['y'].min()

    # Get column name by position
    column_name = df.columns[score_col_num]
    print(f'Column {column_name} has been selected for scores.')    
    #  Select columns and normalize the values as score Si
    selected_columns = df[['x', 'y', column_name]].rename(columns={
    column_name: 'Si'
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
    parser = argparse.ArgumentParser(description='Processes CSV science data file from scientists and generates a projected output file with vertices and their scores.')
    parser.add_argument('sci_raw_data_fp', type=str, help='Absolute path to the input science file."') 
    parser.add_argument('result_fp', type=str, help='Absolute path to the output file (e.g., "Fracking_25_scival.opc") containing the vertices and their scores.') 
    args = parser.parse_args()
    
    main(args.sci_raw_data_fp, args.result_fp)