# Orienteering Problem with Charging Stations

## Scripts

`analysis_cases.py`

Script to do batch processing on a several different cases (vary scenario parameters, science data, etc).

`main_driver.py`

Script that orchestrates the MIP planner flow. The input science data file is expected to be in the 'data' folder.
The results are written as a pickle file in the 'results' folder.

`make_hotels.py`

This script strategically places H hotels using K-Means clustering based on input data of vertex locations and scores.

`process_input_data_ver2.py`

Processes a CSV file containing science data (from scientists) and projects coordinates using a Lambert Conformal Conic projection.
Saves the output as a tab-delimited file with projected coordinates and normalized yearly averages as scores.

`scip/uav_charging_ver5.py`

MIP planner script with SCIP solver for planning UAV flight paths in the presence of multiple charging stations.


`print_results_ver3.py`

Script to load results, print them, and plot (Cartesian or geographic coordinates) if necessary.
Plotting on geographic coordinates with a river map in the background requires the `river_data/river_map.pkl` file or corresponding shape files, which have not been synchornized with Github since the file size is large.

